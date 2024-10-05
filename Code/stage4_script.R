######################################
#### HackBio internship - stage 4 #### 
######################################

# Biomarkers analysis #### 

## load necessary libraries 

library("TCGAbiolinks")
library(SummarizedExperiment)
library(edgeR)
library(gplots)
library(ggplot2)
library(biomaRt)
library(ggrepel)

# project information 
getProjectSummary("TCGA-LGG")
getProjectSummary("TCGA-GBM")

# download LGG dataset 
tcga_lgg <- GDCquery(project = "TCGA-LGG",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification")
GDCdownload(tcga_lgg) 
lgg_data <- GDCprepare(tcga_lgg) 

# download GBM dataset
tcga_gbm <- GDCquery(project = "TCGA-GBM",
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification")
GDCdownload(tcga_gbm) 
gbm_data <- GDCprepare(tcga_gbm) 

# get raw counts from lgg_data
lgg_rawdata <- assays(lgg_data) 
dim(lgg_rawdata$unstranded) # 534 lgg samples

# get raw counts from gbm_data
gbm_rawdata <- assays(gbm_data) 
dim(gbm_rawdata$unstranded) # 176 gbm samples

# prepare metadata for lgg
lgg_data$barcode
lgg_data$paper_IDH.status
table(lgg_data$paper_IDH.status) # 94 wt and 419 mut

meta_lgg <- data.frame("barcode" = lgg_data$barcode,
                       "IDH status" = lgg_data$paper_IDH.status)
View(meta_lgg)
table(is.na(meta_lgg)) # 21 NAs -> eliminate samples with no IDH status
meta_lgg <- na.omit(meta_lgg)

# prepare metadata for gbm 
gbm_data$barcode
gbm_data$paper_IDH.status
table(gbm_data$paper_IDH.status) # 143 wt and 11 mut 

meta_gbm <- data.frame("barcode" = gbm_data$barcode,
                       "IDH status" = gbm_data$paper_IDH.status)
View(meta_gbm)
table(is.na(meta_gbm)) # 22 NAs -> eliminate samples with no IDH status
meta_gbm <- na.omit(meta_gbm)

# final dataset for lgg (513 samples)
lgg <- lgg_rawdata$unstranded
samples_lgg <- meta_lgg$barcode
lgg <- lgg[, colnames(lgg) %in% samples_lgg]
all(colnames(lgg) %in% samples_lgg) # TRUE

# final dataset for gbm (154 samples)
gbm <- gbm_rawdata$unstranded
samples_gbm <- meta_gbm$barcode
gbm <- gbm[, colnames(gbm) %in% samples_gbm]
all(colnames(gbm) %in% samples_gbm)

# merge the two datasets (both counts and metadata)
glioma <- cbind(lgg, gbm)
meta_glioma <- rbind(meta_lgg, meta_gbm)

# export the final dfs
write.csv(glioma, "glioma_rawcounts.csv", row.names = TRUE)
write.csv(meta_glioma, "glioma_metadata.csv", row.names = FALSE)

# clean and preprocess the data
table(is.na(glioma)) # no NAs
glioma <- TCGAanalyze_Normalization(tabDF = glioma, geneInfo = geneInfoHT, method = "geneLength")
glioma <- TCGAanalyze_Normalization(tabDF = glioma, geneInfo = geneInfoHT, method = "gcContent")
glioma <- TCGAanalyze_Filtering(tabDF = glioma,
                                method = "quantile",
                                qnt.cut = 0.25)

## differential expression analysis 
group.wt <- meta_glioma$barcode[meta_glioma$IDH.status == "WT"]
group.mut <- meta_glioma$barcode[meta_glioma$IDH.status == "Mutant"]

DEA <- TCGAanalyze_DEA(mat1 = glioma[ , group.wt], 
                       mat2 = glioma[ , group.mut],
                       Cond1type = "WT",
                       Cond2type = "Mutant",
                       pipeline = "edgeR")
DEA.Level <- 
  TCGAanalyze_LevelTab(DEA, "WT", "Mutant",
                       glioma[ , group.wt],
                       glioma[ , group.mut])

# annotation
gene_names <- rownames(DEA.Level)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
annot <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = gene_names,
  mart = mart)
gene_with_symbols <- merge(DEA.Level, annot, 
                           by.x = "mRNA", 
                           by.y = "ensembl_gene_id", 
                           all.x = TRUE)

# volcano plot
genes_to_label <- subset(gene_with_symbols, abs(logFC) > 5 & FDR < 0.05)

ggplot(DEA.Level, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = ifelse(FDR < 0.05 & abs(logFC) > 1,
                                ifelse(logFC > 1, "Upregulated", "Downregulated"), "Not significant")), size = 2) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") + 
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") + 
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
  theme_minimal() +
  labs(x = "Log2 Fold Change", y = "-Log10 FDR", title = "Volcano Plot of WT vs Mut", color = "Gene Regulation") +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14), 
        plot.title = element_text(hjust = 0.5, size = 16),
        legend.position = "top") +
  ggrepel::geom_label_repel(data = genes_to_label, aes(label = hgnc_symbol))


# heatmap of DEGs
# top differentially expressed genes 
DEGs <- subset(DEA.Level, abs(logFC) > 3 & FDR < 0.01) # 905 DEGs   
heat.DEGs <- glioma[rownames(DEGs), ]

# color code
ccodes <- ifelse(colnames(heat.DEGs) %in% group.wt, "red", "blue") # wt in red, mut in blue

# heatmap
red_green_palette <- colorRampPalette(c("green", "red"))(10)
heatmap.2(x = as.matrix(heat.DEGs),
          col = red_green_palette,
          Rowv = F, Colv = T,
          scale = 'row',
          sepcolor = 'black',
          trace = "none",
          key = TRUE,
          dendrogram = "col",
          cexRow = 0.9, cexCol = 0.7,
          main = "Heatmap of WT vs Mut",
          na.color = 'black',
          ColSideColors = ccodes,
          margins = c(11,10))
legend("left",                       
       legend = c("WT", "Mut"),     
       col = c("red", "blue"),           
       lty = 1,                          
       lwd = 4,                         
       cex = 0.7,
       xpd = TRUE,
       inset = c(-0.1, 0))


# Machine Learning analysis #### 
# load necessary packages
library(caret)
library(DALEX)
library(pROC)
library(iml)
set.seed(1234)

# Load the main and meta data
glioma.data <- read.csv(file = "glioma_rawcounts.csv", header = TRUE)
glioma.meta <- read.csv(file = "glioma_metadata.csv", header = TRUE)
dim(glioma.data)
rownames(glioma.data) <- glioma.data$X

glioma.data$X <- NULL

colnames(glioma.data) <- gsub("\\.", "-", colnames(glioma.data))
glioma.data <- log10(glioma.data + 1)

glioma.data <- data.frame(t(glioma.data))


SDs = apply(glioma.data, 2, sd)
# This line calculates the standard deviation for each column (each gene) in lgg.data and stores the results in the vector SDs. 
# A higher standard deviation indicates greater variability in the gene expression levels.

topPredicts = order(SDs, decreasing = T)[1:3000]
#After calculating the standard deviation, this lines ensures that the top 2000 genes with the highest SD are selected, and in descending order.


glioma.data = glioma.data[, topPredicts]
#The trans_data is now reduced to only include the top 3000 genes that have the highest variability, 
# which are often more informative for subsequent analyses like classification

rownames(glioma.meta) <- glioma.meta$barcode

glioma.meta$barcode <- NULL

merge.data <- merge(glioma.data, glioma.meta, by = "row.names")

rownames(merge.data) <- merge.data$Row.names

merge.data$Row.names <- NULL

dim(merge.data)

anyNA(merge.data)
sum(is.na(merge.data))



# To remove near zero variation

zero <- preProcess(merge.data, method = "nzv", uniqueCut = 15)
merge.data <- predict(zero, merge.data)


# center. Centering helps to stabilize numerical computations 
#and can be particularly important for algorithms sensitive to the scale of the data (like PCA, KNN, etc.).
center <- preProcess(merge.data, method = "center")
merge.data <- predict(center, merge.data)



# to remove highly correlated values.
#Reduce Multicollinearity: High correlations between features can lead to issues in model training, 
#such as inflated standard errors and difficulties in interpreting coefficients. 
#Removing redundant features helps create a more stable model.
# Improve Model Performance: By reducing the number of features, 
# you can improve the efficiency and performance of certain machine learning algorithms, especially those sensitive to multicollinearity.

corr <-preProcess(merge.data, method = "corr", cutoff = 0.5)
merge.data <- predict(corr, merge.data)

dim(merge.data)

glioma.train <- createDataPartition(y = merge.data$IDH.status, p = 0.7)[[1]]

train.glioma <- merge.data[glioma.train,]
test.glioma <- merge.data[-glioma.train,]


ctrl.glioma <- trainControl(method = "cv", number = 10)


knn.glioma <- train(IDH.status~., #the levels of classification
                 data = train.glioma, #the training dataset
                 method = "knn", # the knn method
                 trControl = ctrl.glioma, #the training control
                 tuneGrid = expand.grid(k = seq(1, 50, by = 2)))

#the best K is:
knn.glioma$bestTune

#predict
trainPred.glioma <- predict(knn.glioma, newdata = train.glioma)
testPred.glioma <- predict(knn.glioma, newdata = test.glioma)

trainPred.glioma <- as.factor(trainPred.glioma)
train.glioma$IDH.status <- as.factor(train.glioma$IDH.status)
testPred.glioma <- as.factor(testPred.glioma)
test.glioma$IDH.status <- as.factor(test.glioma$IDH.status)

levels(trainPred.glioma) <- levels(train.glioma$IDH.status)
levels(testPred.glioma) <- levels(test.glioma$IDH.status)

# Interpretation
# confusion matrix
confusionMatrix(trainPred.glioma, train.glioma$IDH.status)
confusionMatrix(testPred.glioma, test.glioma$IDH.status)


explainer.glioma <- explain(knn.glioma, 
                         label = "knn",
                         data = train.glioma,
                         y = as.numeric(train.glioma$IDH.status))


importance.glioma <- feature_importance(explainer.glioma, n_sample = 30)
top10 <- head(importance.glioma$variable, 11)
last10 <- tail(importance.glioma$variable, 12)

# annotation
top10_genes <- head(importance.glioma$variable, 11)
last10_genes <- tail(importance.glioma$variable, 12)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
annot_top <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = top10_genes,
  mart = mart)
annot_top

annot_last <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol'),
  filters = 'ensembl_gene_id',
  values = last10_genes,
  mart = mart)
annot_last
