# Integrating Machine Learning and Differential Expression Analysis for Biomarker Discovery in Glioma 🧬💻🤖

 **Authors:** Manal Agdada (@Manal), Ojiaku Confidence Chisom (@Areta), Abdulrahman Walid Elbagoury (@Willeau), Rahma Nabil Sallam (@rahmanabil2002), Pascal Onaho (@PascalOnaho), Hagar Haitham Elazab (@HBONH33), Mercy Francis (@Mercylee), Ariyo Adesokan (@Adesokan_ariyo1)

**R script:** https://github.com/manal-agdada/Stage-4---HackBio-internship/blob/main/Code/stage4_script.R

**Figures:** https://github.com/manal-agdada/Stage-4---HackBio-internship/tree/main/Figures

## Introduction 
Gliomas are the most common malignant brain tumors, classified into low-grade gliomas (LGG) and glioblastoma (GBM) [1]. Isocitrate dehydrogenase (IDH) mutations are key prognostic markers, as IDH-mutant gliomas have better prognoses and distinct molecular characteristics compared to IDH-wildtype tumors [2]. IDH-mutant gliomas are associated with a hypermethylation phenotype and tend to have more favorable outcomes [2]. In contrast, IDH-wildtype gliomas exhibit more aggressive behavior [2]. Therefore, methylation patterns in gliomas can serve as a valuable tool for molecular and clinical classification of gliomas [2].

## Dataset Description and preprocessing steps
For this analysis, we used TCGA LGG (534 samples) and GBM (176 samples) datasets. RNA-seq data were processed in RStudio using the ‘TCGAbiolinks’ package. After filtering samples by IDH status(wt/mut), and normalizing expression data, differential expression analysis was performed using the edgeR protocol, identifying 9899 differentially expressed genes.

![heatmap](https://github.com/user-attachments/assets/a51c0122-917a-44dd-945e-eb8dabf76f47)
![caption1](https://github.com/user-attachments/assets/2d5293b1-2bf2-47de-afc1-cd80f48f2c7d)

![volcanoplot](https://github.com/user-attachments/assets/9a4e6128-cbdf-460e-b080-e0e4bc8840fe)
![caption2](https://github.com/user-attachments/assets/dc06f511-48ce-4a57-83e1-dfbeccab70ed)

## Methodology and Results

## Conclusion

## References
1. Louis DN, Perry A, Reifenberger G, von Deimling A, Figarella-Branger D, Cavenee WK, Ohgaki H, Wiestler OD, Kleihues P, Ellison DW. The 2016 World Health Organization Classification of Tumors of the Central Nervous System: a summary. Acta Neuropathol. 2016 Jun;131(6):803-20. 
2. Ceccarelli M, Barthel FP, Malta TM, Sabedot TS, Salama SR, Murray BA, Morozova O, Newton Y, Radenbaugh A, Pagnotta SM, Anjum S, Wang J, Manyam G, Zoppoli P, Ling S, Rao AA, Grifford M, Cherniack AD, Zhang H, Poisson L, Carlotti CG Jr, Tirapelli DP, Rao A, Mikkelsen T, Lau CC, Yung WK, Rabadan R, Huse J, Brat DJ, Lehman NL, Barnholtz-Sloan JS, Zheng S, Hess K, Rao G, Meyerson M, Beroukhim R, Cooper L, Akbani R, Wrensch M, Haussler D, Aldape KD, Laird PW, Gutmann DH; TCGA Research Network; Noushmehr H, Iavarone A, Verhaak RG. Molecular Profiling Reveals Biologically Discrete Subsets and Pathways of Progression in Diffuse Glioma. Cell. 2016 Jan 28;164(3):550-63. 
