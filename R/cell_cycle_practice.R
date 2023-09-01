#单细胞测序 周期分析
library(scRNAseq)
BiocManager::install("scRNAseq")
sce.416b <- LunSpikeInData(which="416b") 
sce.416b$block <- factor(sce.416b$block)
