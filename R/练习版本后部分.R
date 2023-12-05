
# 基于 spike-in coverage 的标准化（可供选择的方法） --------------------------------------



counts <- read.table("GitHub/LearningRNA-seq/extdata/GSE29087_L139_expression_tab.txt.gz",
                     colClasses=c(list("character",NULL, NULL, NULL, NULL, NULL, NULL), rep("integer", 96)), skip=6, sep='\t', row.names=1)


library(SingleCellExperiment)

is.spike <- grep("SPIKE", rownames(counts)) 

sce.islam <- SingleCellExperiment(list(counts=as.matrix(counts)))

sce.islam <- splitAltExps(sce.islam, ifelse(is.spike, "ERCC", "gene"))

isSpike(sce.islam, "spike") <- is.spike

dim(sce.islam)

#利用函数去除低质量细胞，识别每种细胞中的异常值

library(scater)

sce.islam <- perCellQCMetrics(sce.islam)

sce.islam$grouping <- rep(c("mESC", "MEF", "Neg"), c(48, 44, 4))

ibsize.drop <- isOutlier(sce.islam$total, nmads=3, type="lower", 
                         log=TRUE, batch=sce.islam$grouping)

feature.drop <- isOutlier(sce.islam$detected, nmads=3, type="lower", 
                          log=TRUE, batch=sce.islam$grouping)

spike.drop <- isOutlier(sce.islam$grouping, nmads=3, type="higher", 
                        batch=sce.islam$grouping)

sce.islam <- sce.islam[,!(libsize.drop | feature.drop | 
                            spike.drop | sce.islam$grouping=="Neg")]

data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           BySpike=sum(spike.drop), Remaining=ncol(sce.islam))





sce <- SingleCellExperiment(list(counts = as.matrix(counts)))

#sce <- newSCESet(countData=counts)

sce$grouping <- rep(c("mESC", "MEF", "Neg"), c(48, 44, 4))

sce1 <- sce[,sce$grouping!="Neg"] # Removing negative control wells.

#sce <- calculateQCMetrics(sce, feature_controls=list(spike=grep("SPIKE", rownames(counts))))

sce2 <- perCellQCMetrics(sce1,subsets=list(spike=grep("SPIKE", rownames(counts))))

#isSpike(sce) <- "spike"
is.spike <- grepl("^ERCC", rownames(sce1))

sce3 <- splitAltExps(sce, ifelse(is.spike, "ERCC", "gene"))

library(scuttle)

sce4 <- computeSpikeFactors(sce3,assay.type = "SPIKE")


