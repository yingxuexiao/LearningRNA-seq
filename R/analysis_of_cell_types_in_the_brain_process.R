readFormat <- function(infile) {
  # First column is empty.
  metadata <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, nrow=10)[,-1]
  rownames(metadata) <- metadata[,1]
  metadata <- metadata[,-1]
  metadata <- as.data.frame(t(metadata))
  # First column after row names is some useless filler.
  counts <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, row.names=1, skip=11)[,-1]
  counts <- as.matrix(counts)
  return(list(metadata=metadata, counts=counts))
}


##读取文件来自网络可采用网络调取文件的方式解决路径（由于下载时间以及下载rnames不唯一，采用直接下载数据到本地）
library(BiocFileCache)

bfc <- BiocFileCache("raw_data", ask = FALSE)


base.url <- file.path("https://storage.googleapis.com",
                      "linnarsson-lab-www-blobs/blobs/cortex")

mRNA.path <- bfcrpath(bfc, file.path(base.url, 
                                     "expression_mRNA_17-Aug-2014.txt"))
mito.path <- bfcrpath(bfc, file.path(base.url, 
                                     "expression_mito_17-Aug-2014.txt"))
spike.path <- bfcrpath(bfc, file.path(base.url, 
                                      "expression_spikes_17-Aug-2014.txt"))


##从网络得到需要的mRNA,spike,mito数据(已下载文件，使用绝对路径)

library(hexView)
endo.data <- readFormat("GitHub/LearningRNA-seq/extdata/expression_mRNA_17-Aug-2014.txt")

spike.data <- readFormat("GitHub/LearningRNA-seq/extdata/expression_spikes_17-Aug-2014.txt")

mito.data <- readFormat("GitHub/LearningRNA-seq/extdata/expression_mito_17-Aug-2014.txt")


##需要对线粒体数据进行格式的处理，使具有相同data格式
m <- match(endo.data$metadata$cell_id, mito.data$metadata$cell_id)

mito.data$metadata <- mito.data$metadata[m,]

mito.data$counts <- mito.data$counts[,m]


##由于原函数newSCESet已无法使用，所以使用函数SingleCellExperiment 处理，将这些计数矩阵合并成一个单独矩阵
library(SingleCellExperiment)

all.counts <- rbind(endo.data$counts, mito.data$counts, spike.data$counts)

sce <- SingleCellExperiment(list(counts=all.counts), colData=endo.data$metadata)

dim(sce)

##对矩阵数据进行基于基因的注释来标注每一行：
nrows <- c(nrow(endo.data$counts), nrow(mito.data$counts), nrow(spike.data$counts))

is.spike <- rep(c(FALSE, FALSE, TRUE), nrows)

is.mito <- rep(c(FALSE, TRUE, FALSE), nrows)

isSpike(sce, "Spike") <- is.spike


#### QC on the cell 
##在作者原始数据已经过滤低质量细胞的基础上，再次进行质控：
library(scater)
#sce <- calculateQCMetrics(sce, feature_controls=list(Spike=is.spike, Mt=is.mito))

sce <- perCellQCMetrics(sce,subsets=list(Spike=is.spike, Mt=is.mito))

isSpike(sce) <- "Spike"

##绘图进行检查：

par(mfrow=c(1,2))

hist(sce$total/1e3, xlab="Library sizes (thousands)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")

hist(sce$detected, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")

hist(sce$subsets_Mt_percent, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")

hist(sce$subsets_Spike_percent, xlab="ERCC proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")


##删除异常表达的值：
libsize.drop <- isOutlier(sce$total, nmads=3, type="lower", log=TRUE)

feature.drop <- isOutlier(sce$detected, nmads=3, type="lower", log=TRUE)

mito.drop <- isOutlier(sce$subsets_Mt_percent, nmads=3, type="higher")

spike.drop <- isOutlier(sce$subsets_Spike_percent, nmads=3, type="higher")


##然后通过组合所有指标的过滤器来移除低质量的细胞
#出问题
#sce <- sce[,!(libsize.drop | feature.drop | spike.drop | mito.drop)]

data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           ByMito=sum(mito.drop), BySpike=sum(spike.drop), Remaining=ncol(sce))


#### 细胞周期分类cell cycle classification

#BiocManager::install("org.Mm.eg.db")

library(org.Mm.eg.db)

anno <- select(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")

ensembl <- anno$ENSEMBL[match(rownames(sce), anno$SYMBOL)]

assignments <- cyclone(sce, mm.pairs, gene.names=ensembl)

plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)

library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))

assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)
table(assignments$phase)