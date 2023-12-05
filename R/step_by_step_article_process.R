
#BiocManager::install("scater")
library(scater)
library(R.utils)


##对所需数据文件（绝对路径）进行解压缩：
gunzip("GitHub/LearningRNA-seq/extdata/GSE61533_HTSEQ_count_results.xls.gz",
       remove=FALSE, 
       overwrite=TRUE)

##install.packages("gdata")

library(gdata)
library(readxl)

all.counts <- read_excel("GitHub/LearningRNA-seq/extdata/GSE61533_HTSEQ_count_results.xls",
                         col_names=TRUE)
all.counts <- as.data.frame(all.counts)
rownames(all.counts) <- all.counts[,1]

all.counts <- all.counts[,-1]

class(all.counts)
#利用函数读取矩阵：
library(SingleCellExperiment)


sce <- SingleCellExperiment(list(counts = as.matrix(all.counts)))


rownames(sce)
dim(sce)

is.spike <- grepl("^ERCC", rownames(sce))
is.mito <- grepl("^mt-", rownames(sce))

summary(is.spike)
Summary(is.mito)
            
sce <- perCellQCMetrics(sce)

#sce <- calculateQCMetrics(sce, feature_controls=list(ERCC=is.spike, Mt=is.mito))

head(colnames(pData(sce)))
     
     
library(scran)
isSpike(sce) <- "ERCC"   

par(mfrow=c(1,2))
hist(sce$total_counts/1e6, xlab="Library sizes (millions)", main="",
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features, xlab="Number of expressed genes", main="",
     breaks=20, col="grey80", ylab="Number of cells")


libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)










# analyzing sc-RNA-seq data containing UMI counts -------------------------

####analyzing single-cell RNA-seq data containing UMI counts

readFormat <- function(infile) { 
  # First column is empty.
  metadata <- read.delim(infile, stringsAsFactors=FALSE, header=FALSE, nrow=10)[,-1] 
  rownames(metadata) <- metadata[,1]
  metadata <- metadata[,-1]
  metadata <- as.data.frame(t(metadata))
  
  # First column after row names is some useless filler.
  counts <- read.delim(infile, stringsAsFactors=FALSE, 
                       header=FALSE, row.names=1, skip=11)[,-1] 
  counts <- as.matrix(counts)
  return(list(metadata=metadata, counts=counts))
}

##读取文件来自网络可采用网络调取文件的方式解决路径（由于下载时间以及下载rnames不唯一，采用直接下载数据到本地）
##从网络得到需要的mRNA,spike,mito数据(已下载文件，使用绝对路径)
library(hexView)
endo.data <- readFormat("GitHub/LearningRNA-seq/extdata/expression_mRNA_17-Aug-2014.txt")

spike.data <- readFormat("GitHub/LearningRNA-seq/extdata/expression_spikes_17-Aug-2014.txt")

mito.data <- readFormat("GitHub/LearningRNA-seq/extdata/expression_mito_17-Aug-2014.txt")


##需要对线粒体数据进行格式的处理，使具有相同data格式
m <- match(endo.data$metadata$cell_id, mito.data$metadata$cell_id)

mito.data$metadata <- mito.data$metadata[m,]

mito.data$counts <- mito.data$counts[,m]

##我们将对应于单个基因的所有行计数相加。
raw.names <- sub("_loc[0-9]+$", "", rownames(endo.data$counts))

new.counts <- rowsum(endo.data$counts, group=raw.names, reorder=FALSE)

endo.data$counts <- new.counts

##由于原函数newSCESet已无法使用，所以使用函数SingleCellExperiment 处理，将这些计数矩阵合并成一个单独矩阵

library(SingleCellExperiment)

all.counts <- rbind(endo.data$counts, mito.data$counts, spike.data$counts)

sce <- SingleCellExperiment(list(counts=all.counts), colData=endo.data$metadata)

dim(sce)

##对矩阵数据进行基于基因的注释来标注每一行：
# nrows <- c(nrow(endo.data$counts), nrow(mito.data$counts), nrow(spike.data$counts))
nrows <- c(nrow(endo.data$counts), nrow(mito.data$counts), nrow(spike.data$counts))

is.spike <- rep(c(FALSE, FALSE, TRUE), nrows)

is.mito <- rep(c(FALSE, TRUE, FALSE), nrows)

# isSpike(sce, "Spike") <- is.spike

#初问题报错 未能找到isspike功能：采取换代码：
#******is.spike <- grepl("^ERCC", rownames(sce))
#******sce3 <- splitAltExps(sce1, ifelse(is.spike, "ERCC", "gene"))


## Adding Ensembl IDs.
library(org.Mm.eg.db)

ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")

rowData(sce)$ENSEMBL <- ensembl

sce


## Quality control on the cells
##在作者原始数据已经过滤低质量细胞的基础上，再次进行质控
library(scater)

# sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito)) 

sce <- perCellQCMetrics(sce,subsets=list(Spike=is.spike, Mt=is.mito))


##绘图进行检查：
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))

hist(sce$total/1e3, xlab="Library sizes (thousands)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")

hist(sce$detected, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")

hist(sce$subsets_Mt_percent, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")

hist(sce$subsets_Spike_percent, xlab="ERCC proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")

## 移除异常值（包括文库大小，特征表达的数量，尖端峰值）
libsize.drop <- isOutlier(sce$total, nmads=3, type="lower", log=TRUE)

feature.drop <- isOutlier(sce$detected, nmads=3, type="lower", log=TRUE)

spike.drop <- isOutlier(sce$subsets_Spike_percent, nmads=3, type="higher")

##从sce数据框或矩阵中选择那些在libsize.drop、feature.drop和spike.drop中对应的
##布尔值为FALSE的样本和特征，即保留那些需要保留的样本和特征。

##sce <- sce[,!(libsize.drop | feature.drop | spike.drop)] 

data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), 
           BySpike=sum(spike.drop), Remaining=ncol(sce))


library(scran)

mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))

assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)

table(assignments$phase)

plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)


ave.counts <- rowMeans(counts(sce))
keep <- rowMeans(counts(sce)) >= 0.2

fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))

#plotQC函数已被弃用：plotQC(sce, type = "highest-expression", n=50) + fontsize

plotHighestExprs(sce)

ave.counts <- calcAverage(sce, use_size_factors=FALSE)

hist(log10(ave.counts), breaks=100, main="", col="grey",
     xlab=expression(Log[10]~"average count"))


clusters <- quickCluster(sce, min.mean=0.1, method="igraph")

sce <- computeSumFactors(sce, cluster=clusters, min.mean=0.1)

summary(sizeFactors(sce))

plot(sizeFactors(sce), sce$total_counts/1e3, log="xy",
     ylab="Library size (thousands)", xlab="Size factor")








