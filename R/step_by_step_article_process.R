
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


library(hexView)
endo.data <- readFormat("GitHub/LearningRNA-seq/extdata/expression_mRNA_17-Aug-2014.txt")

spike.data <- readFormat("GitHub/LearningRNA-seq/extdata/expression_spikes_17-Aug-2014.txt")

mito.data <- readFormat("GitHub/LearningRNA-seq/extdata/expression_mito_17-Aug-2014.txt")


m <- match(endo.data$metadata$cell_id, mito.data$metadata$cell_id)
mito.data$metadata <- mito.data$metadata[m,]
mito.data$counts <- mito.data$counts[,m]

raw.names <- sub("_loc[0-9]+$", "", rownames(endo.data$counts))
new.counts <- rowsum(endo.data$counts, group=raw.names, reorder=FALSE)
endo.data$counts <- new.counts

library(SingleCellExperiment)
all.counts <- rbind(endo.data$counts, mito.data$counts, spike.data$counts)
sce <- SingleCellExperiment(list(counts=all.counts), colData=endo.data$metadata)
dim(sce)


nrows <- c(nrow(endo.data$counts), nrow(mito.data$counts), nrow(spike.data$counts))
is.spike <- rep(c(FALSE, FALSE, TRUE), nrows)
is.mito <- rep(c(FALSE, TRUE, FALSE), nrows)
isSpike(sce, "Spike") <- is.spike

# Adding Ensembl IDs.
library(org.Mm.eg.db)
ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
rowData(sce)$ENSEMBL <- ensembl

sce

library(scater)
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito)) 

sce <- perCellQCMetrics(sce,subsets=list(Spike=is.spike, Mt=is.mito))


par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))

hist(sce$total/1e3, xlab="Library sizes (thousands)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")

hist(sce$detected, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")

hist(sce$subsets_Mt_percent, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")

hist(sce$subsets_Spike_percent, xlab="ERCC proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")


libsize.drop <- isOutlier(sce$total, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$detected, nmads=3, type="lower", log=TRUE)
spike.drop <- isOutlier(sce$subsets_Spike_percent, nmads=3, type="higher")

##sce <- sce[,!(libsize.drop | feature.drop | spike.drop)]

data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), 
           BySpike=sum(spike.drop), Remaining=ncol(sce))


library(scran)

mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))

#assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)

#table(assignments$phase)

#plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)


ave.counts <- rowMeans(counts(sce))
keep <- rowMeans(counts(sce)) >= 0.2

fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))

#plotQC函数已被弃用：plotQC(sce, type = "highest-expression", n=50) + fontsize

plotHighestExprs(sce)

ave.counts <- calcAverage(sce, use_size_factors=FALSE)

hist(log10(ave.counts), breaks=100, main="", col="grey",
     xlab=expression(Log[10]~"average count"))











