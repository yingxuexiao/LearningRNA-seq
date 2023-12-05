
# R:  4.3.1




# analyzing sc-RNA-seq data containing UMI counts -------------------------

####analyzing single-cell RNA-seq data containing UMI counts


# 数据准备 --------------------------------------------------------------------

#设置需要使用的函数形式
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
is.spike <- grepl("^ERCC", rownames(sce))

sce1 <- splitAltExps(sce, ifelse(is.spike, "ERCC", "gene"))

## Adding Ensembl IDs.
library(org.Mm.eg.db)

ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce1), keytype="SYMBOL", column="ENSEMBL")

rowData(sce1)$ENSEMBL <- ensembl

sce1


# 细胞质控 --------------------------------------------------------------------

## Quality control on the cells
##在作者原始数据已经过滤低质量细胞的基础上，再次进行质控

#注：下游需要保持is.spike  is.mito和sce的长度相同，去除两者后端57个数据：
x <- is.spike[1:19839]
y <- is.mito[1:19839]


library(scater)
#sce2 <- perCellQCMetrics(sce,subsets=list(Spike=is.spike,Mt=is.mito))

sce2 <- perCellQCMetrics(sce1,subsets=list(Spike=x, Mt=y))

##绘图进行检查：
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))

hist(sce2$total/1e3, xlab="Library sizes (thousands)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")

hist(sce2$detected, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")

hist(sce2$subsets_Mt_percent, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")

hist(sce2$altexps_ERCC_percent, xlab="ERCC proportion (%)",
     ylab="Number of cells", breaks=20, main="", col="grey80")




## 移除异常值（包括文库大小，特征表达的数量，尖端峰值）
libsize.drop <- isOutlier(sce2$total, nmads=3, type="lower", log=TRUE)

feature.drop <- isOutlier(sce2$detected, nmads=3, type="lower", log=TRUE)

spike.drop <- isOutlier(sce2$altexps_ERCC_percent, nmads=3, type="higher")

##从sce数据框或矩阵中选择那些在libsize.drop、feature.drop和spike.drop中对应的
##布尔值为FALSE的样本和特征，即保留那些需要保留的样本和特征。

sce3 <- sce[,!(libsize.drop | feature.drop | spike.drop)] 

data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), 
           BySpike=sum(spike.drop), Remaining=ncol(sce3))




# 细胞周期分类 ------------------------------------------------------------------
#cell cycle classification

#BiocManager::install("org.Mm.eg.db")
library(scran)

mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))

assignments <- cyclone(sce1, mm.pairs, gene.names=rowData(sce1)$ENSEMBL)

table(assignments$phase)

plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)



# Examining gene-level metrics --------------------------------------------

fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))

#*******画不出图像来plotHighestExprs(sce)+fontsize

#利用细胞平均计数来量化基因丰度：
library(scater)

calculateAverage(sce)

ave.counts <- calculateAverage(sce)

hist(log10(ave.counts), breaks=100, main="", col="grey",
     xlab=expression(Log[10]~"average count"))

##在此基础上进行过滤，去除0值无表达基因

rowData(sce)$ave.count <- ave.counts

to.keep <- ave.counts > 0

sce <- sce[to.keep,]

summary(to.keep)

# 标准化细胞特异性偏差 --------------------------------------------------------------
##聚类细胞进行标准化，利用去卷积的方法，通过减少相同类细胞的差异基因的数量来提高标准化准确性
##然后进行缩放以确保不同簇中细胞的大小因子具有可比性。
library(scran)

clusters <- quickCluster(sce, min.mean=0.1, method="igraph")

sce4 <- computeSumFactors(sce, cluster=clusters, min.mean=0.1)

summary(sizeFactors(sce4))

plot(sizeFactors(sce4), sce2$total/1e3, log="xy",
     ylab="Library size (thousands)", xlab="Size factor")

##同样可以计算内参的大小因子
library(scuttle)
sce5 <- computeSpikeFactors(sce1,"ERCC")

#sce6 <- normalize(sce5)   

sce6 <- normalizeCounts(sce5)
sce6 <- normalizeCounts(sce)

# 消除技术噪音 ------------------------------------------------------------------

##我们通过拟合峰值转录本的均值方差趋势来模拟技术噪声

sce7 <- logNormCounts(sce)

means <- rowMeans(logcounts(sce7))

vars <- rowVars(logcounts(sce7))

var.fit <- fitTrendVar(means, vars)

#代码过时 var.out <- decomposeVar(sce7, var.fit)

var.out <- modelGeneVar(sce7)

plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
     ylab="Variance of log-expression")

points(var.fit$mean, var.fit$var, col="red", pch=16)

o <- order(var.out$mean)

lines(var.out$mean[o], var.out$tech[o], col="red", lwd=2)

curve(var.fit$trend(x), col="dodgerblue", add=TRUE, lwd=2)

##高变基因的筛选
hvg.out <- var.out[which(var.out$FDR <= 0.05 & var.out$bio >= 0.5),]

hvg.out <- hvg.out[order(hvg.out$bio, decreasing=TRUE),]

nrow(hvg.out)

write.table(file="brain_hvg.tsv", hvg.out, sep="\t", quote=FALSE, col.names=NA)

head(hvg.out)

chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]

#plotExpression(sce7, rownames(var.out)[chosen.genes], 
              # alpha=0.05, jitter="jitter") + fontsize

plotExpression(sce7, rownames(var.out)[chosen.genes], 
               point_size = 0.05, jitter="jitter") + fontsize

set.seed(100)

var.cor <- correlatePairs(sce7, subset.row=rownames(hvg.out))

write.table(file="brain_cor.tsv", var.cor, sep="\t", quote=FALSE, row.names=FALSE)

head(var.cor)


##install.packages("irlba")
library(irlba)
sce8 <- denoisePCA(sce7, technical=var.fit$trend)


sce8 <- denoisePCA(sce, technical=var.fit$trend)

ncol(reducedDim(sce8, "PCA"))

library(limma)

adj.exprs <- exprs(sce8)
adj.exprs <- removeBatchEffect(adj.exprs, batch=sce8$sex)
norm_exprs(sce8) <- adj.exprs

chosen <- unique(c(var.cor$gene1, var.cor$gene2))

top.hvg <- rownames(hvg.out)[1]

# TSNE在s4文件中不能找到：tsne1 <- plotTSNE(sce, exprs_values="norm_exprs", colour_by=top.hvg,
                  #perplexity=10, rand_seed=100, feature_set=chosen) + fontsize
#tsne2 <- plotTSNE(sce, exprs_values="norm_exprs", colour_by="Mog",
                  #perplexity=10, rand_seed=100, feature_set=chosen) + fontsize
#multiplot(tsne1, tsne2, cols=2)

sce <- runTSNE(sce8, dimred="PCA", perplexity=50)

tsne1 <- plotTSNE(sce, colour_by="Neurod6") + fontsize
tsne2 <- plotTSNE(sce, colour_by="Mog") + fontsize

par(mfrow=c(1,2))

pca1 <- plotReducedDim(sce8, dimred="PCA", colour_by="Neurod6") + fontsize
pca2 <- plotReducedDim(sce8, dimred="PCA", colour_by="Mog") + fontsize

#pca1 <- plotPCA(sce8, exprs_values="norm_exprs", colour_by=top.hvg) + fontsize

#pca2 <- plotPCA(sce8, exprs_values="norm_exprs", colour_by="Mog") + fontsize

#multiplot(pca1, pca2, cols=2)



# 聚类 ----------------------------------------------------------------------

chosen.exprs <- norm_exprs(sce)[chosen,]

my.dist <- dist(t(chosen.exprs))

my.tree <- hclust(my.dist, method="ward.D2")
# 用树来聚类
library(dynamicTreeCut)

my.clusters <- unname(cutreeDynamic(my.tree, distM=as.matrix(my.dist), verbose=0))

heat.vals <- chosen.exprs - rowMeans(chosen.exprs)

clust.col <- rainbow(max(my.clusters))

library(gplots)

heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=0.3,
          ColSideColors=clust.col[my.clusters], Colv=as.dendrogram(my.tree))


# 检测标记基因： -----------------------------------------------------------------

cluster <- factor(my.clusters)

de.design <- model.matrix(~0 + cluster + sce$sex)

head(colnames(de.design))

library(edgeR)

y <- convertTo(sce, type="edgeR")

y <- estimateDisp(y, de.design)

fit <- glmFit(y, de.design)

summary(y$tagwise.dispersion)

##每个基因都在选择的集群和数据集中的每个其他集群之间进行DE测试

result.logFC <- result.PValue <- list()

chosen.clust <- which(levels(cluster)=="1") # character, as 'cluster' is a factor.

for (clust in seq_len(nlevels(cluster))) {
  if (clust==chosen.clust) { next }
  contrast <- numeric(ncol(de.design))
  contrast[chosen.clust] <- 1
  contrast[clust] <- -1
  res <- glmLRT(fit, contrast=contrast)
  con.name <- paste0('vs.', levels(cluster)[clust])
  result.logFC[[con.name]] <- res$table$logFC
  result.PValue[[con.name]] <- res$table$PValue
}


##为了从每个比较的前10个基因中构建一个标记集，我们将筛选标记：

collected.ranks <- lapply(result.PValue, rank, ties="first")

min.rank <- do.call(pmin, collected.ranks)

marker.set <- data.frame(Top=min.rank, Gene=rownames(y),
                         logFC=do.call(cbind, result.logFC), stringsAsFactors=FALSE)

marker.set <- marker.set[order(marker.set$Top),]

head(marker.set, 10)

write.table(marker.set, file="brain_marker_1.tsv", sep="\t", quote=FALSE, col.names=NA)

top.markers <- marker.set$Gene[marker.set$Top <= 10]

top.exprs <- norm_exprs(sce)[top.markers,,drop=FALSE]

heat.vals <- top.exprs - rowMeans(top.exprs)

heatmap.2(heat.vals, col=bluered, symbreak=TRUE, trace='none', cexRow=0.6,
          ColSideColors=clust.col[my.clusters], Colv=as.dendrogram(my.tree), dendrogram='none')

legend("bottomleft", col=clust.col, legend=sort(unique(my.clusters)), pch=16)




