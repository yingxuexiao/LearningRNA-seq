##利用DESeq2考虑协变量影响(这是知道变量是什么的情况使用)
#使用同样的例子进行说明（例子来源为绝对路径，由于安装包compGenomRData找不到）

counts_file <- ("C://Users/Administrator/Documents/GitHub/LearningRNA-seq/extdata/SRP021193/SRP021193.raw_counts.tsv")
colData_file <- ("C://Users/Administrator/Documents/GitHub/LearningRNA-seq/extdata/SRP021193/SRP021193.colData.tsv")

counts <- read.table(counts_file)
colData <- read.table(colData_file, header = T, sep = '\t',
                      stringsAsFactors = TRUE)
#通过计算TPM来观察样本聚类
library(pheatmap)

#找到基因长度标准化值
geneLengths <-counts$width

rpk <- apply( subset(counts, select = c(-width)), 2,
              function(x) x/(geneLengths/1000))

#利用RPK值对样本大小进行标准化
tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)
#（sort（）用于排序，这里使用前100降序排列）
selectedGenes <- names(sort(apply(tpm, 1, var),
                            decreasing = T)[1:100])

pheatmap(tpm[selectedGenes,],
         scale = 'row',
         annotation_col = colData,
         show_rownames = FALSE)

#使用DESeq2来考虑、分析由其上可知的变量“library selection”：
library(DESeq2)

# 将'width'从矩阵移除
countData <- as.matrix(subset(counts, select = c(-width)))

# 设置一个DESeqDataSet对象
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~LibrarySelection+group)

#运行差异基因表达的结果。
dds <- DESeq(dds)

DEresults <- results(dds, contrast = c('group', 'CASE', 'CTRL'))




###使用RUVSeq来估计协变量的方法（这种针对不清楚潜在变量是什么的情况）
counts_file <- ("C://Users/Administrator/Documents/GitHub/LearningRNA-seq/extdata/SRP049988/SRP049988.raw_counts.tsv")

colData_file <- ("C://Users/Administrator/Documents/GitHub/LearningRNA-seq/extdata/SRP049988/SRP049988.colData.tsv")

counts <- read.table(counts_file)
colData <- read.table(colData_file, header = T,
                      sep = '\t', stringsAsFactors = TRUE)


colData$source_name <- ifelse(colData$group == 'CASE',
                              'EHF_overexpression', 'Empty_Vector')
#利用TPM值绘制热图
geneLengths <-counts$width
rpk <- apply( subset(counts, select = c(-width)), 2,
              function(x) x/(geneLengths/1000))

tpm <- apply(rpk, 2, function(x) x / sum(as.numeric(x)) * 10^6)

selectedGenes <- names(sort(apply(tpm, 1, var),
                            decreasing = T)[1:100])

pheatmap(tpm[selectedGenes,],
         scale = 'row',
         annotation_col = colData,
         cutree_cols = 2,
         show_rownames = FALSE)

#这一步的热图并不能帮助我们查看出确切的变量，所以需要使用RUVSeq去预测可能的协变量
library(EDASeq)
# remove'width'column from counts
countData <- as.matrix(subset(counts, select = c(-width)))

# create a seqExpressionSet object using EDASe qpackage
set <- newSeqExpressionSet(counts = countData,
                           phenoData = colData)

#对原始数据进行RLE 图像和主成分分析并且对分组颜色标记
par(mfrow = c(1,2))

plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=as.numeric(colData$group))

plotPCA(set, col = as.numeric(colData$group), adj = 0.5,
        ylim = c(-0.7, 0.5), xlim = c(-0.5, 0.5))

#对TPM 标准化的数据进行RLE 和主成分分析
par(mfrow = c(1,2))

plotRLE(tpm, outline=FALSE, ylim=c(-4, 4), col=as.numeric(colData$group))

plotPCA(tpm, col=as.numeric(colData$group), adj = 0.5,
        ylim = c(-0.3, 1), xlim = c(-0.5, 0.5))

#利用管家基因作为参考移除不想要的变量
library(RUVSeq)

#选取HK_genes作为管家基因（采取绝对路径）
HK_genes <- read.table(file = "C://Users/Administrator/Documents/GitHub/LearningRNA-seq/extdata/HK_genes.txt",
                       header = FALSE)

#查看管家基因和可用的基因之间的交集
house_keeping_genes <- intersect(rownames(set), HK_genes$V1)

#使用house_keeping_genes作为经验基因集，在RUVg中尝试不同的K值观察主成分分析图像
par(mfrow = c(2, 2))
for(k in 1:4) {
  set_g <- RUVg(x = set, cIdx = house_keeping_genes, k = k)
  plotPCA(set_g, col=as.numeric(colData$group), cex = 0.9, adj = 0.5,
          main = paste0('with RUVg,k=',k),
          ylim = c(-1, 1), xlim = c(-1, 1), )
}

#上述过程的图像可以帮助筛选出K=1的plot，接着再次用RUVg，此时K=1做进一步的诊断：
set_g <- RUVg(x = set, cIdx = house_keeping_genes, k = 1)

##比较在有无RUVg对数据的标准化和移除变量的情况下，RLE、PCA图像的区别
#比较RLE图像
par(mfrow = c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4),
        col=as.numeric(colData$group), main = 'without RUVg')

plotRLE(set_g, outline=FALSE, ylim=c(-4, 4),
        col=as.numeric(colData$group), main = 'with RUVg')

#比较PCA图像
par(mfrow = c(1,2))

plotPCA(set, col=as.numeric(colData$group), adj = 0.5,
        main = 'without RUVg',
        ylim = c(-1, 0.5), xlim = c(-0.5, 0.5))

plotPCA(set_g, col=as.numeric(colData$group), adj = 0.5,
        main = 'with RUVg',
        ylim = c(-1, 0.5), xlim = c(-0.5, 0.5))




##除了使用RUVg还可以使用RUVs:(此方法不需要使用管家基因)
differences <- makeGroups(colData$group)

par(mfrow = c(2, 2))

#用全部基因的信息中寻找两组不同的不想要的变量
for(k in 1:4) {
  set_s <- RUVs(set, unique(rownames(set)), 
                k=k, differences) 
  plotPCA(set_s, col=as.numeric(colData$group),
          cex = 0.9, adj = 0.5,
          main = paste0('with RUVs,k=',k),
          ylim = c(-1, 1), xlim = c(-0.6, 0.6))
}

#根据PCA图确定K=2时，明显分离，于是选择k=2做后续的RUVs函数分析
set_s <- RUVs(set, unique(rownames(set)), k=2, differences)

#比较有无RUVs处理分析后的RLE图像
par(mfrow = c(1,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4),
        col=as.numeric(colData$group),
        main = 'without RUVs')

plotRLE(set_s, outline=FALSE, ylim=c(-4, 4),
        col=as.numeric(colData$group),
        main = 'with RUVs')

#PCA图像区别
par(mfrow = c(1,2))

plotPCA(set, col=as.numeric(colData$group),
        main = 'without RUVs', adj = 0.5,
        ylim = c(-0.75, 0.75), xlim = c(-0.75, 0.75))

plotPCA(set_s, col=as.numeric(colData$group),
        main = 'with RUVs', adj = 0.5,
        ylim = c(-0.75, 0.75), xlim = c(-0.75, 0.75))


#上述使用RUVg和RUVs两种策略处理，现比较原始数据、和分别使用两种策略的PCA图像：
par(mfrow = c(1,3))

plotPCA(countData, col=as.numeric(colData$group),
        main = 'without RUV-rawcounts', adj = 0.5,
        ylim = c(-0.75, 0.75), xlim = c(-0.75, 0.75))

plotPCA(set_g, col=as.numeric(colData$group),
        main = 'with RUVg', adj = 0.5,
        ylim = c(-0.75, 0.75), xlim = c(-0.75, 0.75))

plotPCA(set_s, col=as.numeric(colData$group),
        main = 'with RUVs', adj = 0.5,
        ylim = c(-0.75, 0.75), xlim = c(-0.75, 0.75))

#可以根据plot图像得出使用RUVs效果更好，因此接着使用RUVs处理后的原始数据做热图分析：
library(EDASeq)
library(pheatmap)

# extract normalized counts that are cleared from unwanted variation using RUVs
normCountData <- normCounts(set_s)
selectedGenes <- names(sort(apply(normCountData, 1, var),
                            decreasing = TRUE))[1:500]

pheatmap(normCountData[selectedGenes,],
         annotation_col = colData,
         show_rownames = FALSE,
         cutree_cols = 2,
         scale = 'row')


#再次运行DESeq2使用已经计算后的协变量。
#利用RUVs处理后的数据和DESeq2，可以再次进行差异表达分析
library(DESeq2)

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ group)

#过滤低计数基因
dds <-dds[rowSums(DESeq2::counts(dds)) > 10]

# insert the covariates W1 and W2 computed using RUVs into DESeqDataSet object
colData(dds) <- cbind(colData(dds),
                      pData(set_s)[rownames(colData(dds)),
                                   grep('W_[0-9]',
                                        colnames(pData(set_s)))])

#
design(dds) <- ~ W_1 + W_2 + group
# repeat the analysis
dds <- DESeq(dds)

# extract deseq results
res <- results(dds, contrast = c('group', 'CASE', 'CTRL'))
res <-res[order(res$padj),]




