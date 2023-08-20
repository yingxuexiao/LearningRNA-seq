#download data with R
#A.use github to download SRP029880 data

#B  读取数据

counts_file <- system.file("extdata/rna-seq/SRP029880.raw_counts.tsv",
                           package = "compGenomRData")

#counts_file <- ("C://Users/Administrator/R.work/SRP029880.raw_counts.tsv")

coldata_file <- system.file("extdata/rna-seq/SRP029880.colData.tsv",
                            package = "compGenomRData")

#coldata_file <- ("C://Users/Administrator/R.work/SRP029880.colData.tsv")

#C.将数据转换为矩阵

counts <- as.matrix(read.table(counts_file,header=T,sep='\t'))

###########Computing  CPM

summary(counts[,1:3])

cpm <- apply(subset(counts,select=c(-width)),2,
             function(x) x/sum(as.numeric(x))*10^6)
                    
colSums(cpm)

###########Computing RPKM

#A.创建一个基因长度的矩阵

geneLengths <- as.vector(subset(counts,select=c(width)))
#计算RPKM

rpkm <- apply(X =subset(counts,select=c(-width)),
              MARGIN=2,
              FUN=function(x){10^9*x/geneLengths/sum(as.numeric(x))})
rpkm

colSums(rpkm)

########Computing TPM

#查找基因长度归一化值
rpk <- apply(subset(counts,select=c(-width)),2,
             function(x)x/(geneLengths/1000))

#按照基因大小进行RPK标准化
tpm <- apply(rpk,2,function(x)x/sum(as.numeric(x))*10^6)

tpm
 
colSums(tpm)


###clustering 
#A.数据方差分析
v <- apply(tpm,1,var)
v

#B.利用标准方差整理结果及降序命令挑选前100的基因
selectedgenes <- names(v[order(v,decreasing=T)][1:100])
selectedgenes

#C. 将结果绘制成热图
install.packages("pheatmap")
library(pheatmap)
pheatmap(tpm[selectedgenes,],scale='row',show_rownames = FALSE)

#D. 将结果coldata_file绘制成热图
coldata <- read.table(coldata_file,header=T,sep='\t',stringsAsFactors = TRUE)
pheatmap(tpm[selectedgenes,],scale='row',show_rownames = FALSE)

######PCA主成分分析
library(stats)
library(ggplot2)

#A. 将得到的数据矩阵进行转置（转置矩阵）
M <- t(tpm[selectedgenes,])

#将数据进行对数转换
M <- log2(M+1)

#B. 计算主成分分析
pcaResults <- prcomp(M)

#利用ggplot2中的autoplot进行主成分分析可视化
autoplot(pcaResults,data=coldata,colour=('group'))

summary(pcaResults)

##相关的Plot
library(stats)
correlationMatrix <- cor(tpm)

install.packages("corrplot")
library(corrplot)
corrplot(correlationMatrix,order='hclust',
         addrect=2,addCoef.col='white',
         number.cex=0.7)

library(pheatmap)

##根据聚类相似度将聚类分成两个
pheatmap(correlationMatrix,
         annotation_col = coldata,
         cutree_cols=2)

##### 差异性分析：

countData <- as.matrix(subset(counts,select=c(-width)))

colData <- read.table(coldata_file,header=T,sep='\t',
                      stringsAsFactors = TRUE)

designFormula <- "~group"

#进行DESeq 进行基因差异性表达分析
library(DESeq2)
library(stats)

#获取一个有关于以上获得的count data及coldata的DESeq数据集
dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData = colData,
                              design=as.formula(designFormula))

#打印dds查看其内容
print(dds)

#将基因中的read进行筛选，将read数少于1的移除
dds <- dds[rowSums(DESeq2::counts(dds))>1,]

#利用DESeq处理数据
dds <- DESeq(dds)

#计算差异性，其中“CTRL”作为对照组
DEresults=results(dds,contrast = c("group",'CASE','CTRL'))

#通过提高p值来整理结果
DEresults <- DEresults[order(DEresults$pvalue),]
print(DEresults)

###诊断性图像：诊断性检测帮助提升数据的质量、实验系统
#利用MA plot进行诊断检测，用于观察标准化结果是否良好：大部分散点集中于0轴
library(DESeq2)
DESeq2::plotMA(object=dds,ylim=c(-5,5))

##p值分布
#利用图像观察P值分布，观察到峰值和一致的p值分布
library(ggplot2)
ggplot(data=as.data.frame(DEresults),aes(x=pvalue))+
  geom_histogram(bins=100)

##PCA plot
#我们需要调取标准化的数据进行绘图，并且对不同的实验使用不同颜色观察聚类效果
library(DESeq2)
countsNormalized <- DESeq2::counts(dds,normalized=TRUE)

#获取前500的可变基因
selectedGenes <- names(sort(apply(countsNormalized,1,var),
                            decreasing=TRUE)[1:500])

plotPCA(countsNormalized[selectedGenes,],
        col=as.numeric(colData$group),adj=0.5,
        xlim=c(-0.5,0.5),ylim=c(-0.5,0.6))

#除了使用plotPCA(),还可以选择使用DESeq2::rlog及DESeq2::plotPCA()来绘出PCA的图像。
rld <- rlog(dds)
DESeq2::plotPCA(rld,ntop=500,intgroup='group')+
  ylim(-50,50)+theme_bw()

#除了PCA之外，还可以进行RLE帮助找出还需要进行标准化的数据，快速对原始或者
#标准化的数据进行诊断，查看是否需要进一步的处理

library(DESeq)
par(mfrow=c(1,2))
plotRLE(countData,outline=FALSE,ylim=c(-4,4),
        col=as.numeric(colData$group),
        main='Raw Counts')
plotRLE(DESeq2::counts(colData$group),
        mian='Normalized Counts')





