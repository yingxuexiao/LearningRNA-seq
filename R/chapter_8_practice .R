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

#A.create a vector of gene lengths

geneLengths <- as.vector(subset(counts,select=c(width)))


rpkm <- apply(X =subset(counts,select=c(-width)),
              MARGIN=2,
              FUN=function(x){10^9*x/geneLengths/sum(as.numeric(x))})
rpkm

colSums(rpkm)

########Computing TPM

rpk <- apply(subset(counts,select=c(-width)),2,
             function(x)x/(geneLengths/1000))

tpm <- apply(rpk,2,function(x)x/sum(as.numeric(x))*10^6)

tpm
 
colSums(tpm)



###clustering 
#A.数据方差分析
v <- apply(tpm,1,var)
v

#B. selecting genes
selectedgenes <- names(v[order(v,decreasing=T)][1:100])
selectedgenes

#C. produce a heatmap
install.packages("pheatmap")
library(pheatmap)
pheatmap(tpm[selectedgenes,],scale='row',show_rownames = FALSE)

coldata <- read.table(coldata_file,header=T,sep='\t',stringsAsFactors = TRUE)


######PCA
library(stats)
library(ggplot2)

#A. transpose the matrix
M <- t(tpm[selectedgenes,])

M <- log2(M+1)
#B. compute PCA
pcaResults <- prcomp(M)

autoplot(pcaResults,data=coldata,colour=('group'))

summary(pcaResults)

####correlation plots
library(stats)
correlationMatrix <- cor(tpm)

install.packages("corrplot")
library(corrplot)
corrplot(correlationMatrix,order='hclust',
         addrect=2,addCoef.col='white',
         number.cex=0.7)

library(pheatmap)

pheatmap(correlationMatrix,
         annotation_col = coldata,
         cutree_cols=2)

##### differential expression analysis practicing

countData <- as.matrix(subset(counts,select=c(-width)))

colData <- read.table(coldata_file,header=T,sep='\t',
                      stringsAsFactors = TRUE)

designFormula <- "~group"

library()
