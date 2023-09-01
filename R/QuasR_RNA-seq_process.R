getwd()
###[1] "C:/Users/Administrator/Documents"

library(QuasR)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           

###Spliced alignment of RNA-seq reads

##单细胞测序中reads比对需要使用Rhisat2包
#BiocManager::install("Rhisat2")

sampleFile <- "GitHub/LearningRNA-seq/extdata/samples_rna_paired.txt"
genomeFile <- "GitHub/LearningRNA-seq/extdata/hg19sub.fa"

proj2 <- qAlign(sampleFile, genome = genomeFile,
                splicedAlignment = TRUE, aligner = "Rhisat2")
proj2

##在进行比对过程中使用splicedAlignment = TRUE，
##会导致交叉外显子进行对齐，将和参考基因组比较有较大的缺失
##所以使用splicedAlignment = FALSE

proj2unspl <- qAlign(sampleFile, genome = genomeFile,
                     splicedAlignment = FALSE)

proj2unspl 

#生成比对PDF报表
qQCReport(proj2unspl,"GitHub/LearningRNA-seq/extdata/qc_report1.pdf")

#根据比对结果输出（用来绘制基因组上的图像轨迹）
qExportWig(proj2unspl, 
           binsize = 100L, 
           scaling = TRUE, 
           collapseBySample = TRUE)

##可以通过数据比较可明显发现两者区别
alignmentStats(proj2)
alignmentStats(proj2unspl)


###对基因及外显子表达进行定量分析
##extracting gene regions from TxDb,比对计数、
geneLevels <- qCount(proj2, txdb, reportLevel = "gene")
geneLevels

##reportLevel=“..”可以用于控制注释外显子区域是否单独进行定量
exonLevels <- qCount(proj2, txdb, reportLevel = "exon")
exonLevels

head(geneLevels)
head(exonLevels)


###标准化，用RPKM expression values进行计算
##由qcount得到比对的数量，但是往往还需要进一步标准化：
##区域的长度，文库的大小。这里使用RPKM值进行标准化
geneRPKM <- t(t(geneLevels[,-1] / geneLevels[,1] * 1000)
              / colSums(geneLevels[,-1]) * 1e6)

geneRPKM


###Analysis of alternative splicing: Quantification of exon-exon junctions
##reportLevel=“junction”进行外显子—外显子连接质控
exonJunctions <- qCount(proj2, NULL, reportLevel = "junction")
exonJunctions

knownIntrons <- unlist(intronsByTranscript(txdb))

isKnown <- overlapsAny(exonJunctions, knownIntrons, type = "equal")

table(isKnown)

tapply(rowSums(as.matrix(mcols(exonJunctions))),
       isKnown, summary)

exonBodyLevels <- qCount(proj2, txdb, reportLevel = "exon",
                         includeSpliced = FALSE)
exonBodyLevels
summary(exonLevels - exonBodyLevels)




