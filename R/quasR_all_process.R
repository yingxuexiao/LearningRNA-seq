####setwd("C:/Users/Administrator/Documents/GitHub/LearningRNA-seq")

###QuasR: short for Quantify and annotate short reads in R
###QuasR: starting from the raw sequence reads, over pre-processing and alignment, up to quantification.

###installation
BiocManager::install("QuasR")
library(QuasR)

###需要对应的一些包：
suppressPackageStartupMessages({
  library(QuasR)
  library(BSgenome)
  library(Rsamtools)
  library(rtracklayer)
  library(GenomicFeatures)
  library(Gviz)
})
###quasR所用的例子需要安装下面的包：
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

library(BSgenome.Hsapiens.UCSC.hg19)

##确定文件路径
sampleFile <- "LearningRNA-seq/extdata/samples_chip_single.txt"
genomeFile <- "LearningRNA-seq/extdata/hg19sub.fa"

##进行比对
proj <- qAlign(sampleFile, genomeFile)

proj
#对比对结果生成PDF报告
qQCReport(proj, "extdata/qc_report.pdf")

## creating QC plots生产指控报告
library(rtracklayer)
library(GenomicFeatures)

annotFile <- "extdata/hg19sub_annotation.gtf"

txStart <- import.gff(annotFile, format = "gtf", feature.type = "start_codon")

promReg <- promoters(txStart, upstream = 500, downstream = 500)

names(promReg) <- mcols(promReg)$transcript_name

promCounts <- qCount(proj, query = promReg)


###选择参考基因组使用BSgenome
library(BSgenome)
available.genomes()

genomeName <- "BSgenome.Hsapiens.UCSC.hg19"

###选择合适的基因组。需要避免冗长

###测序数据的预处理：利用preprocessReads()处理，移除3末端，移除接头，过滤低质量reads
td <- tempdir()

infiles <- system.file(package = "QuasR", "extdata",
                       c("rna_1_1.fq.bz2","rna_2_1.fq.bz2"))

outfiles <- file.path(td, basename(infiles))

res <- preprocessReads(filename = infiles,
                       outputFilename = outfiles,
                       truncateEndBases = 3,
                       Lpattern = "AAAAAAAAAA",
                       minLength = 14, 
                       nBases = 2)
unlink(outfiles)

###双末端测序采用严格的过滤，即使一条符合要求，那么也会两条链对应的部分都移除。
td <- tempdir()

infiles1 <- system.file(package = "QuasR", "extdata", "rna_1_1.fq.bz2")

infiles2 <- system.file(package = "QuasR", "extdata", "rna_1_2.fq.bz2")

outfiles1 <- file.path(td, basename(infiles1))

outfiles2 <- file.path(td, basename(infiles2))

res <- preprocessReads(filename = infiles1,
                       filenameMate = infiles2,
                       outputFilename = outfiles1,
                       outputFilenameMate = outfiles2,
                       nBases = 0)









