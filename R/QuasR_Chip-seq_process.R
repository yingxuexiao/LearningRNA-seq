


###确定当前工作目录路径

getwd()
###[1] "C:/Users/Administrator/Documents"

###Align reads using the qAlign function
##此步骤，必须确定reads已经进行预处理，样本文件已经准备好，确保有充足的磁盘内存。
##处理的文档来自于QuasR包自带的数据集（使用时路径尽量完整）
library(QuasR)

sampleFile <- "GitHub/LearningRNA-seq/extdata/samples_chip_single.txt"

auxFile <- "GitHub/LearningRNA-seq/extdata/auxiliaries.txt"  #准备的辅助文件

genomeFile <- "GitHub/LearningRNA-seq/extdata/hg19sub.fa"

proj1 <- qAlign(sampleFile, 
                genome = genomeFile, 
                auxiliaryFile = auxFile)
proj1

##比对后每个输入序列都应得到一个bam文件，同时auxfile也应该得到一个bam文件，
##标准格式：list.files("extdata", pattern = ".bam$")
list.files("C:/Users/Administrator/Documents/GitHub/LearningRNA-seq/extdata", 
           pattern = ".bam$")


##bam文件名由序列文件的基本名称和添加的随机字符串组成。
##随机后缀确保新生成的比对文件不会覆盖现有的比对文件
list.files("GitHub/LearningRNA-seq/extdata", 
           pattern = "^chip_1_1_")[1:3]


###Create a quality control report
##QuasR能生成质量控制报告（qQCReport内部使用ShortRead、FastQC质控工具）
qQCReport(proj1, 
          pdfFilename = "GitHub/LearningRNA-seq/extdata/qc_report.pdf")
##creating QC plots：reads质量，碱基含量等生成PDF报告。


###Alignment statistics
##读取bam文件的 mapped reads并生成表格。
alignmentStats(proj1)


###Export genome wig file from alignments
##根据比对结果输出（用来绘制基因组上的图像轨迹）
qExportWig(proj1, 
           binsize = 100L, 
           scaling = TRUE, 
           collapseBySample = TRUE)

###Count alignments using qCount
##1）在Chip-seq中使用qcount需要加载包GenomicFeatures、
library(GenomicFeatures)


##2）首先从带有基因注释的.gtf文件中创建一个TxDb文件
annotFile <- "GitHub/LearningRNA-seq/extdata/hg19sub_annotation.gtf"

chrLen <- scanFaIndex(genomeFile)
#scanFaIndex()依赖于包Rsamtools.

chrominfo <- data.frame(chrom = as.character(seqnames(chrLen)),
                        length = width(chrLen),
                        is_circular = rep(FALSE, length(chrLen)))

txdb <- makeTxDbFromGFF(file = annotFile, format = "gtf",
                        chrominfo = chrominfo,
                        dataSource = "Ensembl",
                        organism = "Homo sapiens")

txdb

promReg <- promoters(txdb, upstream = 1000, downstream = 500,
                     columns = c("gene_id","tx_id"))

gnId <- vapply(mcols(promReg)$gene_id,
               FUN = paste, FUN.VALUE = "",
               collapse = ",")

promRegSel <- promReg[ match(unique(gnId), gnId) ]

names(promRegSel) <- unique(gnId)

promRegSel

cnt <- qCount(proj1, promRegSel)
cnt

##输出为原始比对数据，结果未经任何标准化。
##可用Gviz（）对数据进行可视化，观察到CpG岛启动子上的H3K4me3
#gr1 <- import("Sample1.wig.gz")

gr1 <- "Sample1.wig.gz"
               

gr2 <- "Sample2.wig.gz"

library(Gviz)

axisTrack <- GenomeAxisTrack()

dTrack1 <- DataTrack(range = gr1, 
                     name = "Sample 1", 
                     type = "h")

dTrack2 <- DataTrack(range = gr2, 
                     name = "Sample 2",
                     type = "h")

txTrack <- GeneRegionTrack(txdb, 
                           name = "Transcripts", 
                           showId = TRUE)

plotTracks(list(axisTrack, dTrack1, dTrack2, txTrack),
           chromosome = "chr3", extend.left = 1000)
