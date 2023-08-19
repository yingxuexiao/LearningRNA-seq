
#quality check on sequencing reads
library(ShortRead)
library(BiocParallel)

#预准备确定文件具体位置
setwd("C:/Users/11/Rwork/scRNA-seq/bigdata/fastq/E-MTAB-1147/")
dirPath <- "C:/Users/11/Rwork/scRNA-seq/bigdata/fastq/E-MTAB-1147" 


#A  利用rqc()抓取对应的fastq file 
library(BiocManager)
BiocManager::install("Rqc")
library(Rqc)

qcRes=rqc(path=dirPath,pattern=".fastq.gz",openBrowser = FALSE)

#B  Per base sequence quality boxplot.
rqcCycleQualityBoxPlot(qcRes)

#C  Percentage of nucleotide bases per position
rqcCycleBaseCallsLinePlot(qcRes)

#D  Check the percent of different duplication levels in FASTQ files
rqcReadFrequencyPlot(qcRes)





#filter and trim raeds
#读取一个fastq文件，并过滤每个质量分数低于20的读数。

library(ShortRead)
#A 读取fastq文件
fastqFile <-"C:/Users/11/Rwork/scRNA-seq/bigdata/fastq/ERR127302_1_subset.fastq.gz" 

#B  read fastq file
fq <- readFastq(fastqFile)

#C  get quality scores per base as a matrix
qPerBase = as(quality(fq), "matrix")

#D  get number of bases per read that have quality score below 20
qcount = rowSums( qPerBase <= 20)

#E  Number of reads where all Phred scores >= 20
fq[qcount == 0]

#F  write out fastq file with only reads where all quality
#   scores per base are above 20
writeFastq(fq[qcount == 0],
           paste(fastqFile, "Qfiltered", sep="_"))

#G  fastq文件很大，常用的方法是通过片段化逐条读取
## set up streaming with block size 1000
f <- FastqStreamer(fastqFile,readerBlockSize=1000)

#H   we set up a while loop to call yield() function to
# go through the file
while(length(fq <- yield(f))) {
  qPerBase = as(quality(fq), "matrix")
  qcount = rowSums( qPerBase <= 20)
  writeFastq(fq[qcount == 0],
             paste(fastqFile, "Qfiltered", sep="_"),
             mode="a")
}




#Mapping/aligning reads to the genome举例
library(QuasR)

#A  copy example data to current working directory
file.copy(system.file(package="QuasR", "extdata"),
          ".", recursive=TRUE)

#B  genome file in fasta format
genomeFile <- "extdata/hg19sub.fa"

#C  text file containing sample names and fastq file paths
sampleFile <- "extdata/samples_chip_single.txt"

#D create alignments
proj <- qAlign(sampleFile, genomeFile)


#当fastqcr需要用MAC|Linux时，需要采用以下算法帮助fastqcr使用。
.check_if_unix <<- function() {
  return(NULL)
}
assignInNamespace(".check_if_unix", .check_if_unix, ns = "fastqcr")

wsl.exe /path_to_wsl_fastqc_install/fastqc

fastqcr::fastqc("fastq_path", fastqc.path = "wsl.exe /path_to_wsl_fastqc_install/fastqc")
