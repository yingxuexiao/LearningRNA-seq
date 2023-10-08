##quality check on sequencing reads
#预准备确定文件具体位置
dirPath <- "C:/Users/Administrator/Documents/GitHub/LearningRNA-seq/data/E-MTAB-1147" 

#1)  利用rqc()抓取对应的fastq文件并且得到质量检测报告html版

library(Rqc)
qcRes <- rqc(path=dirPath,
             pattern=".fastq.gz",
             openBrowser = FALSE)

#得到html报告：'C:\Users\Administrator\AppData\Local\Temp\
#RtmpM5vw6e/rqc_report.html' has been created.

# Showing input files
knitr::kable(perFileInformation(qcRes))

#2) 利用对应的函数得到报告中对应的数据
# 获取单碱基测序的质量报告，帮助考虑是否Trim得分低于20的reads，
#也能帮助检测是否存在接头污染
rqcCycleQualityBoxPlot(qcRes)

#3)获取每个位置核苷酸碱基的百分比的图像
rqcCycleBaseCallsLinePlot(qcRes)

#4)  检查FASTQ文件中不同重复层次的百分比，一般的reads只复制一次
rqcReadFrequencyPlot(qcRes)

#5)  获取GC含量图像
rqcCycleGCPlot(qcRes)

#6) Heatmap of top represented reads
rqcFileHeatmap(qcRes[[1]])

#7) Average Quality
rqcReadQualityPlot(qcRes)

#8) Read Length Distribution
rqcReadWidthPlot(qcRes)





###利用ShortRead函数进行数据处理：读取一个fastq文件，并过滤每个质量分数低于20的读数。
library(ShortRead)

#1) 读取fastq文件
fastqFile <-system.file(package = "ShortRead",
                        "extdata/E-MTAB-1147",
                        "ERR127302_1_subset.fastq.gz") 

#2) 读取fastq file
fq <- readFastq(fastqFile)
fq

#3)  利用矩阵得到每个碱基的质量得分
qPerBase = as(quality(fq), "matrix")
qPerBase
dim(qPerBase)

#4)  获取碱基质量得分<20的数量
qcount = rowSums(qPerBase <= 20)

qcount
#5)  得到菲尔德质量评分>=20的reads 数量
fq[qcount == 0]

#6)  读写出每个碱基中所有read得分高于20的fastq 文件
writeFastq(fq[qcount == 0],
           paste(fastqFile, "Qfiltered", sep="_"))

#7)  fastq文件很大，常用的方法是通过片段化逐条读取
## set up streaming with block size 1000，帮助我们顺利读取片段化的reads
f <- FastqStreamer(fastqFile,readerBlockSize=1000)
f

while(length(fq <- yield(f))) {
  qPerBase = as(quality(fq), "matrix")
  qcount = rowSums( qPerBase <= 20)
  writeFastq(fq[qcount == 0],
             paste(fastqFile, "Qfiltered", sep="_"),
             mode="a")
}

#将ShortRead每一步的过程合起来，设置While循环，调用yield()来遍历所有片段化文件。
#去除<20 reads的基因。