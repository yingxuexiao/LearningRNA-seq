
##quality check on sequencing reads
#预准备确定文件具体位置
dirPath <- "C:/Users/Administrator/Documents/GitHub/LearningRNA-seq/data/E-MTAB-1147" 

#A  利用rqc()抓取对应的fastq文件并且得到质量检测报告html版

library(Rqc)
qcRes=rqc(path=dirPath,pattern=".fastq.gz",openBrowser = FALSE)
#得到html报告：'C:\Users\Administrator\AppData\Local\Temp\
#RtmpM5vw6e/rqc_report.html' has been created.

#B 利用对应的函数得到报告中对应的数据
# 1.获取单碱基测序的质量报告，帮助考虑是否Trim得分低于20的reads，
#也能帮助检测是否存在接头污染
rqcCycleQualityBoxPlot(qcRes)

# 2.获取每个位置核苷酸碱基的百分比的图像
rqcCycleBaseCallsLinePlot(qcRes)

#D  检查FASTQ文件中不同重复层次的百分比，一般的reads只复制一次
rqcReadFrequencyPlot(qcRes)




##较为广泛进行测序质量控制的工具为FASTQC，R中使用fastqcr
library(fastqcr)
#由于不是mac linux系统，所以需进行操作，让fastqcr用于windows系统：
.check_if_unix <<- function() {
  return(NULL)
}

assignInNamespace(".check_if_unix", .check_if_unix, ns = "fastqcr")

wsl.exe /path_to_wsl_fastqc_install/fastqc

fastqcr::fastqc("fastq_path", fastqc.path = "wsl.exe /path_to_wsl_fastqc_install/fastqc")

#安装fastqcr相关的java工具
fastqc_install()

#获取处理后的结果并将其放在对应的文件夹中(新创一个叫results的文件夹)
fastqc(fq.dir=dirPath,qc.dir="fastqc_results")

#进行fastqc分析，读取已处理的报告
qc_report(qc.path="fastqc_results",
          result.file = "reportFile",preview = TRUE)

qc <- qc_read("fastqc_results/ERR127302_1_subset.fastq.gz")
qc

qc_plot(qc,"Per base sequence quality")




#####过滤、去除reads
library(QuasR)
#获取fastq 文件的路径
fastqFiles <- system.file(package = "ShortRead",
                          "extdata/E-MTAB-1147",
                          c("fastqc_results/ERR127302_1_subset.fastq.gz",
                            "fastqc_results/ERR127302_2_subset.fastq.gz"))
#
outfiles <- paste(tempfile(pattern=c("processed_1_","processed_2_")),
                  ".fastq",sep="")
preprocessReads(fastqFiles,outfiles,nBases = 1,
                truncateStartBases = 3,
                Lpattern = "ACCCGGGA",
                minLength = 40)
#报错 找不到处理过的数据文件，并且最后进行过滤也报错。




####利用ShortRead函数进行数据处理：读取一个fastq文件，并过滤每个质量分数低于20的读数。
library(ShortRead)
#A 读取fastq文件
fastqFile <-system.file(package = "ShortRead","extdata/E-MTAB-1147",
                        "ERR127302_1_subset.fastq.gz") 

#B 读取fastq file
fq <- readFastq(fastqFile)

#C  利用矩阵得到每个碱基的质量得分
qPerBase = as(quality(fq), "matrix")

#D  获取碱基质量得分<20的数量
qcount = rowSums( qPerBase <= 20)

#E  得到菲尔德质量评分>=20的reads 数量
fq[qcount == 0]

#F  读写出每个碱基中所有read得分高于20的fastq 文件
writeFastq(fq[qcount == 0],
           paste(fastqFile, "Qfiltered", sep="_"))

#G  fastq文件很大，常用的方法是通过片段化逐条读取
## set up streaming with block size 1000，帮助我们顺利读取片段化的reads
f <- FastqStreamer(fastqFile,readerBlockSize=1000)

while(length(fq <- yield(f))) {
  qPerBase = as(quality(fq), "matrix")
  qcount = rowSums( qPerBase <= 20)
  writeFastq(fq[qcount == 0],
             paste(fastqFile, "Qfiltered", sep="_"),
             mode="a")
}
#将ShortRead每一步的过程合起来，设置While循环，调用yield()来遍历所有片段化文件。
#去除<20 reads的基因。




####映射Mapping/aligning reads to the genome举例
library(QuasR)

#A  复制示例的数据到当前的工作目录路径
#文件本身在QuasR包中，从C:/program file放至当前目录dirPath中
file.copy(system.file(package="QuasR", "extdata"),
          ".", recursive=TRUE)

#B  从文件中打开hg19sub.fa并命名
genomeFile <- "extdata/hg19sub.fa"

#C  样本文件中包含了实例名称和fastq文件路径
sampleFile <- "extdata/samples_chip_single.txt"

#D  
proj <- qAlign(sampleFile, genomeFile)
proj




