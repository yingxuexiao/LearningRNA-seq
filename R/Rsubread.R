### R code from vignette source 'Rsubread.Rnw'


###################################################
library(Rsubread)
ref <- system.file("extdata","reference.fa",package="Rsubread")
buildindex(basename="reference_index",reference=ref)




###################################################
reads <- system.file("extdata","reads.txt.gz",package="Rsubread")
align.stat <- align(index="reference_index",
                    readfile1=reads,
                    output_file="alignResults.BAM",
                    phredOffset=64)



###################################################
reads1 <- system.file("extdata","reads1.txt.gz",package="Rsubread")
reads2 <- system.file("extdata","reads2.txt.gz",package="Rsubread")

align.stat2 <- align(index="reference_index",
                     readfile1=reads1,
                     readfile2=reads2,
                     
output_file="alignResultsPE.BAM",phredOffset=64)



###################################################
ann <- data.frame(
GeneID=c("gene1","gene1","gene2","gene2"),
Chr="chr_dummy",
Start=c(100,1000,3000,5000),
End=c(500,1800,4000,5500),
Strand=c("+","+","-","-"),
stringsAsFactors=FALSE)
ann

fc_SE <- featureCounts("alignResults.BAM",annot.ext=ann)
fc_SE


###################################################
fc_PE <- featureCounts("alignResultsPE.BAM",
                       annot.ext=ann,
                       isPairedEnd=TRUE)
fc_PE


###################################################
if(grepl("linux", R.version$os)) {
  md5.zip <- "ffd5036b36e25e9b61efc412e71820dd"
  URL <- "https://shilab-bioinformatics.github.io/cellCounts-Example/cellCounts-Example.zip"
  temp.file <- tempfile()
  temp.dir <- tempdir()
  downloaded <- tryCatch({
      download.file(URL, destfile = temp.file)
      tools::md5sum(temp.file) %in% md5.zip
    },
    error = function(cond){
      return(FALSE)
    }
  )
  if(!downloaded) cat("Unable to download the file.\n") 
} else downloaded <- FALSE

if(downloaded){
  unzip(temp.file, exdir=paste0(temp.dir,"/cellCounts-Example"))
  library(Rsubread)
  buildindex(paste0(temp.dir,"/chr1"), 
    paste0(temp.dir,"/cellCounts-Example/hg38_chr1.fa.gz"))
}


###################################################
if(downloaded){
  sample.sheet <- data.frame(
    BarcodeUMIFile = paste0(temp.dir,"/cellCounts-Example/reads_R1.fastq.gz"), 
    ReadFile = paste0(temp.dir,"/cellCounts-Example/reads_R2.fastq.gz"),
    SampleName="Example", stringsAsFactors=FALSE
  )
  counts <- cellCounts(paste0(temp.dir,"/chr1"), sample.sheet,  nthreads=1,
    input.mode="FASTQ", annot.inbuilt="hg38")
}


###################################################
if(downloaded) print(counts$sample.info)


###################################################
if(downloaded) print(dim(counts$counts$Example))


###################################################
x <- qualityScores(filename=reads,offset=64,nreads=1000)
x[1:10,1:10]


###################################################
propmapped("alignResults.BAM")


