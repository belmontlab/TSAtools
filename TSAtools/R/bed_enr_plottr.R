#' BED Enrichment Plotter
#'
#' This function generates an enrichment plot of ChIP-seq called peaks centered on specific genomic locations. 
#' It requires a BED file containing the ChIP peaks, a genome assembly, and a center file with genomic locations.
#'
#' @param BED A string specifying the path to the BED file containing the ChIP peaks.
#' @param genome A string specifying the genome assembly (e.g., "hg38" or "hg19"). Default is "hg38".
#' @param center_file A string specifying the path to the BED file containing center genomic locations.
#' @param bed_title A string representing the title for the BED dataset in the plot.
#' @param center_name A string representing the name of the center locations used in the plot.
#' @return A ggplot object representing the peak enrichment plot.
#'
#' @export
#'
#' @examples
#' CTCF1 <- bed_enr_plotr(
#'   BED = "gnmc/CTCF_all_ENCFF175EOP_hg38.bed",
#'   genome = "hg38",
#'   center_file = "gnmc/HSlp_center_refc1r1_hg38.bed",
#'   bed_title = "CTCF Peaks",
#'   center_name = "HSlp Center"
#' )
#'


bed_enr_plotr<- function(BED,genome="hg38",center_file,bed_title,center_name){
  library(BSgenome)
  # library(BiocInstaller)
  # biocLite("BSgenome.Hsapiens.UCSC.hg38")
  library(rtracklayer)
  library(plyr)
  library(ggplot2)
  
  #import centers
  loci<- data.frame(loci=read.table(file = center_file),type="actual")
  loci<-loci[,c(1,3,ncol(loci))]
  colnames(loci)[1]<-"chr"
  colnames(loci)[2]<-"center"
  loci$chr<-as.character(loci$chr)
  loci$type<-as.factor(loci$type)
  
  genome<-tileGenome(seqlengths(getBSgenome(genome))[1:23],tilewidth = 1e4,cut.last.tile.in.chrom = TRUE)
  # genome<-tileGenome(seqinfo(c1r1),tilewidth = 1e4,cut.last.tile.in.chrom = TRUE)
  
  ol<-as.data.frame(findOverlaps(genome,
                                 import(BED,format = "bed"),
                                 # maxgap = 0L,
                                 # minoverlap = 1L,
                                 type = "any",
                                 select ="all"))
  
  ol.count<-count(df = ol,vars = "queryHits")
  ChIP.bindingst.10KB<-genome
  ChIP.bindingst.10KB$count<-0
  ChIP.bindingst.10KB$count[ol.count$queryHits]<-ol.count$freq
  
  #location of centers on BigWig
  #loci$BWind<-0
  
  doMC::registerDoMC(cores = 6)
  loci<-adply(loci,.margins = 1,.fun = function(x){
    which((x[1,1]==seqnames(ChIP.bindingst.10KB))&
            (x[1,2]<=end(ranges(ChIP.bindingst.10KB))))[1]
  },.parallel = TRUE)
  colnames(loci)[4]<-"BWind"
  loci<- loci[!is.na(loci$BWind),] #remove NAs
  
  
  seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to")) #vectorize seq function 
  
  ind.temp<-c(seq2((loci$BWind-50),(loci$BWind+50),1))
  ind.temp[ind.temp<1]<-1
  ind.temp[ind.temp>length(ChIP.bindingst.10KB)]<-1
  actual.bins <- (data.frame(matrix(ChIP.bindingst.10KB$count[ind.temp],nrow=dim(loci)[1],ncol=101,byrow = TRUE)))
  data<-data.frame(High.Slope.avg=apply(actual.bins,MARGIN = 2,mean,na.rm=TRUE)
                   ,bin=seq(-500,500,10),Feature=bed_title)
  
  
  
  plot.info <- data.frame(x=paste("Distance from",center_name,"Kb",sep=" "),y="Peaks/10Kb", 
                          title=paste (bed_title,"vs.",center_name))
  peak.count.plot<-ggplot(data, aes(y=High.Slope.avg ,x=bin,color=Feature)) +
    geom_line(size=1)+theme_bw(base_size = 8, base_family = "")+
    xlab(plot.info$x)+ylab(plot.info$y)+
    ggtitle(plot.info$title)+
    theme(plot.title = element_text(hjust = 0.5,size=10),legend.position = "bottom")
  
  return(peak.count.plot)
}