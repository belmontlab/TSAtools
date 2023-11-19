#' Slope Calculation of TSA Signals
#'
#' This function calculates the slope of TSA signals (or any other continuous genomic signals)
#' using a cubic polynomial. It works by moving a window along the genome, 
#' fitting a cubic polynomial within that window, and computing the slope at the central bin.
#'
#' @param BW Input BigWig file represented as a GRanges object.
#' @param bin Bin window size for initial data binning. 
#'           If not specified (NA), the original data binning is used. 
#'           For example, for 20KB binning, specify 2e4. Default is NA.
#' @param slop_wind The number of bins to include in the slope calculation window.
#' @param deg Degree of the polynomial used for slope calculation. 
#'            Typically, deg=3 is used for cubic polynomial fitting.
#'
#' @return A GRanges object with calculated slopes for each bin.
#'
#' @export
#'
#' @example
#' x <- TSA_slopr_var(bigwig, bin=2e4, slop_wind=20, deg=3)
#'
BW_slopr<-function(BW,bin,slop_wind,deg){
  #bin= How to bin the original data, if you want to bin every 20KB it would be
  #2e4. slop_wind=
  require(rtracklayer)
  require(GenomicRanges)
  require(parallel)
  
  BW_bnr<-function(BW,bin){
    bins<-tileGenome(seqinfo(BW),tilewidth = bin,cut.last.tile.in.chrom = TRUE)
    BW_Rle<-coverage(BW,weight="score")
    BW_bin<-binnedAverage(bins,BW_Rle,"score")
    return(BW_bin)
  }
  
  if (!is.na(bin)){
    BW_bn<-BW_bnr(BW,bin)
  } else{
    BW_bn<-BW
  }
  
  # Calculate the number of cores and leave 1 alone for computer to breath!
  no_cores <- detectCores() - 1
  BW_df<- mcols(BW_bn) #extract score and store it in dataframe
  l<-length(BW_bn) # number of bins for calculation
  wind<<-slop_wind # window to fit polynomial. It is defined as a global variable because it is passed to slop function. See below.
  deg<<-deg
  ind<-matrix(1:l,nrow=l,ncol=1) #have an index
  
  seq2 <- Vectorize(seq.default, vectorize.args = c("from", "to")) #vectorize seq function 
  
  
  indm <-c(seq2(ind-(wind/2),ind+(wind/2),1)) #create indexes for polynomial fitting bins
  indm[indm<1]<-NA #remove negative values
  indm[indm>l]<-NA # remove (indexes > length)
  #extract scores for polynomial fitting
  fit_bns<-matrix(BW_df$score[indm],nrow=l,ncol=wind+1,byrow = TRUE) 
  
  
  # Initiate cluster
  cl <- makeCluster(no_cores, type="FORK")
  
  #using parApply to both vectorize and parallelize my code.
  slop<-function(x){
    #fits a degree 3 polynomial
    if (any(is.na(x))|any(x==0)){
      return (0)
    } else{
      c<-lm(x~poly(1:(wind+1),degree = deg,raw = TRUE))$coefficients
      return(abs(c[2]))
    }
    
  }
  fits<-parApply(cl,fit_bns,MARGIN = 1,FUN=slop)
  
  stopCluster(cl)
  
  
  temp<-BW_bn
  temp$score<-fits
  
  return(temp)
}