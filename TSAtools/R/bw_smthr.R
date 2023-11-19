#' Polynomial Smoothing of TSA Signals
#'
#' This function smoothens TSA signals (or any other continuous genomic signals)
#' using polynomial approximation. It operates by moving a window along the genome,
#' fitting a polynomial within that window, and reading the value at the central bin.
#' The degree of the polynomial is user-defined.
#'
#' @param BW Input BigWig file represented as a GRanges object.
#' @param bin Bin window size for initial data binning. If not specified (NA), 
#'           the original data binning is used. Default is NA.
#' @param smth_wind The number of bins to include in the polynomial smoothing window.
#' @param deg Degree of the polynomial used for smoothing. 
#'            For example, use deg=1 for linear, deg=3 for cubic smoothing.
#'
#' @return A GRanges object with smoothed scores.
#'
#' @export
#'
#' @example x <- TSA_smthr_var(bigwig, bin=NA, smth_wind=20, deg=1)
#' @example x <- TSA_smthr_var(bigwig, bin=2e4, smth_wind=20, deg=3)
#'
#' @keywords genomics, TSA-signal, smoothing
BW_smthr<-function(BW,bin,smth_wind,deg){
  #bin= How to bin the original data, if you want to bin every 20KB it would be 2e4
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
  wind<-smth_wind # window to fit polynomial
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
  smth<-function(x){
    #fits a degree 3 polynomial
    if (any(is.na(x))|any(x==0)){
      return (0)
    } else{
      return(predict(lm(x~poly(1:(wind+1),degree = deg,raw = TRUE)),data.frame(1:(wind+1)))[wind/2+1])
    }

  }

  fits<-parApply(cl,fit_bns,MARGIN = 1,FUN=smth)

  stopCluster(cl)


  temp<-BW_bn
  temp$score<-fits

  return(temp)
}
