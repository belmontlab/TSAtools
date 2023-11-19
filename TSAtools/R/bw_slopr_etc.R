#'TSA slopr
#'
#'This Function approximates TSA signal slope (or any other continous genomic signal) by a cubic polynomial. 
#'The way it works is that it move a window along the genome and fit cubic polynomial in that window and reads the slope at the centeral bin.
#'
#' @param BW input BigWig file - GRanges object
#' @param bin bin window. Defaults to 20000 bp
#' @return binned GRanges
#' @seealso 
#' @export
#' @example 
#' x <- BW_slopr(bigwig, bin=2e4)
#' 
#' 
BW_slopr<-function(BW,bin,slop_wind,abs_value=TRUE){
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
      c<-lm(x~poly(1:(wind+1),3,raw = TRUE))$coefficients
      return ((c[2]+2*c[3]*((wind/2)+1)+3*c[4]*(((wind/2)+1)^2)))
    }
    
  }
  fits<-parApply(cl,fit_bns,MARGIN = 1,FUN=slop)
  
  stopCluster(cl)
  
  
  temp<-BW_bn
  temp$score<-fits
  if(abs_value){temp$score<-abs(temp$score)}
  
  return(temp)
}