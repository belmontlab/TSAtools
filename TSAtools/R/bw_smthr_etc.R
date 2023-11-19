#'TSA smthr
#'
#'This Function approximates TSA signal (or any other continous genomic signal) by a cubic polynomial.
#'The way it works is that it move a window along the genome and fit cubic polynomial in that window and reads the value at the centeral bin.
#' @param BW input BigWig file - GRanges object
#' @param bin bin window. Defaults to 20000 bp
#' @param smth_wind number of bins to include in polynomial smoothing
#' @return smooth GRanges
#' @seealso
#' @export
#' @example x <- BW_smthr(bigwig, bin=2e4, smth_wind=20)
#'
BW_smthr<-function(BW,bin,smth_wind){
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
      return(predict(lm(x~poly(1:(wind+1),3,raw = FALSE)),data.frame(1:(wind+1)))[wind/2+1])
    }

  }

  fits<-parApply(cl,fit_bns,MARGIN = 1,FUN=smth)

  stopCluster(cl)


  temp<-BW_bn
  temp$score<-fits

  return(temp)
}
