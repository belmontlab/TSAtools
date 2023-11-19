#' BigWig File Binning
#'
#' This function takes a BigWig file represented as a GRanges object and bins it into specified intervals. 
#' It is primarily used to reduce the resolution of the data for easier analysis and visualization. 
#' The default bin size is set to 20,000 base pairs (20 kb).
#'
#' @param BW A GRanges object representing the input BigWig file.
#' @param bin The size of the binning window in base pairs (bp). 
#'           Defaults to 20,000 bp (20 kb). This parameter determines the resolution of the binning.
#' @return A GRanges object representing the binned data. Each range in the returned object corresponds 
#'         to a bin, and the associated values represent the average score of the BigWig data within that bin.
#'
#' @export
#'
#' @examples
#' # Assuming 'bigwig' is a GRanges object
#' # Example usage of BW_bnr with default bin size
#' binned_data <- BW_bnr(bigwig)
#' # Example with a custom bin size of 10,000 bp
#' binned_data_custom <- BW_bnr(bigwig, bin=1e4)
#'
BW_bnr<-function(BW,bin=2e4){
  #gets a BW file and bin it over fixed distances
  require(rtracklayer)
  require(GenomicRanges)
  bins<-tileGenome(seqinfo(BW),tilewidth = bin,cut.last.tile.in.chrom = TRUE)
  BW_Rle<-coverage(BW,weight="score")
  BW_bin<-binnedAverage(bins,BW_Rle,"score")
  return(BW_bin)
}

