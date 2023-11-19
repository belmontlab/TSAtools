#' Convert BigWig Values to Percentiles
#'
#' This function transforms the score values in a BigWig file into their corresponding percentiles. 
#' Users have the option to apply a cutoff threshold and to export the resulting data back into a BigWig file.
#'
#' @param bigwig A BigWig file path or a GRanges object representing the BigWig data.
#' @param cut_off A numeric value representing the percentile cutoff threshold. 
#'                Scores below this percentile will be excluded. If FALSE (default), no cutoff is applied.
#' @param bigwig_export A logical flag indicating whether to export the result as a BigWig file. 
#'                      If TRUE, the resulting data is exported. If FALSE (default), the function returns a GRanges object.
#'
#' @return If `bigwig_export` is FALSE, returns a GRanges object with scores transformed to percentiles.
#'         If `bigwig_export` is TRUE, exports a BigWig file and returns NULL.
#'
#' @keywords genomics, BigWig, percentile, transformation
#'
#' @examples
#' # Convert BigWig scores to percentiles without a cutoff and return as a GRanges object
#' bw_perc("path/to/bigwig.bw")
#'
#' # Apply a 50th percentile cutoff and return as a GRanges object
#' bw_perc("path/to/bigwig.bw", cut_off = 0.50)
#'
#' # Convert to percentiles and export as a new BigWig file
#' bw_perc("path/to/bigwig.bw", bigwig_export = TRUE)
#'

bw_perc<-function(bigwig,cut_off=FALSE,bigwig_export=FALSE){
  library(rtracklayer)
  library(GenomeInfoDb)
  bw<-as.gr(bigwig, genome=NULL)
  perc.rank <- function(x) trunc(rank(x))/length(x)
  temp<-bw[bw$score!=0]
  temp$score<-perc.rank(temp$score)
  cutoff_file_ext<-""
  if (cut_off){
    temp<-temp[temp$score>cut_off]
    cutoff_file_ext<-trunc(cut_off*100)
  }
  if (bigwig_export){
      export(object = temp,con = paste(tools::file_path_sans_ext(bigwig_path),"_",
                                       cutoff_file_ext,"pct.bw",sep='')
      )
    }else return(temp)


  }
