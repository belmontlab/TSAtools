#' Convert TSA-Score to Distance Using Specific Constants
#'
#' This function takes a GRanges instance containing TSA-scores from a BigWig file 
#' and converts these scores into distances based on provided constants A, B, and R. 
#' The formula used for conversion is derived from specific TSA-based calculations.
#' For more details see:
#' Chen, Yu et. al “Mapping 3D Genome Organization Relative to Nuclear Compartments 
#' Using TSA-Seq as a Cytological Ruler.” The Journal of Cell Biology 217,
#' no. 11 (November 5, 2018): 4025–48. https://doi.org/10.1083/jcb.201807108.
#'
#' @param BW_GRanges A GRanges object representing the TSA-score BigWig file. 
#'                   It should contain a score column with TSA-scores.
#' @param A A constant used in the TSA to distance conversion formula. 
#'          This parameter impacts the baseline adjustment in the calculation.
#' @param B Another constant used in the conversion formula, 
#'          contributing to the normalization of the TSA-score.
#' @param R The rate constant in the conversion formula, 
#'          which affects the rate of change from TSA-score to distance.
#'
#' @return A GRanges object similar to the input but with TSA-scores 
#'         converted to distances using the specified constants.
#'
#' @export
#'
#' @examples
#' # Assuming BW_GRanges is a GRanges object with TSA-scores:
#' converted_distances <- TSA_dis_calc(BW_GRanges, A = 0.5, B = 2, R = 1.5)
#' # This will return a GRanges object where the scores are now distances.
#'
#' @keywords genomics, TSA-score, distance conversion
#' 
TSA_dis_calc<- function(BW_GRanges,A,B,R){
  TSA_score<-BW_GRanges$score
  dist<-(log((2^TSA_score)-A, base = exp(1))/R)-(log(B,base=exp(1))/R)
  dist_bw<-BW_GRanges
  dist_bw$score<-dist
  return(dist_bw)
}
