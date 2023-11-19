#' Average Score Calculation for Genomic Ranges
#'
#' This function computes the average score for each feature in a set of query genomic ranges (`loc.gr`) based on the scores available in a set of subject genomic ranges (`score.gr`). It uses a specified score column in `score.gr` for the calculation.
#'
#' @param loc.gr A GRanges object representing the query genomic ranges.
#' @param score.gr A GRanges object representing the subject genomic ranges, which contains the scores.
#' @param score.col A string specifying the column name in `score.gr` from which the average value is to be calculated. The default is "score".
#' @return A GRanges object (`loc.gr`) with an additional column representing the average score for each range.
#'
#' @export
#'
#' @examples
#' # Example usage of bed.avg.score
#' avg_scores <- bed.avg.score(loc.gr = x, score.gr = y, score.col = "score")
#'
bed.avg.score<-function(loc.gr,score.gr,score.col="score"){
	library(data.table)
  # 
  ov <- as.matrix(findOverlaps(loc.gr,score.gr))
  
  dt=data.table(id=ov[,1],score=values(score.gr)[ov[,2],which(names(values(score.gr))==score.col)])
  dt=dt[,list(av=mean(score)),by=id]
  names(dt)[2]<-paste0(score.col)
  res=loc.gr[dt$id,]
  mcols(res)=DataFrame(mcols(res),dt[,2])
  res
}
