#' High Slope Center Finder
#'
#' This function identifies the center of high-slope regions in genomic data. 
#' Two methods are available for finding the center: 'simple_center' and 'local_center'. 
#' The latter, uses the second derivative of smoothed TSA-seq data to find these regions, 
#' previously identified using the BW_slopr function. 
#'
#' @param bw_slp_path Path to the input file with slope data, typically generated 
#'        by the BW_slopr function.
#' @param quantile_cutoff The quantile used to determine high slope regions. 
#'        Defaults to 0.9, meaning the top 10% of slopes are considered high.
#' @param method Method for identifying the center of high-slope regions. 
#'        Options are 'local_center' (finds local maxima within high slope regions) 
#'        and 'simple_center' (takes the central position of high slope regions). 
#'        Default is 'local_center'.
#' @param center_width Width of the region around the center to be considered. 
#'        Defaults to 20,000 bp (2e4).
#' @param export Boolean indicating whether to export the results. 
#'        If TRUE, results are saved as a .bed file. Defaults to TRUE.
#'
#' @return A GRanges object representing the center of high-slope regions.
#'
#' @examples centers <- slp_center_find("example.bw", quantile_cutoff=0.9, method="local_center", center_width=2e4)
#'
#' @export
#' @seealso `BW_slopr`
#' @keywords genomics, TSA-seq
slp_center_find<-function(bw_slp_path,quantile_cutoff=0.9,method="local_center",center_width=2e4,export=TRUE){
  #

  library(rtracklayer)

  slp<-import(bw_slp_path)
  cut_off<-quantile(slp$score,quantile_cutoff)
  slp.high<-slp[slp$score>cut_off]
  slp.high.bed<-reduce(slp.high)

  if (method=="local_center"){
    slp$max<-0
    slp$max[(which(diff(sign(diff(slp$score)))==-2)+1)]<-1
    slp.sub<-slp[(slp$score>cut_off)&slp$max==1]

    ol<-as.data.frame(findOverlaps(slp.high.bed,
                                   slp.sub,
                                   minoverlap=1L,
                                   type="any",
                                   select ="all" ))

    ol$subscore<-slp.sub$score[ol$subjectHits] # enquire the slope values
    ol<-ol[order(ol$queryHits,-ol$subscore),]
    ol<-ol[!duplicated(ol$queryHits),]
    slp.center<-resize(slp.sub[ol$subjectHits],width = center_width, fix="center")
  }
  if (method=="simple_center"){
    slp.center<-resize(slp.high.bed,width = center_width, fix="center")
  }

  cutoff_file_ext<-paste("top",(1-quantile_cutoff)*100,sep="")

  if (export){
    export(object=slp.center, con=paste(tools::file_path_sans_ext(basename(bw_slp_path)),"_",
                                        cutoff_file_ext,method,".bed",sep=''))
    export(object = slp.high,con=paste(tools::file_path_sans_ext(basename(bw_slp_path)),"_",
                                       cutoff_file_ext,".bw",sep=''))
  }
  return(slp.center)


}
