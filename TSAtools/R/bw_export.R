#' Export GRanges Metadata Columns as BigWig Files
#'
#' This function facilitates exporting variables stored in the metadata columns (mcols) of a GRanges object as BigWig files. 
#' Users can choose to export either all variables or a specified subset.
#'
#' @param grange A GRanges object whose metadata columns are to be exported.
#' @param prefix A string representing the prefix to be used for naming the exported BigWig files. 
#'               This prefix is combined with the variable names to create the file names.
#' @param vars A character vector specifying which variables (columns) from the metadata of the GRanges 
#'             object should be exported. If set to "all" (default), all variables are exported.
#'
#' @return This function does not return a value but writes BigWig files to the current working directory
#'         or specified path.
#'
#' @keywords genomics, BigWig, GRanges, export
#' @export
#'
#' @examples
#' # Assuming 'gr' is a GRanges object with metadata columns
#' # Export all variables in metadata columns as BigWig files with the prefix "sample"
#' bw_export(gr, prefix = "sample")
#'
#' # Export only specific variables
#' bw_export(gr, prefix = "sample", vars = c("var1", "var2"))
#'
bw_export<-function(grange,prefix,vars="all"){
  # exports the mcols of GRANES as bw
  if (vars=="all"){
    mcol_names<-names(mcols(grange))
  } else {
    mcol_names<-names(mcols(grange))[which(names(mcols(grange))%in%vars)]
  }
  for (name in mcol_names){
    exprt<-grange
    mcols(exprt)<-mcols(grange)[,name]
    colnames(mcols(exprt))<-"score"
    export(exprt,paste0(prefix,"-",name,".bw"))
  }
}
