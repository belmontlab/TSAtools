#'as.gr
#'
#' This function imports/trim a genomic file/granges and
#' make it comply to match certain genome chrom.sizes.
#' Convert BigWig File or GRanges Object to Genome-Compliant GRanges
#'
#' `as.gr` imports a genomic file (e.g., BigWig) or a GRanges object and modifies it to 
#' comply with a specified genome chromosomal sizes. It allows the option to keep only 
#' standard chromosomes and trim the rest according to the specified genome assembly version.
#'
#' @param input Either a path to a BigWig file or a GRanges object.
#' @param keepSTDchroms Logical; if TRUE, only standard chromosomes are kept and others are trimmed. 
#'                     Default is FALSE.
#' @param genome A string specifying the genome assembly version (e.g., "hg38").
#' @return A GRanges object that is compliant with the specified genome chromosomal sizes.
#'
#' @seealso \code{\link{GRanges}}, \code{\link{import}}
#'
#' @export
#'
#' @examples
#' # Importing a BigWig file according to hg38 chromosomal sizes
#' bw <- as.gr(bigwig_path, keepSTDchroms = FALSE, genome = "hg38")
#'
#' # Trimming a GRanges object to only standard chromosomes according to hg38
#' bw_trimmed <- as.gr(bw, keepSTDchroms = TRUE, genome = "hg38")
#'

as.gr <- function(input, keepSTDchroms = FALSE, genome = NULL) {
  library(regioneR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(rtracklayer)
  
  if (is.character(input)) {
    granged <- toGRanges(import(input), genome = genome)
  } else if (is(input, "GRanges")) {
    granged <- toGRanges(input, genome = genome)
  }
  
  if (keepSTDchroms) {
    granged <- keepStandardChromosomes(granged, species = "Homo sapiens", pruning.mode = "coarse")
    granged <- dropSeqlevels(granged, "chrM", pruning.mode = "coarse")
  }
  
  granged <- trim(granged)
  return(granged)
}

