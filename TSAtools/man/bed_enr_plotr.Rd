% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bed_enr_plottr.R
\name{bed_enr_plotr}
\alias{bed_enr_plotr}
\title{BED Enrichment Plotter}
\usage{
bed_enr_plotr(BED, genome = "hg38", center_file, bed_title, center_name)
}
\arguments{
\item{BED}{A string specifying the path to the BED file containing the ChIP peaks.}

\item{genome}{A string specifying the genome assembly (e.g., "hg38" or "hg19"). Default is "hg38".}

\item{center_file}{A string specifying the path to the BED file containing center genomic locations.}

\item{bed_title}{A string representing the title for the BED dataset in the plot.}

\item{center_name}{A string representing the name of the center locations used in the plot.}
}
\value{
A ggplot object representing the peak enrichment plot.
}
\description{
This function generates an enrichment plot of ChIP-seq called peaks centered on specific genomic locations. 
It requires a BED file containing the ChIP peaks, a genome assembly, and a center file with genomic locations.
}
\examples{
CTCF1 <- bed_enr_plotr(
  BED = "gnmc/CTCF_all_ENCFF175EOP_hg38.bed",
  genome = "hg38",
  center_file = "gnmc/HSlp_center_refc1r1_hg38.bed",
  bed_title = "CTCF Peaks",
  center_name = "HSlp Center"
)

}
