% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare.R
\name{get.exons.geuvadis}
\alias{get.exons.geuvadis}
\title{Get exon data (Geuvadis)}
\usage{
get.exons.geuvadis(path = getwd())
}
\arguments{
\item{path}{data directory}
}
\description{
This function transforms exon data (virtual machine):
retains exons on the autosomes,
labels exons in the format "chromosome_start_end",
extracts corresponding gene names
}
\examples{
path <- "/virdir/Scratch/arauschenberger/trial"

}
