% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test.R
\name{test.single}
\alias{test.single}
\title{Conduct single test}
\usage{
test.single(Y, X, map, i, limit = NULL, steps = NULL, W = NULL,
  w.type = NULL)
}
\arguments{
\item{Y}{exon expression\strong{:}
matrix with \eqn{n} rows (samples) and \eqn{p} columns (exons)}

\item{X}{SNP genotype\strong{:}
matrix with \eqn{n} rows (samples) and \eqn{q} columns (SNPs)}

\item{map}{list with names "genes", "exons", and "snps"
(output from \code{\link{map.genes}}, \code{\link{map.exons}},
and \code{\link{map.snps}})}

\item{i}{gene index\strong{:}
integer between \eqn{1} and \code{nrow(map$genes)}}

\item{limit}{cutoff for rounding \code{p}-values}

\item{steps}{size of permutation chunks\strong{:}
integer vector}

\item{W}{numeric matrix\strong{:}
a square matrix with as many rows as there are covariates in the independent data set
  It represents the correlation structure expected in the independent data}

\item{w.type}{string\strong{:}
if \code{W=NULL} and w.type is a given string, take the type: "cov" is the only one allowed. 
Then the inner product of the indep data \code{X} is taken}
}
\description{
This function tests for alternative splicing.
}
\details{
The maximum number of permutations equals \code{sum(steps)}. Permutations is
interrupted if at least \code{limit} test statistics for the permuted data
are larger than the test statistic for the observed data.
}
\examples{
NA

}
