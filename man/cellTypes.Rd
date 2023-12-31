% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cellTypes}
\alias{cellTypes}
\title{Single-cell RNA-seq data from peripheral blood mononuclear cells (PBMCs)
from a healthy donor that was published May 26th, 2016. The name of the
dataset is 3k PBMCs from a Healthy Donor.}
\format{
A data frame with 9 observations and 2 variables
\describe{
\item{Cell.Type}{The name of the cell type}
\item{Markers}{The genes that were used as markers to identify a cell
as belonging to the corresponding cell type}
}
}
\source{
10X Genomics
}
\usage{
cellTypes
}
\description{
extdata/filtered_gene_bc_matrices contains the matrices from this experiment
}
\details{
The raw data was downloaded from the 10X website, and the cell type lables
were taken from a Seurat tutorial in the Seurat package.
}
\examples{
\dontrun{
cellTypes
}
}
\references{
10X Genomics. 2016. “3k PBMCs from a Healthy Donor.”
data \href{https://www.10xgenomics.com/resources/datasets}{Link}.

Butler, A. (2015). Seurat: Tools for Single Cell Genomics
R package version 4.4.0
data \href{https://cran.r-project.org/web/packages/Seurat/}{Link}.
}
\keyword{datasets}
