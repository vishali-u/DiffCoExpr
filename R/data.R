#' Single-cell RNA-seq data from peripheral blood mononuclear cells (PBMCs) 
#' from a healthy donor that was published May 26th, 2016. The name of the
#' dataset is 3k PBMCs from a Healthy Donor.
#' 
#' @source 10X Genomics
#' 
#' TODO: how to make a reference for data downloaded from 10X website
#' 
#' data/filtered_gene_bc_matrices contains all the matrices that are generated
#' during alignment of scRNA-seq using the cellranger pipeline
#' 
#' barcodes.rds contains a list of all the cell barcodes
#' @examples
#' \dontrun{
#'  barcodes
#' }
#' "barcodes"

#' genes.rda contains the genes that were identified during mapping
#' @examples
#' \dontrun{
#' genes
#' }
#' "genes"

#' matrix.mda contains the raw gene counts per cell
#' @examples
#' \dontrun{
#' matrix
#' }
#' "matrix"

#' cell_types.rda maps cell types to marker genes
#' 
#' @references 
#' Butler, A. (2015). Seurat: Tools for Single Cell Genomics
#' R package version 4.4.0
#' data \href{https://cran.r-project.org/web/packages/Seurat/}{Link}.
#' 
#' #' @examples
#' \dontrun{
#' cell_types
#' }
#' "cell_types"
