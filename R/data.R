#' Single-cell RNA-seq data from peripheral blood mononuclear cells (PBMCs) 
#' from a healthy donor that was published May 26th, 2016. The name of the
#' dataset is 3k PBMCs from a Healthy Donor. 
#' 
#' extdata/filtered_gene_bc_matrices contains the matrices from this experiment
#' 
#' @source 10X Genomics
#' 
#' @format A data frame with 9 observations and 2 variables
#' \describe{
#'  \item{Cell.Type}{The name of the cell type}
#'  \item{Markers}{The genes that were used as markers to identify a cell
#'    as belonging to the corresponding cell type}
#' }
#' 
#' @references 
#' Butler, A. (2015). Seurat: Tools for Single Cell Genomics
#' R package version 4.4.0
#' data \href{https://cran.r-project.org/web/packages/Seurat/}{Link}.
#' 
#' @examples
#' \dontrun{
#' cellTypes
#' }
"cellTypes"

#' Example expression matrix that was generated using the same data as 
#' described above. To generate the data, we ran prepareData using the example
#' data in the package, and then ran getExpressionMatrix with the cell type
#' set to "platelet". This object was included to make the examples in the
#' function descriptions run faster.
#' 
#' @format A dataframe with 2000 observations and 13 variables
#' \describe{
#'  each observation is a gene, and each variable is a cell barcode    
#' }
#' 
#' @examples
#' \dontrun{
#' expressionMatrixPlatelet
#' }
"expressionMatrixPlatelet"

#' Another example expression matrix that was generated using the same PBMC
#' data. To generate the data, we ran prepareData using the example
#' data in the package, and then ran getExpressionMatrix with the cell type
#' set to "DC". This object was included to make the examples in the
#' function descriptions run faster.
#' 
#' @format A dataframe with 2000 observations and 344 variables
#' \describe{
#'  each observation is a gene, and each variable is a cell barcode    
#' }
#' 
#' @examples
#' \dontrun{
#' expressionMatrixDC
#' }
"expressionMatrixDC"

#' An example correlation matrix that was generated using the same PBMC
#' data. To generate the data, we ran prepareData using the example
#' data in the package, getExpressionMatrix with the cell type
#' set to "Platelet", and finally ran getCorrelationMatrix with this expression
#' matrix (other arguments were left as the default). This object was included 
#' to make the examples in the function descriptions run faster.
#' 
#' @format A dataframe with 646 observations and 646 variables
#' \describe{
#'  each observation and variable are genes, and each cell stores the 
#'  correlation between the genes
#' }
#' 
#' @examples
#' \dontrun{
#' correlationMatrixPlatelet
#' }
"correlationMatrixPlatelet"

#' Another example correlation matrix that was generated using the same PBMC
#' data. To generate the data, we ran prepareData using the example
#' data in the package, getExpressionMatrix with the cell type
#' set to "DC", and finally ran getCorrelationMatrix with this expression
#' matrix (other arguments were left as the default). This object was included 
#' to make the examples in the function descriptions run faster.
#' 
#' @format A dataframe with 646 observations and 646 variables
#' \describe{
#'  each observation and variable are genes, and each cell stores the 
#'  correlation between the genes
#' }
#' 
#' @examples
#' \dontrun{
#' correlationMatrixDC
#' }
"correlationMatrixDC"

#' An example coexpression network that was generated using the same PBMC
#' data. To generate the data, we ran prepareData using the example
#' data in the package, getExpressionMatrix with the cell type
#' set to "Platelet", getCorrelationMatrix with the expression, and finally 
#' getCoexpressionNetwork with the correlation matrix. This object was included 
#' to make the examples in the function descriptions run faster.
#' 
#' @format A dataframe with 47895 observations and 3 variables
#' \describe{
#'  \item{Gene1}{One gene in the pair.}
#'  \item{Gene2}{Second gene in the pair.}
#'  \item{Correlation}{The correlation between gene1 and gene2.}
#' }
#' 
#' @examples
#' \dontrun{
#' coexpressionNetworkPlatelet
#' }
"coexpressionNetworkPlatelet"

#' Another example coexpression network that was generated using the same PBMC
#' data. To generate the data, we ran prepareData using the example
#' data in the package, getExpressionMatrix with the cell type
#' set to "DC", getCorrelationMatrix with the expression, and finally 
#' getCoexpressionNetwork with the correlation matrix. This object was included 
#' to make the examples in the function descriptions run faster.
#' 
#' @format A dataframe with 47895 observations and 3 variables
#' \describe{
#'  \item{Gene1}{One gene in the pair.}
#'  \item{Gene2}{Second gene in the pair.}
#'  \item{Correlation}{The correlation between gene1 and gene2}
#' }
#' 
#' @examples
#' \dontrun{
#' coexpressionNetworkDC
#' }
"coexpressionNetworkDC"



