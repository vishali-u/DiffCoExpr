% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{coexpressionNetworkPlatelet}
\alias{coexpressionNetworkPlatelet}
\title{An example coexpression network that was generated using the same PBMC
data. To generate the data, we ran prepareData using the example
data in the package, getExpressionMatrix with the cell type set to
"Platelet", getCorrelationMatrix with the expression maxtris, and finally
getCoexpressionNetwork with the correlation matrix. This object was included
to make the examples in the function descriptions run faster.}
\format{
A dataframe with 47895 observations and 3 variables
\describe{
\item{Gene1}{One gene in the pair.}
\item{Gene2}{Second gene in the pair.}
\item{Correlation}{The correlation between gene1 and gene2.}
}
}
\usage{
coexpressionNetworkPlatelet
}
\description{
An example coexpression network that was generated using the same PBMC
data. To generate the data, we ran prepareData using the example
data in the package, getExpressionMatrix with the cell type set to
"Platelet", getCorrelationMatrix with the expression maxtris, and finally
getCoexpressionNetwork with the correlation matrix. This object was included
to make the examples in the function descriptions run faster.
}
\examples{
\dontrun{
coexpressionNetworkPlatelet
}
}
\keyword{datasets}
