#' Prepare the scRNA-seq data for analysis by normalizing, clustering, and
#' assigning cell type to cluster if applicable 
#'
#' @param geneMatrixPath A path to filtered gene barcode matrices produced
#'    using the cellranger pipeline.
#' @param cellTypesPath A path to a csv file that maps marker genes to cell 
#'    type; default value is NULL (the clusters are not assigned a cell type).
#'    
#' @return Returns an object of class SeuratObject that stores the data for
#'    a single cell RNA-seq experiment
#' 
#' @examples
#' # Using pbmc dataset available with package
#' geneMatrixPath <- system.file('extdata/pbmc/filtered_gene_bc_matrices')
#' cellTypes <- system.file('extdata/pbmc/cell_types.csv')
#' 
#' pbmc <- prepareData(geneMatrixPath = geneMatrixPath,
#'                     cellTypesPath = cellTypes)
#' pbmc
#' 
#' @references
#' Butler A, Hoffman P, Smibert P, Papalexi E, Satija R (2018). 
#' “Integrating single-cell transcriptomic data across different conditions, 
#' technologies, and species.” Nature Biotechnology, 36, 411-420. 
#' \href{https://doi.org/10.1038/nbt.4096}{Link}.
#' 
#' Butler, A. (2015). Seurat: Tools for Single Cell Genomics
#' R package version 4.4.0
#' data \href{https://cran.r-project.org/web/packages/Seurat/}{Link}.
#' 
#' @export
#' @import Seurat
#' 
prepareData <- function(geneMatrixPath, 
                         cellTypesPath=NULL) {
  
  # Load the data and store the non-normalized counts in a Seurat object
  srat.data <- Seurat::Read10X(data.dir = geneMatrixPath)
  srat <- Seurat::CreateSeuratObject(counts = srat.data,
                                     min.cells = 3,
                                     min.features = 200)
  
  # Calculate percentage of all counts belonging to a subset of the features 
  # for each cell (as a quality control). Filter results based on this.
  srat[["percent.mt"]] <- Seurat::PercentageFeatureSet(srat, 
                                                       pattern = "^MT-")
  srat <- subset(srat, 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                   percent.mt < 5)
  
  # Normalize expression measurements for each cell (using default
  # normalization methods from the Seurat package)
  srat <- Seurat::NormalizeData(srat)
  
  # Scale the data and then perform dimensional reduction (pca)
  srat <- Seurat::FindVariableFeatures(srat,
                                       selection.method = "vst",
                                       nfeatures = 2000)
  all.genes <- rownames(srat)
  srat <- Seurat::ScaleData(srat,
                            features = all.genes)
  features <- Seurat::VariableFeatures(object = srat)
  srat <- Seurat::RunPCA(srat,
                         features = features)
  
  # Cluster the cells
  srat <- Seurat::FindNeighbors(srat, dims = 1:10)
  srat <- Seurat::FindClusters(srat, resolution = 0.5)
  
  # Perform non-linear dimensional reduction
  srat <- Seurat::RunUMAP(srat, dims = 1:10)
  
  # If cellTypesPath is NULL, we can return here (the cell types will not be
  # annotated and downstream analyses will be done based on cluster)
  if (is.null(cellTypesPath)) {
    return(srat)
  }
  
  # Find all marker genes for each cluster (to annotate with cell identity)
  srat.markers <- Seurat::FindAllMarkers(srat, only.pos = TRUE)
  
  # Read in the cellTypes file and store data in a list
  originalCellTypes <- read.csv(file = cellTypesPath)
  cellTypes <- strsplit(as.character(originalCellTypes$Markers), ",\\s*")
  names(cellTypes) <- originalCellTypes$Cell.Type
  
  # Create a data frame mapping cluster ID to a comma separated string of all
  # marker genes in that string
  markersSubset <- srat.markers[, c("cluster", "gene")]
  clusterToGene <- aggregate(gene ~ cluster,
                             data = markersSubset, 
                             FUN = function(x) { list(x) })
  
  # Initialize values to store cell identities in the cluster data frame
  clusterToGene$cellIdentity <- NA
  numCellTypes <- length(cellTypes)
  numClusters <- nrow(clusterToGene)
  
  # Iterate over each cluster
  for (i in seq_len(numClusters)) {
    currClusterGenes <- clusterToGene$gene[[i]]
    
    # Iterate over each cell type
    for (j in seq_len(numCellTypes)) {
      currType <- names(cellTypes)[j]
      currMarkers <- cellTypes[[j]]
      
      # Check if all markers for current cell type are in the cluster's genes
      # Assign the cell type if all markers are present
      if (all(currMarkers %in% currClusterGenes)) {
        clusterToGene$cellIdentity[j] <- currType
        # Stop checking after the first match, assuming one identity per cluster
        break 
      }
    }
  }
  
  # Make sure the data frame is ordered from smallest to largest cluster number
  clusterToGene <- clusterToGene[order(clusterToGene$cluster), ]
  
  # Rename the identies to the cell types in the Seurat object
  newClusterIds <- clusterToGene$cellIdentity
  names(newClusterIds) <- levels(srat)
  srat <- Seurat::RenameIdents(srat, 
                               newClusterIds)
  
  return(srat)
}





