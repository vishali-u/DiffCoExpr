#' Prepare the scRNA-seq data for analysis by normalizing, clustering, and
#' assigning cell type to cluster if applicable
#'
#' @param geneMatrixPath A path to the directory containing filtered gene 
#'    barcode matrices produced using the cellranger pipeline (as a reminder, 
#'    cellranger will produce a directory containing 3 files: matrix.mtx, 
#'    genes.tsv/features.tsv, and barcodes.tsv)
#'    
#' @param cellTypesPath A path to a csv file that maps marker genes to cell 
#'    type; default value is NULL (the clusters are not assigned a cell type).
#'    The first column must be the marker genes and the second column must be
#'    cell type labels
#'    
#' @return Returns an object of class SeuratObject that stores the data for
#'    a single cell RNA-seq experiment
#' 
#' # Using pbmc dataset available with package
#' # Takes too long to run
#' @examples
#' \dontrun{
#' geneMatrixPath <- system.file("extdata",
#'                               "filtered_gene_bc_matrices",
#'                                package =  "DiffCoExpr")
#' cellTypes <- system.file("extdata", 
#'                          "cellTypes.csv", 
#'                           package = "DiffCoExpr")
#' 
#' pbmc <- prepareData(geneMatrixPath = geneMatrixPath,
#'                     cellTypesPath = cellTypes)
#' pbmc
#' }
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
#' @importFrom utils read.csv
prepareData <- function(geneMatrixPath, 
                        cellTypesPath=NULL) {
  
  # --- Check some conditions to stop early if necessary --------------------
  
  if (!dir.exists(geneMatrixPath)) {
    stop("The gene matrices path you have provided does not exist or is not a
         directory.")
  }
  
  if (is.null(cellTypesPath)) {
    warning("You have not provided a file mapping cell type markers to cell
             types. Cell types will not be added, and cluster names c0, c1,..., 
             will be used instead.")
  }
  
  if (!is.null(cellTypesPath) && !file.exists(cellTypesPath)) {
    stop("You have tried to provide a file mapping cell type markers to cell 
          types. However, the file does not exist. Please re-run with the a
          valid file path or without adding a file path.")
  }
  
  # Check if the needed files are in the provided directory
  # Note: depending on which version of cellranger was used, there output will 
  # have either barcodes.tsv or features.tsv
  matrixPath <- file.path(geneMatrixPath, "matrix.mtx")
  genesFile <- file.path(geneMatrixPath, "genes.tsv")
  barcodesFile <- file.path(geneMatrixPath, "barcodes.tsv")
  featuresFile <- file.path(geneMatrixPath, "features.tsv")
  
  matrixExists <- file.exists(matrixPath)
  genesExists <- file.exists(genesFile)
  barcodesOrFeaturesExists <- file.exists(barcodesFile) || 
    file.exists(featuresFile)
  
  if (!(matrixExists && genesExists && barcodesOrFeaturesExists)) {
    stop("At least one of the needed matrices is missing. Ensure the directory
         you provide contains matrix.mtx, genes.tsv, and either barcodes.tsv
         or features.tsv.")
  }
  
  # --- Start running Seurat pipeline ---------------------------------------
  
  # Load the data and store the non-normalized counts in a Seurat object
  srat.data <- Seurat::Read10X(data.dir = geneMatrixPath)
  srat <- Seurat::CreateSeuratObject(counts = srat.data,
                                     min.cells = 3,
                                     min.features = 200)
  
  # Calculate percentage of all counts belonging to a subset of the features 
  # for each cell (as a quality control). Filter results based on this.
  srat[["percent.mt"]] <- Seurat::PercentageFeatureSet(srat, 
                                                       pattern = "^MT-")
  srat <- 
    subset(srat, 
           subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  
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
                         features = features,
                         verbose = FALSE)
  
  # Cluster the cells
  srat <- Seurat::FindNeighbors(srat, dims = 1:10)
  srat <- Seurat::FindClusters(srat, resolution = 0.5, verbose = FALSE)
  
  # Perform non-linear dimensional reduction
  srat <- Seurat::RunUMAP(srat, dims = 1:10)
  
  # If cellTypesPath is NULL, we can return here (the cell types will not be
  # annotated and downstream analyses will be done based on cluster)
  if (is.null(cellTypesPath)) {
    return(srat)
  }
  
  # --- Add cell type names based on cellTypes file ---------------------------
  
  # Find all marker genes for each cluster (to annotate with cell identity)
  srat.markers <- Seurat::FindAllMarkers(srat, 
                                         only.pos = TRUE, 
                                         verbose = FALSE)
  
  # Read in the cellTypes file and store data in a list
  originalCellTypes <- utils::read.csv(file = cellTypesPath)
  cellTypes <- strsplit(as.character(originalCellTypes[,1]), ",\\s*")
  names(cellTypes) <- originalCellTypes[,2]
  
  # Print the cell types to the console
  cellTypesString <- paste(names(cellTypes), collapse = ", ")
  message("Cell Types: ", cellTypesString)
  
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

# [END]