#' Prepare the scRNA-seq data for analysis by normalizing, clustering, and
#' assigning cell type to cluster if applicable 
#'
#' @param gene_bc_matrix_path A path to filtered gene barcode matrices produced
#'    using the cellranger pipeline.
#' @param cell_types_path A path to a csv file that maps marker genes to cell 
#'    type; default value is NULL (the clusters are not assigned a cell type).
#'    
#' @return Returns an object of class SeuratObject that stores the data for
#'    a single cell RNA-seq experiment
#' 
#' @examples
#' # Using raw pbmc dataset available with package
#' TODO: use system.file
#' gene_bc_matrix <- file.path('inst/extdata/pbmc/filtered_gene_bc_matrices')
#' cell_types <- file.path('inst/extdata/pbmc/cell_types.csv')
#' 
#' pbmc <- prepare_data(gene_bc_matrix_path = gene_bc_matrix,
#'                      cell_types_path = cell_types)
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
#' @import scCustomize
#' 
prepareData <- function(gene_bc_matrix_path, 
                         cell_types_path=NULL) {
  
  # Load the data and store the non-normalized counts in a Seurat object
  srat.data <- Seurat::Read10X(data.dir = gene_bc_matrix_path)
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
  
  # If cell_types is NULL, we can return here (the cell types will not be
  # annotated and downstream analyses will be done based on cluster)
  if (is.null(cell_types_path)) {
    return(srat)
  }
  
  # Find all marker genes for each cluster (to annotate with cell identity)
  srat.markers <- Seurat::FindAllMarkers(srat, only.pos = TRUE)
  
  # Read in the cell_types file and store data in a list
  original_cell_types <- read.csv(file = cell_types_path)
  cell_types <- strsplit(as.character(original_cell_types$Markers), ",\\s*")
  names(cell_types) <- original_cell_types$Cell.Type
  
  # Create a data frame mapping cluster ID to a comma separated string of all
  # marker genes in that string
  markers_subset <- srat.markers[, c("cluster", "gene")]
  cluster_to_gene <- aggregate(gene ~ cluster, 
                               data = markers_subset, 
                               FUN = function(x) { list(x) })
  
  # Initialize values to store cell identities in the cluster data frame
  cluster_to_gene$cell_identity <- NA
  num_cell_types <- length(cell_types)
  num_clusters <- nrow(cluster_to_gene)
  
  # Iterate over each cluster
  for (i in seq_len(num_clusters)) {
    curr_cluster_genes <- cluster_to_gene$gene[[i]]
    
    # Iterate over each cell type
    for (j in seq_len(num_cell_types)) {
      curr_type <- names(cell_types)[j]
      curr_markers <- cell_types[[j]]
      
      # Check if all markers for current cell type are in the cluster's genes
      # Assign the cell type if all markers are present
      if (all(curr_markers %in% curr_cluster_genes)) {
        cluster_to_gene$cell_identity[j] <- curr_type
        # Stop checking after the first match, assuming one identity per cluster
        break 
      }
    }
  }
  
  # Make sure the data frame is ordered from smallest to largest cluster number
  cluster_to_gene <- cluster_to_gene[order(cluster_to_gene$cluster), ]
  
  # Rename the identies to the cell types in the Seurat object
  new_cluster_ids <- cluster_to_gene$cell_identity
  names(new_cluster_ids) <- levels(srat)
  srat <- Seurat::RenameIdents(srat, 
                               new_cluster_ids)
  # TODO: rename the clusters maybe
  # srat <- scCustomize::Rename_Clusters(srat, 
  #                                      new_cluster_ids,
  #                                      meta_col_name = "Original Clusters")
  
  return(srat)
}





