#!/usr/bin/env Rscript

#This script is to calculate shannon entropy over an .rds object, AFTER batch correction. 
#Entropy is calculated depending on the corrected space that enters the script.

#libraries
library(scran)

args <-R.utils::commandArgs(asValues=TRUE)
if (is.null(args[["input"]])) {
  print("Provide a valid input file name (Batch corrected object) --> RDS file")
}
if (is.null(args[["output_entropy"]])) {
  print("Provide a valid output file name (Per batch and cell type entropy) --> CSV file")
}
if (is.null(args[["output_matrix"]])) {
  print("Provide a valid output file name for the matrix. RDS file")
}
#input file
object <- readRDS(args[["input"]])
#entropy function
shannon_entropy <- function(x, batch_vector, N_batches) {
  freq_batch = table(batch_vector[x ==1])/length(batch_vector[x == 1])
  freq_batch_positive = freq_batch[freq_batch > 0]
  return(-sum(freq_batch_positive * log(freq_batch_positive))/log(N_batches))
}

#specify the arguments of the inside function in the outside function!!!!!!!!!!!
compute_entropy <- function(corrected_space, bool, x, batch_vector, N_batches, cell_type_vector, N_cell_types){
  knn_graph <- as.matrix(buildKNNGraph(corrected_space, k=30, d=50, transposed = bool)[])
  batch_entropy <- apply(knn_graph, 1, function(x) {shannon_entropy (x, batch_vector, N_batches)})
  celltype_entropy <- apply(knn_graph, 1, function(x) {shannon_entropy (x, cell_type_vector, N_cell_types)})
  entropy <- cbind(batch_entropy,celltype_entropy)
  names(entropy) <- c("batch_entropy", "cell_type_entropy")
  return(entropy)
}

save_results <- function(x, col_names){
  write.table(x, file = args[["output_entropy"]], sep = "\t", row.names = FALSE, col.names = col_names)
}

# 1) SCE objects 
if (class(object) == "SingleCellExperiment"){
  print("The input object is a Single Cell Experiment class object")
  
  #Calculate entropy
  batch_vector <- colData(object)$Batch
  N_batches <- length(unique(batch_vector))
  if ("cell_type1" %in% colnames(colData(object))){
    print("The object has cell type annotation")
    cell_type_vector <- colData(object)$cell_type1
    N_cell_types <- length(unique(cell_type_vector))
  }
  # COMPUTE ENTROPY 
  # 1.1) tools that correct low_D embeddings
  if ("corrected_embedding" %in% reducedDimNames(object)){
    corrected_space <- as.matrix(object@reducedDims@listData[["corrected_embedding"]])
    col_names <- c("PCA_batch_entropy", "PCA_cell_type_entropy")
    save_results(compute_entropy(corrected_space, bool = TRUE, x, batch_vector, N_batches, cell_type_vector, N_cell_types), col_names) # <------ FUNCTION
    print("Entropy calculated in PCA space!")
    saveRDS(corrected_space, file = args[["output_matrix"]])
    }
  
  # 1.2) methods that correct expression matrix
  if ("corrected" %in% names(assays(object))){
    corrected_space <- assay(object, "corrected")
    col_names <- c("Counts_batch_entropy", "Counts_cell_type_entropy")
    save_results(compute_entropy(corrected_space, bool = FALSE, x, batch_vector, N_batches, cell_type_vector, N_cell_types), col_names) # <------ FUNCTION
    print("Entropy calculated over counts matrix!")
    saveRDS(corrected_space, file = args[["output_matrix"]])
  }
}
#2) Seurat objects
# 2.1) Seurat v2 multiCCA
if (class(object) == "seurat"){
  library(Seurat)
  print("The input object is a Seurat class object")
  
  batch_vector <- object@meta.data$Batch
  N_batches <- length(unique(batch_vector))
  if ("cell_type1" %in% colnames(object@meta.data)){
    cell_type_vector <- object@meta.data$cell_type1
    N_cell_types <- length(unique(cell_type_vector))
  }
    
    if( "cca" %in% names(object@dr)){
        print("Object corrected by Seurat_v2_multiCCA")
        #entropy before
        cca_matrix <- attr(object@dr[["cca"]], "cell.embeddings")
        entropy_before_cca <- compute_entropy(corrected_space = cca_matrix, bool = TRUE, x, batch_vector, N_batches, cell_type_vector, N_cell_types) # <------ 
        #entropy after
        aligned_cca_matrix <- attr(object@dr[["cca.aligned"]], "cell.embeddings") 
        entropy_after_cca <- compute_entropy(corrected_space = aligned_cca_matrix, bool = TRUE, x, batch_vector, N_batches, cell_type_vector, N_cell_types) # <------ FUNCTION
  
        entropy_df <- cbind(entropy_before_cca, entropy_after_cca)
        print(class(entropy_df))
        col_names <- c("bef_cca_batch_entropy", "bef_cca_cell_type_entropy", "aft_cca_batch_entropy", "aft_cca_cell_type_entropy")
        #colnames(entropy_df) <- col_names
  
        print("Entropy calculated over CCA space of Seurat object!")
        save_results(entropy_df, col_names)
        saveRDS(aligned_cca_matrix, file = args[["output_matrix"]])
        print("Congratulations, entropy calculated over Seurat_v2_multiCCA object")
    }
}   

#2.2) Seurat v3 anchors   
if (class(object) == "Seurat"){
  library(Seurat)
  print("The input object is a Seurat class object")
  
  batch_vector <- object@meta.data$Batch
  N_batches <- length(unique(batch_vector))
  if ("cell_type1" %in% colnames(object@meta.data)){
    cell_type_vector <- object@meta.data$cell_type1
    N_cell_types <- length(unique(cell_type_vector))
  }  

    if ("integrated" %in% names(object@assays)){
        
    print("Object corrected by Seurat_v3_anchors")
        
    space <- as.matrix(object@assays[["integrated"]]@data)
    col_names <- c("Seurat_v3_Batch_entropy", "Seurat_v3_Cell_type_entropy")
    save_results(compute_entropy(corrected_space = space, bool = FALSE, x, batch_vector, N_batches, cell_type_vector, N_cell_types), col_names)# <------ 
    print("Congratulations, entropy calculated over Seurat_v3_anchors object")
    saveRDS(space, file = args[["output_matrix"]])
    } 
}
