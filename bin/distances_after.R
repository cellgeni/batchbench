#!/usr/bin/env Rscript

#libraries
library(scran)
library(SingleCellExperiment)

#BiocManager::install("lsa")
library(lsa)


#arguments
args <- R.utils::commandArgs(asValues=TRUE)
if (is.null(args[["input"]])) {
  print("Provide a valid input file name (Batch corrected object) --> RDS file")
}
if (is.null(args[["output"]])) {
  print("Provide a valid output file name (Per batch and cell type entropy) --> CSV file")
}

#input file
object <- readRDS(args[["input"]])

#subset space by batch and by cell type
subset_space <- function(space, b_vector, N_b){
  colnames(space) <- b_vector
  subset_list <- lapply(1:N_b, function(x) {space[, colnames(space) == names(table(b_vector)[x])]}) #careful with comma!
  return(subset_list)
}
#compute similarity. This function returns a list of cosine similarity symmetric matrices.
compute_similarity <- function(subset_list){
  sim_list <- list()
  for (i in 1: length(subset_list)){
    sim_list[[i]] <- cosine(x = subset_list[[i]])
  }
  return(sim_list)
}
#extract values from symmetric matrices
extract_values <- function(sim_list){
    vals <- list()
    for(n in 1: length(sim_list)){
      vals[[n]] <-  c(sim_list[[n]][lower.tri(sim_list[[n]], diag = FALSE)])
    }
    return(vals)
}

#save results function
save_results <- function(x){
  saveRDS(x, file = args[["output"]])
}


# 1) SCE objects 
if (class(object) == "SingleCellExperiment"){
    print("The object is a SingleCellExperiment class object")
  
  #arguments for similarity calculation
  batch_vector <- as.character(object$Batch)
  
  if ("cell_type1" %in% colnames(colData(object))){
    print("The object has cell type annotation")
    cell_type_vector <- as.character(object$cell_type1)
  }
  
  #redifine vectors after filtering
  batch_vector <- as.character(object$Batch)
  N_batches <- length(unique(batch_vector))
  cell_type_vector <- as.character(object$cell_type1)
  N_cell_types <- length(unique(cell_type_vector))
  
  #1.1) Distance over PCA graph
  if ("corrected_embedding" %in% reducedDimNames(object)){
  
    #embedding after correction
    space <- object@reducedDims@listData[["corrected_embedding"]]
    
    space <- t(space) #We transpose to avoid code duplication, as pca_graph is cells x PCs and exp_matrix is genes x cells. 
    similarity_batch_list <- compute_similarity(subset_list = subset_space(space = space, b_vector = batch_vector, N_b = N_batches))
    similarity_cell_type_list <- compute_similarity(subset_list = subset_space(space = space, b_vector = cell_type_vector, N_b = N_cell_types))
    similarity_list <- append(similarity_batch_list, similarity_cell_type_list)
    
    
    values_list <- extract_values(similarity_list)
    
    #convert list of lists to dataframe
    values_df <- as.data.frame(t(plyr::ldply(values_list, rbind, .id = NULL)))
    
    colnames(values_df) <- c(paste0("pca_dist_BATCH_aft_", names(table(batch_vector))), 
                     paste0("pca_dist_CELL_aft_", names(table(cell_type_vector))))
    print(summary(values_df))
    print("Distances over the PCA graph succesful!")
    
    #save results
    save_results(values_df)
  }
  
  #1.2) Distance over expression matrix
  if ("corrected" %in% names(assays(object))){
    #expression matrix after correction
    space <- assay(object, "corrected")
    
    similarity_batch_list <- compute_similarity(subset_list = subset_space(space = space, b_vector = batch_vector, N_b = N_batches))
    similarity_cell_type_list <- compute_similarity(subset_list = subset_space(space = space, b_vector = cell_type_vector, N_b = N_cell_types))
    similarity_list <- append(similarity_batch_list, similarity_cell_type_list)
    
    
    values_list <- extract_values(similarity_list)
    
    #convert list of lists to dataframe
    values_df <- as.data.frame(t(plyr::ldply(values_list, rbind, .id = NULL)))
    
    colnames(values_df) <- c(paste0("Counts_dist_BATCH_aft_", names(table(batch_vector))), 
                             paste0("Counts_dist_CELL_aft_", names(table(cell_type_vector))))
    print(summary(values_df))
    print("Distances over the counts matrix succesful!")
    #save results
    save_results(values_df)
  }
}
#2) Seurat objects
if (class(object) == "seurat"){
  library(Seurat)
  #vectors for similarity calculation
  batch_vector <- as.character(object@meta.data[["Batch"]])
  N_batches <- length(unique(batch_vector))
  cell_type_vector <- as.character(object@meta.data[["cell_type1"]])
  N_celltypes <- length(unique(cell_type_vector))

    # 2.1) Seurat v2 multiCCA
    if( "cca" %in% names(object@dr)){
        print("Object corrected by Seurat_v2_multiCCA")
  
  #embedding before correction
  cca_embedding <- attr(object@dr[["cca"]], "cell.embeddings") 
  cca_embedding <- t(cca_embedding)
  #embedding after correction
  al_cca_embedding <- attr(object@dr[["cca.aligned"]], "cell.embeddings") 
  al_cca_embedding <- t(al_cca_embedding)
  
  #join before and after embeddings
  embeddings <- list(cca_embedding, al_cca_embedding)
  b_vectors <- list(batch_vector, cell_type_vector)
  n_bs <- list(N_batches, N_celltypes)
  
  #separate the PCA space by BATCHES and CELL TYPES from both embeddings cca and al_cca
  similarity_list <- list()
  for(i in 1:length(embeddings)){
    for (n in 1:length(b_vectors)){
      similarity <- compute_similarity(subset_list = subset_space(space = embeddings[[i]] , b_vector = b_vectors[[n]], N_b = n_bs[[n]]))
      similarity_list <- append(similarity_list, similarity)
    }
  }
  
  
  values_list <- extract_values(similarity_list)
  #convert list of lists to dataframe
  values_df <- as.data.frame(t(plyr::ldply(values_list, rbind, .id = NULL)))
  
  colnames(values_df) <- c(paste0("bef_cca_dist_BATCH_aft_", names(table(batch_vector))), 
                           paste0("bef_cca_dist_CELL_aft_", names(table(cell_type_vector))), 
                           paste0("aft_cca_dist_BATCH_aft_", names(table(batch_vector))),
                           paste0("aft_cca_dist_CELL_aft_", names(table(cell_type_vector)))
                           )
  print(summary(values_df))
  print("Distances over the Seurat_v2_multiCCA graph succesful")

  #save results
  save_results(values_df)
} 
}

#2) Seurat objects
if (class(object) == "Seurat"){
  library(Seurat)
  #vectors for similarity calculation
  batch_vector <- as.character(object@meta.data[["Batch"]])
  N_batches <- length(unique(batch_vector))
  cell_type_vector <- as.character(object@meta.data[["cell_type1"]])
  N_cell_types <- length(unique(cell_type_vector))

    # 2.1) Seurat v3 anchors
    if( "integrated" %in% names(object@assays)){
    print("Object corrected by Seurat_v3_anchors")

    space <- as.matrix(object@assays[["integrated"]]@data)
    
    similarity_batch_list <- compute_similarity(subset_list = subset_space(space = space, b_vector = batch_vector, N_b = N_batches))
    similarity_cell_type_list <- compute_similarity(subset_list = subset_space(space = space, b_vector = cell_type_vector, N_b = N_cell_types))
    similarity_list <- append(similarity_batch_list, similarity_cell_type_list)
    
    
    values_list <- extract_values(similarity_list)
    
    #convert list of lists to dataframe
    values_df <- as.data.frame(t(plyr::ldply(values_list, rbind, .id = NULL)))
    
    colnames(values_df) <- c(paste0("Counts_dist_BATCH_aft_", names(table(batch_vector))), 
                             paste0("Counts_dist_CELL_aft_", names(table(cell_type_vector))))
    print(summary(values_df))
    print("Distances over the counts matrix succesful!")
    #save results
    save_results(values_df)
    print("Distances over the Seurat_v3_anchors graph succesful") 
 }

