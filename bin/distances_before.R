#!/usr/bin/env Rscript

#This script is for the calculation of the cosine distances of the datasets prior to batch correction.
#As some methods correct the PCA space and others the count matrix, distance over both spaces are calculated.i
#libraries
library(SingleCellExperiment)
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
 # write.table(x, file = args[["output"]], sep = "\t", row.names = FALSE)
  saveRDS(x, args[["output"]])
}

#arguments for similarity calculation
batch_vector <- as.character(object$Batch)
N_batches <- length(unique(batch_vector))

if ("cell_type1" %in% colnames(colData(object))){
  print("The object has cell type annotation")
  cell_type_vector <- as.character(object$cell_type1)
  N_cell_types <- length(unique(cell_type_vector))
}

#1.1) Distance over PCA graph
library(scater)
object <- runPCA(object, ncomponents = 10, method = "prcomp", exprs_values = "logcounts")

space <- object@reducedDims@listData[["PCA"]]

space <- t(space) #We transpose to avoid code duplication, as pca_graph is cells x PCs and exp_matrix is genes x cells. 
similarity_batch_list <- compute_similarity(subset_list = subset_space(space = space, b_vector = batch_vector, N_b = N_batches))
similarity_cell_type_list <- compute_similarity(subset_list = subset_space(space = space, b_vector = cell_type_vector, N_b = N_cell_types))
pca_similarity_list <- append(similarity_batch_list, similarity_cell_type_list)
#extract values
pca_values_list <- extract_values(pca_similarity_list)
#name the list
names(pca_values_list)<- c(paste0("bef_correct_PCA_dist_BATCH_", names(table(batch_vector))),
                         paste0("bef_correct_PCA_dist_CELL_", names(table(cell_type_vector))))


#1.2) Distance over counts matrix
space <- assay(object, "logcounts")

similarity_batch_list <- compute_similarity(subset_list = subset_space(space = space, b_vector = batch_vector, N_b = N_batches))
similarity_cell_type_list <- compute_similarity(subset_list = subset_space(space = space, b_vector = cell_type_vector, N_b = N_cell_types))
counts_similarity_list <- append(similarity_batch_list, similarity_cell_type_list)
#extract values
counts_values_list <- extract_values(pca_similarity_list)
#name the list
names(counts_values_list)<- c(paste0("bef_correct_Counts_dist_BATCH_", names(table(batch_vector))),
                         paste0("bef_correct_Counts_dist_CELL_", names(table(cell_type_vector))))

                                                                           
results_list <- list(pca_values_list, counts_values_list)


print("Distances over the Pca graph, and counts matrix (prior to batch correction) succesful!")
#save results
save_results(results_list)
