#!/usr/bin/env Rscript

#This script is to calculate shannon entropy over an .rds object, before batch correction. 
#To know the degree of mixing of the dataprior to correction. 
#We calculate both batch_entropy and cell_type entropy over the PCA graph and the count_matrix.

#libraries
library(scran)
#value arguments
args <-R.utils::commandArgs(asValues=TRUE)
if (is.null(args[["input"]])) {
  print("Provide a valid input file name (Batch corrected object) --> RDS file")
}
if (is.null(args[["output"]])) {
  print("Provide a valid output file name (Per batch and cell type entropy) --> CSV file")
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
  #names(entropy) <- c("batch_entropy", "cell_type_entropy")
  return(entropy)
}

save_results <- function(x){
  write.table(x, file = args[["output"]], sep = "\t",  row.names = FALSE)
}

#arguments for entropy calculation
batch_vector <- colData(object)$Batch
N_batches <- length(unique(batch_vector))

if ("cell_type1" %in% colnames(colData(object))){
  print("The object has cell type annotation")
  cell_type_vector <- colData(object)$cell_type1
  N_cell_types <- length(unique(cell_type_vector))
}
  
#entropy over the PCA graph 
library(scater)
object <- runPCA(object, method = "prcomp", exprs_values = "logcounts", ncomponents = 10)

corrected_space <- as.matrix(object@reducedDims@listData[["PCA"]])
entropy_pca <- compute_entropy(corrected_space, bool = TRUE, x, batch_vector, N_batches, cell_type_vector, N_cell_types) # <------ FUNCTION
print(class(entropy_pca))
colnames(entropy_pca) <- c("Bef_correct_PCA_batch", "Bef_correct_PCA_cell_type")
print("Entropy calculated in PCA space!")

#entropy over the expression matrix
corrected_space <- assay(object, "logcounts")
entropy_counts <- compute_entropy(corrected_space, bool = FALSE, x, batch_vector, N_batches, cell_type_vector, N_cell_types) # <------ FUNCTION
print(class(entropy_counts))
colnames(entropy_counts) <- c("Bef_correct_Counts_batch", "Bef_correct_Counts_cell_type")
print("Entropy calculated over counts matrix!")
    
result <- cbind(entropy_pca,entropy_counts)

print(summary(result))
#write result to a csv
save_results(result)