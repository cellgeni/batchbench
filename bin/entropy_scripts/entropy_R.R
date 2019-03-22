#libraries
library(scran)

args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}  

#input file
object <- readRDS(args[1])
#space over which entropy has to be calculated
space <- args[2]

#entropy function
shannon_entropy <- function(x, batch_vector, N_batches) {
  freq_batch = table(batch_vector[x ==1])/length(batch_vector[x == 1])
  freq_batch_positive = freq_batch[freq_batch > 0]
  return(-sum(freq_batch_positive * log(freq_batch_positive))/log(N_batches))
}

# 1) SCE objects 
if (class(object) == "SingleCellExperiment"){
  
  library(SingleCellExperiment) 
  print("The input object is a Single Cell Experiment class object")
  
  #Calculate entropy
  batch_vector <- colData(object)$dataset
  N_batches <- length(unique(batch_vector))
  print("batch_vector Done!")

  
  if ("cell_type1" %in% colnames(colData(object))){
      print("The object has cell type annotation")
      cell_type_vector <- colData(object)$cell_type1          
      N_cell_types <- length(unique(cell_type_vector))
  }
  
  # COMPUTE ENTROPY 
  #methods that correct expression matrix
  if (space == "exp_mat"){
      print("Entropy is being calculated over the expression matrix")
      #entropy before 
      exp_matrix <- assay(object, "logcounts")
      knn_before <- as.matrix(buildKNNGraph(exp_matrix, k=30, d=50, transposed=FALSE)[])
      glob_entropy_bef <- apply(knn_before, 1, function(x) {shannon_entropy (x, batch_vector, N_batches)})
      celltype_entropy_before <- apply(knn_before, 1, function(x) {shannon_entropy (x, cell_type_vector, N_cell_types)})
      
      #entropy after 
      corrected_matrix <- assay(object, "corrected")
      knn_graph_after <- as.matrix(buildKNNGraph(corrected_matrix, k=30, d=50, transposed=FALSE)[])
      glob_entropy_after <- apply(knn_graph_after, 1, function(x) {shannon_entropy (x, batch_vector, N_batches)})
      celltype_entropy_after <- apply(knn_graph_after, 1, function(x) {shannon_entropy (x, cell_type_vector, N_cell_types)})
      
      print("Entropy calculated!")
      
      result <- cbind(glob_entropy_bef,celltype_entropy_before, glob_entropy_after, celltype_entropy_after )
      print(summary(result))
      #write result to a csv?
      write.csv(result, file = args[3])
  }
  
  #tools that correct low_D embeddings
  if (space == "pca_graph"){
      print("Entropy is being calculated over the low-D embedding")
      #entropy before 
      PCA_matrix <- as.matrix(object@reducedDims@listData[["PCA"]])
      knn_before <- as.matrix(buildKNNGraph(PCA_matrix, k=30, d=50, transposed=TRUE)[])
      glob_entropy_bef <- apply(knn_before, 1, function(x) {shannon_entropy (x, batch_vector, N_batches)})
      celltype_entropy_before <- apply(knn_before, 1, function(x) {shannon_entropy (x, cell_type_vector, N_cell_types)})
        
      #entropy after 
      corrected_embedding <- as.matrix(object@reducedDims@listData[["corrected_embedding"]])
      kNN_after <- as.matrix(buildKNNGraph(PCA_matrix, k=30, d=50, transposed=TRUE)[])
      glob_entropy_after <- apply(kNN_after, 1, function(x) {shannon_entropy (x, batch_vector, N_batches)})
      celltype_entropy_after <- apply(kNN_after, 1, function(x) {shannon_entropy (x, cell_type_vector, N_cell_types)})
      
      print("Entropy calculated!")
      
      #write result to a csv?
      result <- cbind(glob_entropy_bef,celltype_entropy_before, glob_entropy_after, celltype_entropy_after)
      print(summary(result))
      write.csv(result, file = args[3])
  }
}

#2) Seurat objects
if (class(object) == "seurat"){
  library(Seurat)
  print("The input object is a Seurat class object")
  
  batch_vector <- object@meta.data$dataset
  N_batches <- length(unique(batch_vector))
  
  if ("cell_type1" %in% colnames(object@meta.data)){
    cell_type_vector <- object@meta.data$cell_type1           
    N_cell_types <- length(unique(cell_type_vector))
  }
  
  #entropy before
  CCA_matrix <- GetCellEmbeddings(object = object, reduction.type = "cca")
  knn_before <- as.matrix(buildKNNGraph(CCA_matrix, k=30, d=50, transposed=TRUE)[])
  glob_entropy_bef <- apply(knn_before, 1, function(x) {shannon_entropy (x, batch_vector, N_batches)})
  
  #entropy after
  Aligned_CCA_matrix <- GetCellEmbeddings(object = object, reduction.type = "cca.aligned")
  knn_after <- as.matrix(buildKNNGraph(Aligned_CCA_matrix, k=30, d=50, transposed=TRUE)[])
  glob_entropy_after <- apply(knn_after, 1, function(x) {shannon_entropy (x, batch_vector, N_batches)})
  
  #if there is cell type annotation, calculate cell type entropy
  if(exists(x = "cell_type_vector") & exists(x = "N_cell_types")){
  celltype_entropy_before <- apply(knn_before, 1, function(x) {shannon_entropy (x, cell_type_vector, N_cell_types)})
  celltype_entropy_after <- apply(knn_after, 1, function(x) {shannon_entropy (x, cell_type_vector, N_cell_types)})                            
  }
  
  print("Entropy calculated!")
  
  if (exists(x = "celltype_entropy_before") & exists(x = "celltype_entropy_after")){
  result <- cbind(glob_entropy_bef,celltype_entropy_before, glob_entropy_after, celltype_entropy_after)
  } else {
    result <- cbind(glob_entropy_bef, glob_entropy_after)
    }
  
  print(summary(result))
  write.csv(result, file = args[3])

  }