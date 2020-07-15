#!/usr/bin/env Rscript

# Calculate Shannon entropy over an .rds object 
# Entropy is calculated depending on the corrected space that enters the script.

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to rds input file'
  ),
  make_option(
    c("-a", "--assay_name"),
    action = "store",
    default = "logcounts",
    type = 'character',
    help = 'Counts assay to add to the h5ad object'
  ),
  make_option(
    c("-c", "--corrected_assay"),
    action = "store",
    default = "corrected",
    type = 'character',
    help = 'Corrected counts assay name'
  ),
  make_option(
    c("-m", "--method"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Bacth correction method the input comes from'
  ),
  make_option(
    c("-e", "--corrected_emb"),
    action = "store",
    default = "corrected_emb",
    type = 'character',
    help = 'Name of the fastMNN corrected low-dimensional embedding.'
  ),
  make_option(
    c("-b", "--batch_key"),
    action = "store",
    default = "Batch",
    type = 'character',
    help = 'Batch key in cell metadata'
  ),
  make_option(
    c("-t", "--celltype_key"),
    action = "store",
    default = "cell_type1",
    type = 'character',
    help = 'Cell type key in cell metadata'
  ),
  make_option(
    c("-k", "--k_num"),
    action = "store",
    default = "30",
    type = 'integer',
    help = 'number of nearest neighbors to consider during graph construction.'
  ),
  make_option(
    c("-d", "--dim_num"),
    action = "store",
    default = "50",
    type = 'integer',
    help = ' number of dimensions to use for the search'
  ),
  make_option(
    c("-o", "--output_entropy"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to csv output file'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

## FUNCTIONS ##
# entropy function
shannon_entropy <- function(x, batch_vector, N_batches) {
  freq_batch = table(batch_vector[x >=1])/length(batch_vector[x >= 1])
  freq_batch_positive = freq_batch[freq_batch > 0]
  return(-sum(freq_batch_positive * log(freq_batch_positive))/log(N_batches))
}

# build KNN graph 
build_knn_graph <- function(corrected_space, K_num, dim_num, bool){
  knn_graph <- as.matrix(buildKNNGraph(corrected_space, k=k_num, d=dim_num, transposed = bool)[])
  knn_graph
}

# compute shannon entropy per cell
compute_entropy <- function(knn_graph, batch_vector, N_batches, cell_type_vector, N_cell_types){
  batch_entropy <- apply(knn_graph, 1, function(x) {shannon_entropy (x, batch_vector, N_batches)})
  celltype_entropy <- apply(knn_graph, 1, function(x) {shannon_entropy (x, cell_type_vector, N_cell_types)})
  entropy <- cbind(batch_entropy,celltype_entropy)
  names(entropy) <- c("Batch_entropy", "Cell_type_entropy")
  entropy
}
# save entropy values per cell
save_results <- function(entropy, row_names){
  write.csv(entropy, file = opt$output_entropy, row.names = row_names)
}

##  EXECUTE  ##

# args
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay
corrected_emb <- opt$corrected_emb
method <- opt$method
if(is.null(method) || is.na(method)){ stop("Please provide the batch correction method the input file comes from") }
batch_key <- opt$batch_key
celltype_key <- opt$celltype_key
k_num = opt$k_num
dim_num = opt$dim_num

# input file
dataset <- readRDS(opt$input_object)

batch_vector <- as.character(dataset$Batch)
N_batches <- length(unique(batch_vector))
cell_type_vector <- as.character(dataset@colData[[celltype_key]])
N_cell_types <- length(unique(cell_type_vector))

print(batch_vector)
print(N_batches)
print(cell_type_vector)
print(N_cell_types)
# load scran package to compute KNN graph
suppressPackageStartupMessages(require(scran))

# for graph correcting methods (BBKNN)
if (method %in% c("bbknn", "BBKNN")){
  corrected_graph <- assay(dataset, corrected_assay) # graph stored in assay in SCE for BBKNN
  entropy <- compute_entropy(knn_graph = corrected_graph, batch_vector, N_batches, cell_type_vector, N_cell_types)
  save_results(entropy, row_names = rownames(corrected_graph))
}

# for embedding correcting methods (fastMNN, harmony)
if (method %in% c("fastMNN", "FastMNN", "Harmony", "harmony")){
  corr_embedding <- reducedDim(dataset, corrected_emb)
  print(corr_embedding)
  knn_graph <- build_knn_graph(corrected_space = corr_embedding, K_num = k_num, dim_num = dim_num, bool = TRUE)
  entropy <- compute_entropy(knn_graph = knn_graph, batch_vector, N_batches, cell_type_vector, N_cell_types)
  save_results(entropy, row_names = rownames(corr_embedding))
}

# for counts matrix correcting methods except Seurat (logcounts, mnnCorrect, limma, ComBat, Scanorama)
if (method %in% c("logcounts", "Logcounts", "mnnCorrect", "mnncorrect", "limma", "Limma", "ComBat", "combat", "Scanorama", "scanorama", "Seurat3", "seurat3", "Seurat")){
  if(method %in% c("logcounts", "Logcounts")){
    counts_mat <- assay(dataset, assay_name)
  } else { counts_mat <- assay(dataset, corrected_assay)
  }
  knn_graph <- build_knn_graph(corrected_space = counts_mat, K_num = k_num, dim_num = dim_num, bool = FALSE)
  entropy <- compute_entropy(knn_graph = knn_graph, batch_vector, N_batches, cell_type_vector, N_cell_types)
  save_results(entropy, row_names = colnames(counts_mat))
}

# for Seurat object
#if (method %in% c("Seurat3", "seurat3", "Seurat")){
#  suppressPackageStartupMessages(require(Seurat))
#  
#  batch_vector <- dataset@meta.data[[batch_key]]
#  N_batches <- length(unique(batch_vector))
#  cell_type_vector <- dataset@meta.data[[celltype_key]]
#  N_cell_types <- length(unique(cell_type_vector))
#  
#  counts_mat <- dataset@assays[[corrected_assay]]@data
#  knn_graph <- build_knn_graph(corrected_space = counts_mat, K_num = k_num, dim_num = dim_num, bool = FALSE)
#  entropy <- compute_entropy(knn_graph = knn_graph, batch_vector, N_batches, cell_type_vector, N_cell_types)
#  save_results(entropy, row_names = colnames(counts_mat))
#}
#
print(paste0("Congrats! Entropy computed for method: ", method))
