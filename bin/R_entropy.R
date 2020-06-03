#!/usr/bin/env Rscript

# Calculate Shannon entropy over an .rds object, AFTER batch correction. 
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
        c("-c", "--celltype_key"),
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

suppressPackageStartupMessages(require(scran))

#args <-R.utils::commandArgs(asValues=TRUE)
#if (is.null(args[["input"]])) {
#  print("Provide a valid input file name (Batch corrected object) --> RDS file")
#}
#if (is.null(args[["output"]])) {
#  print("Provide a valid output file name (Per batch and cell type entropy) --> CSV file")
#}
#if (is.null(args[["output_matrix"]])) {
#  print("Provide a valid output file name for the matrix. RDS file")
#}

# args
assay_name <- opt$assay_name
batch_key <- opt$batch_key
celltype_key <- opt$celltype_key
k_num = opt$k_num
dim_num = opt$dim_num
corrected_assay <- opt$corrected_assay
corrected_emb <- opt$corrected_emb
#input file
object <- readRDS(opt$input_object)

# entropy function
shannon_entropy <- function(x, batch_vector, N_batches) {
  freq_batch = table(batch_vector[x ==1])/length(batch_vector[x == 1])
  freq_batch_positive = freq_batch[freq_batch > 0]
  return(-sum(freq_batch_positive * log(freq_batch_positive))/log(N_batches))
}

compute_entropy <- function(corrected_space, K_num, dim_num, bool, x, batch_vector, N_batches, cell_type_vector, N_cell_types){
  knn_graph <- as.matrix(buildKNNGraph(corrected_space, k=k_num, d=dim_num, transposed = bool)[])
  batch_entropy <- apply(knn_graph, 1, function(x) {shannon_entropy (x, batch_vector, N_batches)})
  celltype_entropy <- apply(knn_graph, 1, function(x) {shannon_entropy (x, cell_type_vector, N_cell_types)})
  entropy <- cbind(batch_entropy,celltype_entropy)
  names(entropy) <- c("Batch_entropy", "Cell_type_entropy")
  return(entropy)
}

save_results <- function(x, col_names){
  write.csv(x, file = args[["output"]], sep = "", row.names = FALSE, col.names = col_names)
}

# 1) SCE objects 
if (class(object) == "SingleCellExperiment"){
  print("The input object is a Single Cell Experiment class object")
  
  #Calculate entropy
  batch_vector <- as.character(object[[batch_key]])
  N_batches <- length(unique(batch_vector))
  cell_type_vector <- colData(object)[[celltype_key]]
  N_cell_types <- length(unique(cell_type_vector))
  
  # COMPUTE ENTROPY 
  # 1.1) tools that correct low_D embeddings
  if (corrected_emb %in% reducedDimNames(object)){
    corrected_space <- as.matrix(object@reducedDims@listData[[corrected_emb]])
    col_names <- c("PCA_batch_entropy", "PCA_cell_type_entropy")
    save_results(compute_entropy(corrected_space, k_num = k_num, dim_num = dim_num, bool = TRUE, x, batch_vector, N_batches, cell_type_vector, N_cell_types), col_names) 
    print("Entropy calculated in PCA space!")
    }
  
  # 1.2) methods that correct expression matrix
  if ("corrected" %in% names(assays(object))){
    corrected_space <- assay(object, "corrected")
    col_names <- c("Batch_entropy", "Cell_type_entropy")
    save_results(compute_entropy(corrected_space, k_num = k_num, dim_num = dim_num, bool = FALSE, x, batch_vector, N_batches, cell_type_vector, N_cell_types), col_names) 
    print("Entropy calculated over counts matrix!")
  }
}

#2) Seurat objects
if (class(object) == "Seurat"){
  suppressPackageStartupMessages(require(Seurat))
  print("The input object is a Seurat class object")
  
  batch_vector <- object@meta.data[[batch_key]]
  N_batches <- length(unique(batch_vector))
  cell_type_vector <- object@meta.data[[celltype_key]]
  N_cell_types <- length(unique(cell_type_vector))

    space <- as.matrix(object@assays[[corrected_assay]]@data)
    col_names <- c("Batch_entropy", "Cell_type_entropy")
    save_results(compute_entropy(corrected_space = space, bool = FALSE, x, batch_vector, N_batches, cell_type_vector, N_cell_types), col_names) 
    print("Entropy calculated over Seurat object!")
}
