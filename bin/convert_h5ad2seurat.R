#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to rds Seurat object'
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
    help = 'Batch correction method the object comes from'
  ),
  make_option(
    c("-o", "--output_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the rds SCE object'
  )
)
opt <- parse_args(OptionParser(option_list=option_list))

## FUNCTIONS ##
# extract counts 
extract_counts <- function(h5ad_file){
  # set up condition error as numpy matrices may vary its properties depending on the dataset
  h5ad_counts <- tryCatch({
    pandas$DataFrame(h5ad_file$X, index = h5ad_file$obs_names, columns = h5ad_file$var_names)
  }, error=function(e){pandas$DataFrame(h5ad_file$X$todense(), index = h5ad_file$obs_names, columns = h5ad_file$var_names)})
  # transpose as h5ad matrices are cells x genes
  counts_mat <- t(as.matrix(reticulate::py_to_r(h5ad_counts)))
  # add dimnames
  gene_names <- extract_var_names(h5ad_file)
  cell_names <- extract_obs_names(h5ad_file)
  dimnames(counts_mat) <- list(gene_names, cell_names)
  counts_mat
}

# extract metadata
extract_metadata <- function(h5ad_file){
  meta_data <- reticulate::py_to_r(pandas$DataFrame(h5ad_file$obs, dtype = "object"))
  meta_data
}

# extract embeddings
extract_embeddings <- function(h5ad_file){
  emb_names <- reticulate::py_to_r(h5ad_file$obsm_keys())
  emb_list <- lapply(as.list(emb_names), function(x) reticulate::py_to_r(h5ad_file$obsm[x]))
  names(emb_list) <- emb_names
  emb_list
}

# function to extract obs_names
extract_obs_names <- function(h5ad_file){
  obs_names <- as.vector(h5ad_file$obs_names$values)
  obs_names
}
# function to extract var_names
extract_var_names <- function(h5ad_file){
  obs_names <- as.vector(h5ad_file$var_names$values)
  obs_names
}

# function to extract BBKNN graph
extract_bbknn_graph <- function(h5ad_file){
  graph_conn <- as(reticulate::py_to_r(h5ad_file[["obsp"]][["connectivities"]]$toarray()), "dgCMatrix")
  cell_names <- extract_obs_names(h5ad_file)
  dimnames(graph_conn) <- list(cell_names, cell_names)
  graph_conn
}

# build SCE object from H5ad object parts
build_sce <- function(counts_mat, meta_data, gene_names, emb_list){
  sce <- SingleCellExperiment(assays = list(corrected = counts_mat),
                              colData = meta_data,
                              reducedDims = emb_list)
  # if present , add rownames
  if(!is.null(gene_names)){rowData(sce) <- data.frame("feature_names" = gene_names)}
  sce
}

# build Seurat object
build_seurat <- function(counts_mat, corrected_assay, meta_data, emb_list, add_graph){
  seurat_object <- CreateSeuratObject(counts = counts_mat, assay = corrected_assay, meta.data = meta_data)
  seurat_object@reductions <- emb_list
  # add BBKNN graph 
  if(!(is.null(add_graph))){seurat_object@graphs[["bbknn_connectivities"]] <- add_graph}
  seurat_object
}

## EXECUTE ##

# args 
method <- opt$method
if(is.na(method) || is.null(method)) stop("Please specify the batch correction method the input object comes from")
corrected_assay <- opt$corrected_assay

suppressPackageStartupMessages(require(reticulate))
suppressPackageStartupMessages(require(Seurat))

anndata <- reticulate::import('anndata', convert = F)
# Read h5ad file
h5ad_file <- anndata$read_h5ad(opt$input_object)
# load additional libraries
pandas <- reticulate::import('pandas', convert = F)
numpy <- reticulate::import('numpy', convert = F)

# Scanorama h5ad object conversion
if(method %in% c("scanorama", "Scanorama")){
  # extract counts
  counts_mat <- extract_counts(h5ad_file)
  # extract metadata
  meta_data <- extract_metadata(h5ad_file)
  # extract embeddings
  emb_list <- extract_embeddings(h5ad_file)
  # extract gene names
  gene_names <- extract_var_names(h5ad_file)
  # build SeuratObject
  seurat_obj <- build_seurat(counts_mat = counts_mat, corrected_assay = corrected_assay, meta_data = meta_data, emb_list = emb_list, add_graph = NULL)
  }

# BBKNN h5ad object conversion
if(method %in% c("bbknn", "BBKNN")){
  # extract counts
  counts_mat <- extract_counts(h5ad_file)
  # extract graph
  graph_conn <- extract_bbknn_graph(h5ad_file)
  # extract metadata
  meta_data <- extract_metadata(h5ad_file)
  # extract embeddings
  emb_list <- extract_embeddings(h5ad_file)
  # extract gene names
  cell_names <- extract_var_names(h5ad_file)
  # build SeuratObject
  seurat_obj <- build_seurat(counts_mat = counts_mat, corrected_assay = corrected_assay, meta_data = meta_data, emb_list = emb_list, add_graph = graph_conn)
}

# save converted object
saveRDS(seurat_obj, opt$output_object)
print(paste0(method, " h5ad object succesfully converted to Seurat Object"))
