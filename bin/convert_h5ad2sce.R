#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to h5ad object'
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
  list(corrected = counts_mat)
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
  var_names <- as.vector(h5ad_file$var_names$values)
  var_names
}
# set matrix dimnames
set_dimnames <- function(mat, cell_names){
     dimnames(mat) <- list(cell_names, cell_names)
     mat
}
# function to extract BBKNN graph
extract_bbknn_graphs <- function(h5ad_file){
  connect_graph <- as(reticulate::py_to_r(h5ad_file[["obsp"]][["connectivities"]]$toarray()), "dgCMatrix")
  dist_graph <- as(reticulate::py_to_r(h5ad_file[["obsp"]][["distances"]]$toarray()), "dgCMatrix")
  graph_list <- list(connectivities = connect_graph, distances = dist_graph)
  cell_names <- extract_obs_names(h5ad_file)
  # set matrix dimnames
  graph_list <- lapply(graph_list, function(mat) set_dimnames(mat, cell_names = cell_names))
  graph_list
}

# build SCE object from H5ad object parts
build_sce <- function(assay_list, meta_data, gene_names, emb_list){
  sce <- SingleCellExperiment(assays = assay_list, 
                              colData = meta_data,
                              reducedDims = emb_list)
  # if present , add rownames
  if(!is.null(gene_names)){rowData(sce) <- data.frame("feature_names" = gene_names)}
  sce
}

## EXECUTE ##

# args 
method <- opt$method
if(is.na(method) || is.null(method)) stop("Please specify the batch correction method the input object comes from")

suppressPackageStartupMessages(require(reticulate))
suppressPackageStartupMessages(require(SingleCellExperiment))

anndata <- reticulate::import('anndata', convert = F)
# Read h5ad file
h5ad_file <- anndata$read_h5ad(opt$input_object)
# load additional libraries
pandas <- reticulate::import('pandas', convert = F)
numpy <- reticulate::import('numpy', convert = F)

# Scanorama h5ad object conversion
if(method %in% c("scanorama", "Scanorama")){
  # extract counts
  counts_mat_list <- extract_counts(h5ad_file)
  # extract metadata
  meta_data <- extract_metadata(h5ad_file)
  # extract embeddings
  emb_list <- extract_embeddings(h5ad_file)
  # extract gene names
  gene_names <- extract_var_names(h5ad_file)
  # build SCE
  sce <- build_sce(assay_list = counts_mat_list, meta_data = meta_data, gene_names = gene_names, emb_list = emb_list)
}

# BBKNN h5ad object conversion
if(method %in% c("bbknn", "BBKNN")){
  # extract graph
  graph_list <- extract_bbknn_graphs(h5ad_file)
  # extract metadata
  meta_data <- extract_metadata(h5ad_file)
  # extract embeddings
  emb_list <- extract_embeddings(h5ad_file)
  # extract gene names
  cell_names <- extract_var_names(h5ad_file)
  # build SCE
  sce <- build_sce(assay_list = graph_list, meta_data = meta_data, gene_names = NULL, emb_list = emb_list)
}

# save converted object
saveRDS(sce, opt$output_object)
print(paste0(method, " h5ad object succesfully converted to SCE"))
