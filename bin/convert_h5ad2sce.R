#!/usr/bin/env Rscript

# Convert SCE to H5ad class object
#suppressPackageStartupMessages(library("optparse"))
#
#option_list = list(
#    make_option(
#        c("-i", "--input_object"),
#        action = "store",
#        default = NA,
#        type = 'character',
#        help = 'Path to rds input file' 
#    ),
#    make_option(
#        c("-b", "--batch_key"),
#        action = "store",
#        default = "Batch",
#        type = 'character',
#        help = 'Minimum number of cells for a gene to be expressed in.'
#    ),
#    make_option(
#        c("-c", "--celltype_key"),
#        action = "store",
#
#Convert h5ad objects to SCE.

#TODO
##Fix assayname assignment via assay_name argument!
##NOTE: Currently assay saved as 'corrected'.

suppressPackageStartupMessages(require(reticulate))
suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(Seurat))

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {stop("Provide valid path to input SCE object.")}
#if (is.null(args[["assay_name"]])) {stop("Provide an assay name for to store h5ad counts in the SCE object.")}
if (is.null(args[["output"]])) {stop("Provide valid path to output Seurat object.")}

#function to manually convert h5ad to Seurat 
h5ad2sce <- function(h5ad_path){
  anndata <- reticulate::import('anndata', convert = F)
  pandas <- reticulate::import('pandas', convert = F)
  numpy <- reticulate::import('numpy', convert = F) 
  #Read h5ad file
  h5ad_file <- anndata$read_h5ad(h5ad_path)
  #Extract Counts   
  #set up condition ERROR this is because logcounts and scanorama counts are numpy objects with different properties
  h5ad_counts <- tryCatch({
    pandas$DataFrame(h5ad_file$X, index = h5ad_file$obs_names, columns = h5ad_file$var_names)
  }, error=function(e){pandas$DataFrame(h5ad_file$X$todense(), index = h5ad_file$obs_names, columns = h5ad_file$var_names)})
  
  counts_mat <- t(as.matrix(reticulate::py_to_r(h5ad_counts)))
  #Extract Metadata
  meta_data <- reticulate::py_to_r(pandas$DataFrame(h5ad_file$obs, dtype = "object"))
  #Extract Embeddings
  emb_names <- reticulate::py_to_r(h5ad_file$obsm_keys())
  emb_list <- lapply(as.list(emb_names), function(x) reticulate::py_to_r(h5ad_file$obsm[x]))
  names(emb_list) <- emb_names

  #Build SCE
  sce <- build_sce(counts_mat = counts_mat, meta_data = meta_data, emb_list = emb_list)
  sce
}


#function to build sce
build_sce <- function(counts_mat, meta_data, assay_name, emb_list){
  sce <- SingleCellExperiment(assays = list(corrected = counts_mat),
                            colData = meta_data,
                            rowData = list('feature_names' = rownames(counts_mat)),
                            reducedDims = emb_list)
  sce
}

sce <- h5ad2sce(h5ad_path = args[["input"]])

#save converted object
saveRDS(sce, file = args[["output"]])
print("h5ad successfully converted to SCE object!")
