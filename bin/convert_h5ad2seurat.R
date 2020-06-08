#!/usr/bin/env Rscript

#Convert h5ad objects to Seurat.

#TODO: probably in the future we wont need the manuall conv_h5ad2seurat as we won't have logcounts.h5ad
#Because input will be either a rds or Seurat object and then converted. Anyway that's the future.
  #Add EMbedding
  #Extract Embeddings
  #emb_names <- reticulate::py_to_r(h5ad_file$obsm_keys())
  #emb_list <- lapply(as.list(emb_names), function(x) reticulate::py_to_r(h5ad_file$obsm[x]))
  #names(emb_list) <- emb_names
#NOTE: Output assay name set to 'corrected'

suppressPackageStartupMessages(require(reticulate))
suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(Seurat))

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {stop("Provide valid path to input SCE object.")}
if (is.null(args[["output"]])) {stop("Provide valid path to output Seurat object.")}

#function to manually convert h5ad to Seurat 
conv_h5ad2seurat <- function(h5ad_path){
  anndata <- reticulate::import('anndata', convert = F)
  pandas <- reticulate::import('pandas', convert = F)
  numpy <- reticulate::import('numpy', convert = F) 
  #read h5ad file
  h5ad_file <- anndata$read_h5ad(h5ad_path)
  #Extract Counts   
  #set up condition ERROR this is because logcounts and scanorama counts are numpy objects with different properties
  h5ad_counts <- tryCatch({
    pandas$DataFrame(h5ad_file$X, index = h5ad_file$obs_names, columns = h5ad_file$var_names)
  }, error=function(e){pandas$DataFrame(h5ad_file$X$todense(), index = h5ad_file$obs_names, columns = h5ad_file$var_names)})
  
  seurat_counts <- t(as.matrix(reticulate::py_to_r(h5ad_counts)))
  #Extract Metadata
  meta_data <- reticulate::py_to_r(pandas$DataFrame(h5ad_file$obs, dtype = "object"))
  #Create seurat object
  seurat_h5ad <- CreateSeuratObject(counts = seurat_counts, meta.data = meta_data)
}

#convert with seurat funciton
h5ad2seurat <- tryCatch({
  ReadH5AD(args[["input"]], assay = "corrected")
}, error=function(e){conv_h5ad2seurat(args[["input"]])})

#save seurat object
saveRDS(h5ad2seurat, file = args[["output"]])
print("h5ad successfully converted to Seurat object!")
