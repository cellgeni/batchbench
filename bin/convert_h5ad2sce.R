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
        c("-o", "--output_object"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to the rds SCE object'
    )
)
opt <- parse_args(OptionParser(option_list=option_list))

# args
corrected_assay <- opt$corrected_assay

suppressPackageStartupMessages(require(reticulate))
suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(Seurat))

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
build_sce <- function(counts_mat, meta_data, corrected_assay, emb_list){
  sce <- SingleCellExperiment(assays = list(corrected_assay = counts_mat),
                            colData = meta_data,
                            rowData = list('feature_names' = rownames(counts_mat)),
                            reducedDims = emb_list)
  sce
}

sce <- h5ad2sce(h5ad_path = opt$input_object)

#save converted object
saveRDS(sce, file = opt$output_object)
print("h5ad successfully converted to SCE object!")
