#!/usr/bin/env Rscript

# Convert h5ad objects to Seurat.

suppressPackageStartupMessages(library("optparse"))

option_list = list(
    make_option(
        c("-i", "--input_object"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to h5ad input file' 
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
        help = 'Path to rdd Seurat class output file'
    )
)
opt <- parse_args(OptionParser(option_list=option_list))


suppressPackageStartupMessages(require(reticulate))
suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(Seurat))

# args 
corrected_assay <- opt$corrected_assay

# function to manually convert h5ad to Seurat 
conv_h5ad2seurat <- function(h5ad_path){
  anndata <- reticulate::import('anndata', convert = F)
  pandas <- reticulate::import('pandas', convert = F)
  numpy <- reticulate::import('numpy', convert = F) 
  # read h5ad file
  h5ad_file <- anndata$read_h5ad(h5ad_path)
  # Extract Counts   
  # set up condition ERROR this is because logcounts and scanorama counts are numpy objects with different properties
  h5ad_counts <- tryCatch({
    pandas$DataFrame(h5ad_file$X, index = h5ad_file$obs_names, columns = h5ad_file$var_names)
  }, error=function(e){pandas$DataFrame(h5ad_file$X$todense(), index = h5ad_file$obs_names, columns = h5ad_file$var_names)})
  
  seurat_counts <- t(as.matrix(reticulate::py_to_r(h5ad_counts)))
  # Extract Metadata
  meta_data <- reticulate::py_to_r(pandas$DataFrame(h5ad_file$obs, dtype = "object"))
  # Create seurat object
  seurat_h5ad <- CreateSeuratObject(counts = seurat_counts, meta.data = meta_data)
}

# convert with seurat function
h5ad2seurat <- tryCatch({
  ReadH5AD(opt$input_object, assay = corrected_assay)
}, error=function(e){conv_h5ad2seurat(args[["input"]])})

# save seurat object
# substitute file name suffix
rds_file_name <- gsub(".h5ad", ".rds", opt$output_object)
saveRDS(h5ad2seurat, file = rds_file_name)
print("h5ad successfully converted to Seurat object!")
