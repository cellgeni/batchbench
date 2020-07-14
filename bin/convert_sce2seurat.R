#!/usr/bin/env Rscript

# Convert SCE objects to Seurat.

# TODO
## Specify which low-D embedding to incorporate in the Seurat object. --> For Clustering analysis!

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
    c("-a", "--assay_name"),
    action = "store",
    default = "logcounts",
    type = 'character',
    help = 'Counts assay name' 
  ),
  make_option(
    c("-m", "--method"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Bacth correction method the input comes from'
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

suppressPackageStartupMessages(require(Seurat))

# args
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay
method <- opt$method
# read input object
sce <- readRDS(opt$input_object)
# convert sce to seurat
if(method == "logcounts"){ sce2seurat <- as.Seurat(sce, assay = assay_name, data = assay_name, counts = assay_name) }else{ sce2seurat <- as.Seurat(sce, assay = corrected_assay, data = corrected_assay, counts = corrected_assay)} 
# save object
saveRDS(sce2seurat, file = opt$output_object)
print("SCE successfully converted to Seurat object!")
