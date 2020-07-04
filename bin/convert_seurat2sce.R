#!/usr/bin/env Rscript

# Convert Seurat to SCE

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

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(SingleCellExperiment))
# args
corrected_assay <- opt$corrected_assay
# read Seurat object
seurat_obj <- readRDS(opt$input_object)
if(!(corrected_assay %in% names(seurat_obj@assays))) stop("Corrected assay name provided not in Seurat object")
# Convert Seurat to SCE
seurat2sce <- as.SingleCellExperiment(seurat_obj, assay = corrected_assay)
saveRDS(seurat2sce, opt$output_object)
print("Seurat successfully converted to SCE object!")
