#!/usr/bin/env Rscript

#Convert Seurat to SCE

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
        c("-a", "--assay_name"),
        action = "store",
        default = "logcounts",
        type = 'character',
        help = 'Counts assay considered' 
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
suppressPackageStartupMessages(require(sceasy))
# args
assay_name <- opt$assay_name
# read Seurat object
seurat_obj <- readRDS(opt$input_object)
# Convert Seurat to SCE
seurat2sce <- as.SingleCellExperiment(seurat_obj, assay = assay_name)
saveRDS(seurat2sce, opt$output_object)
print("Seurat successfully converted to SCE object!")
