#!/usr/bin/env Rscript

#Convert SCE to h5ad object

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
        help = 'Path to h5ad output file'
    )
)
opt <- parse_args(OptionParser(option_list=option_list))

#suppressPackageStartupMessages(library(reticulate))  
suppressPackageStartupMessages(library(SingleCellExperiment)) 
suppressPackageStartupMessages(library(loomR))
suppressPackageStartupMessages(library(sceasy))

# read input object
sce <- readRDS(opt$input_object)
# convert sce2h5ad and save
sceasy:::sce2anndata(obj = sce, outFile = opt$output_object, main_layer = opt$corrected_assay, transfer_layers = NULL)
print("SCE successfully converted to H5ad object.")
