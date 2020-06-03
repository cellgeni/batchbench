#!/usr/bin/env Rscript

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
        c("-a", "--assay_name"),
        action = "store",
        default = "logcounts",
        type = 'character',
        help = 'Counts assay to add to the h5ad object'
    ),
    make_option(
        c("-c", "--corrected_assay"),
        action = "store",
        default = "corrected",
        type = 'character',
        help = 'Corrected counts assay name'
    ),
    make_option(
        c("-e", "--corrected_emb"),
        action = "store",
        default = "corrected_emb",
        type = 'character',
        help = 'Name of the fastMNN corrected low-dimensional embedding.'
    ),
    make_option(
        c("-o", "--output_object"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to h5ad output file'
   ) 
)

suppressPackageStartupMessages(library(reticulate)) 
suppressPackageStartupMessages(library(SingleCellExperiment)) 
suppressPackageStartupMessages(library(loomR))
suppressPackageStartupMessages(library(sceasy)

# args 
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay 
corrected_emb <- opt$corrected_emb
# read input file
dataset <- readRDS(opt$input_object)

if (class(dataset) == "SingleCellExperiment"){
  #Those methods that return an embedding, don't have corrected as an assay
  if(corrected_emb %in% reducedDimNames(dataset)){
  sceasy:::sce2anndata(obj = dataset, outFile = opt$output_object, main_layer = assay_name, transfer_layers = NULL)
  }
  #Methods that correct the expression matrix, do have corrected as an assay
  else {
  sceasy:::sce2anndata(obj = dataset, outFile = opt$output_object, main_layer = corrected_assay, transfer_layers = NULL)
  }
}

if (class(dataset) == "Seurat"){
  sceasy:::seurat2anndata(obj = dataset, outFile = opt$output, assay = corrected_assay, main_layer = 'data', transfer_layers = NULL )
}
