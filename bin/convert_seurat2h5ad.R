#!/usr/bin/env Rscript

#Convert Seurat to h5ad object

suppressPackageStartupMessages(library(reticulate))  
suppressPackageStartupMessages(library(SummarizedExperiment)) 
suppressPackageStartupMessages(library(SingleCellExperiment)) 

suppressPackageStartupMessages(library(loomR))
suppressPackageStartupMessages(library(sceasy))

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {stop("Provide valid path to input SCE object.")}
if (is.null(args[["assay_name"]])) {print("No assay name provided, usign 'logcounts' as default.")
  args[["assay_name"]] <- 'RNA'}
if (is.null(args[["output"]])) {stop("Provide valid path to output h5ad object.")}

#read object
sce <- readRDS(args[["input"]])
#convert sce2h5ad
sceasy:::seurat2anndata(obj = sce, outFile = args[["output"]], main_layer = 'data', assay = args[["assay_name"]], transfer_layers = NULL)

print("Seurat successfully converted to H5ad object.")


