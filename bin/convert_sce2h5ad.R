#!/usr/bin/env Rscript

#Convert SCE to h5ad object

suppressPackageStartupMessages(library(reticulate))  
suppressPackageStartupMessages(library(SummarizedExperiment)) 
suppressPackageStartupMessages(library(SingleCellExperiment)) 

suppressPackageStartupMessages(library(loomR))
suppressPackageStartupMessages(library(sceasy))

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {stop("Provide valid path to input SCE object.")}
if (is.null(args[["assay_name"]])) {print("No assay name provided, usign 'logcounts' as default.")
  args[["assay_name"]] <- 'logcounts'}
if (is.null(args[["output"]])) {stop("Provide valid path to output h5ad object.")}

#read object
sce <- readRDS(args[["input"]])
#convert sce2h5ad
sceasy:::sce2anndata(obj = sce, outFile = args[["output"]], main_layer = args[["assay_name"]], transfer_layers = NULL)

print("SCE successfully converted to H5ad object.")

