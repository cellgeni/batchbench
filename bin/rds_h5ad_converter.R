#!/usr/bin/env Rscript

library(reticulate)  
library(SummarizedExperiment) 
library(SingleCellExperiment) 

library(loomR)
library(sceasy)


args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {
    print("Provide a valid input file name --> rds file")
}
if (is.null(args[["output"]])) {
    print("Provide a valid output file name --> h5ad file")
}

dataset <- readRDS(args[["input"]])

if (class(dataset) == "SingleCellExperiment"){
  #Those methods that return an embedding, don't have corrected as an assay
  if("corrected_embedding" %in% reducedDimNames(dataset)){
  sceasy:::sce2anndata(obj = dataset, outFile = args[["output"]], main_layer = "logcounts", transfer_layers = NULL)
  }
  #Methods that correct the expression matrix, do have corrected as an assay
  else {
  sceasy:::sce2anndata(obj = dataset, outFile = args[["output"]], main_layer = "corrected", transfer_layers = NULL)
  }
}

if (class(dataset) == "Seurat"){
  sceasy:::seurat2anndata(obj = dataset, outFile = args[["output"]], assay = 'integrated', main_layer = 'data', transfer_layers = NULL )
}
