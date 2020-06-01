#!/usr/bin/env Rscript
suppressPackageStartupMessages(require(Seurat))

#input
args <-R.utils::commandArgs(asValues = TRUE)

if (is.null(args[["input"]])) {
  print("Provide a valid input file name --> RDS file")
  }
if (is.null(args[["assay_name"]])) {
  print("Assay to use")
  }
if (is.null(args[["batch_key"]])) {
  print("Cell column where batch is defined")
  }
if (is.null(args[["hvg_method"]])) {
  print("Select method to FindVariableFeatures: one of 'vst', 'mean.var.plot' or 'dispersion'")
  }
if (is.null(args[["n_features"]])) {
  print("N of variable features to find in FindVariableFeatures")
  }
if (is.null(args[["n_anchors"]])) {
  print("Number of anchors to find to perform integration")
  }
if (is.null(args[["verbose"]])) {
  print("V")
  }
if (is.null(args[["output"]])) {
  print("Provide a valid outpsut file name --> RDS file")
    }
# input dataset
dataset <- readRDS(args[["input"]])
# convert to seurat object
dataset <- as.Seurat(dataset, counts = args[["assay_name"]], data = args[["assay_name"]])
batch_vector <- as.character(dataset[[args[["batch_key"]]]])
batch_names  <- names(table(batch_vector))
N_batches <- length(batch_names)

# split seurat object into batches
batch_list <- lapply(1:N_batches, function(x) {abc <- dataset[, dataset[[args[["batch_key"]]]] == batch_names[x]]})
# Normalize and find HVG
for (i in 1:N_batches) {
  batch_list[[i]] <- NormalizeData(object = batch_list[[i]], 
						verbose = FALSE)
  batch_list[[i]] <- FindVariableFeatures(object = batch_list[[i]], 
                                          	selection.method = args[["hvg_method"]], 
						nfeatures = args[["n_features"]], 
						verbose = args[["verbose"]])
}

if(any(sapply(batch_list, ncol)) < 200) {
  k_filter <- (min(sapply(batch_list, ncol)))
  }else{k_filter =200}
# Find integration anchors
anchors <- FindIntegrationAnchors(object.list = batch_list, dims = 1:args[["n_anchors"]], k.filter = k_filter)
# Integrate subsets
integrated <- IntegrateData(anchorset = anchors, dims = 1:args[["n_anchors"]])
# save seurat3 corrected object
saveRDS(integrated, args[["output"]])
