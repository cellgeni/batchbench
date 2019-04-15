#!/usr/bin/env Rscript
library(Seurat)

#input
args <- commandArgs(trailingOnly = TRUE)
dataset <- readRDS(args[1])

dataset <- as.Seurat(dataset, data = "logcounts")

batch_vector <- as.character(dataset$Batch)
len <- length(names(table(batch_vector)))

#split object into batches
batch_list <- lapply(1:len, function(x) {abc <- dataset[, dataset$Batch == names(table(batch_vector))[x]]})

#Normalize and find HVG
for (i in 1:length(batch_list)) {
  batch_list[[i]] <- NormalizeData(object = batch_list[[i]], verbose = FALSE)
  batch_list[[i]] <- FindVariableFeatures(object = batch_list[[i]], 
                                          selection.method = "dispersion", nfeatures = 2000, verbose = FALSE)
}

#Find integration anchors
anchors <- FindIntegrationAnchors(object.list = batch_list, dims = 1:30)

#integrate subsets
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)

#Output
saveRDS(integrated, args[2])
