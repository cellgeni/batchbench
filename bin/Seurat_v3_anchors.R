#!/usr/bin/env Rscript
library(Seurat)

#input
args <-R.utils::commandArgs(asValues = TRUE)

if (is.null(args[["input"]])) {
  print("Provide a valid input file name --> RDS file")
  }
if (is.null(args[["output"]])) {
  print("Provide a valid outpsut file name --> RDS file")
    }
dataset <- readRDS(args[["input"]])
print("A")
dataset <- as.Seurat(dataset,counts = "logcounts", data = "logcounts")
print("B")
batch_vector <- as.character(dataset$Batch)
len <- length(names(table(batch_vector)))

#split object into batches
batch_list <- lapply(1:len, function(x) {abc <- dataset[, dataset$Batch == names(table(batch_vector))[x]]})
print("C")
#Normalize and find HVG
for (i in 1:length(batch_list)) {
  batch_list[[i]] <- NormalizeData(object = batch_list[[i]], verbose = FALSE)
  batch_list[[i]] <- FindVariableFeatures(object = batch_list[[i]], 
                                          selection.method = "dispersion", nfeatures = 2000, verbose = FALSE)
}
print("D")

if(any(sapply(batch_list, ncol)) < 200) {
  k_filter <- (min(sapply(batch_list, ncol)))
  }else{k_filter =200}
#Find integration anchors
anchors <- FindIntegrationAnchors(object.list = batch_list, dims = 1:30, k.filter = k_filter)
#integrate subsets
integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
print("E")
#Output
saveRDS(integrated, args[["output"]])
print("F")
