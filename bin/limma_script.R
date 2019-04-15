#!/usr/bin/env Rscript
  
#libraries
library(SingleCellExperiment) #object processing
library(limma)


args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {
  print("Provide a valid input file name --> RDS file")
}
if (is.null(args[["output"]])) {
  print("Provide a valid outpsut file name --> RDS file")
}

#input file
dataset <- readRDS(args[["input"]])

#cell batch label vector
batch_vector <- as.character(dataset$Batch)
#run limma
t1 = Sys.time()

assay(dataset, "corrected") <- removeBatchEffect(x = assay(dataset, "logcounts"), batch = batch_vector)

t2 = Sys.time()
print(t2-t1)

saveRDS(dataset, file = args[["output"]])

print("congratulations, LIMMA worked!!!")
