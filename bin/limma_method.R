#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(limma))


args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {
  print("Provide a valid input file name --> RDS file")
}
if (is.null(args[["output"]])) {
  print("Provide a valid outpsut file name --> RDS file")
}

# input file
dataset <- readRDS(args[["input"]])
# cell batch label vector
batch_vector <- as.character(dataset$Batch)
# Run limma
assay(dataset, "corrected") <- removeBatchEffect(x = assay(dataset, "logcounts"), batch = batch_vector)
# save corrected object
saveRDS(dataset, file = args[["output"]])
print("congratulations, Limma worked!")
