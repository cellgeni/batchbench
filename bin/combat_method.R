#!/usr/bin/env Rscript
  
suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(sva))

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {
  print("Provide a valid input file name --> RDS file")
}
if (is.null(args[["output"]])) {
  print("Provide a valid outpsut file name --> RDS file")
}

# input file
dataset <- readRDS(args[["input"]])
# remove those genes with 0 variance, if not ComBat throws error!
dataset <- dataset[rowSums(logcounts(dataset)) > 0, ] # in principle already removed in QC
# run ComBat
mod_data <- as.data.frame(t(as.matrix(logcounts(dataset))))
# Basic batch removal
mod0 = model.matrix(~ 1, data = mod_data)

assay(dataset, "corrected") <- ComBat(
  dat = t(mod_data),
  batch = as.character(dataset$Batch),
  mod = mod0,
  par.prior = TRUE,
  prior.plots = FALSE
)
# save corrected object
saveRDS(dataset, file = args[["output"]])
print("ComBat worked!")
