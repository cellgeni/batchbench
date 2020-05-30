#!/usr/bin/env Rscript
  
#libraries
library(SingleCellExperiment) #object processing
library(sva) #ComBat


args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {
  print("Provide a valid input file name --> RDS file")
}
if (is.null(args[["output"]])) {
  print("Provide a valid outpsut file name --> RDS file")
}

#input file
dataset <- readRDS(args[["input"]])

#remove those genes with 0 variance, if not ComBat throws error!
dataset <- dataset[rowSums(logcounts(dataset)) > 0, ]

#data
mod_data <- as.data.frame(t(as.matrix(logcounts(dataset))))
# Basic batch removal
mod0 = model.matrix(~ 1, data = mod_data)
#run ComBat
t1 = Sys.time()

assay(dataset, "corrected") <- ComBat(
  dat = t(mod_data),
  batch = as.character(dataset$Batch),
  mod = mod0,
  par.prior = TRUE,
  prior.plots = FALSE
)

t2 = Sys.time()
print(t2-t1)

#Note: an error appears when telling:Error in while (change > conv) { : missing value where TRUE/FALSE needed
#This is fixed, removing those genes with 0 variance across cells!
#This may be also because of -Inf values in logcounts

print("ComBat worked!")

saveRDS(dataset, file = args[["output"]])

print("congratulations, this worked!!!")
