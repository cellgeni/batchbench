#!/usr/bin/env Rscript
  
#libraries
library(SingleCellExperiment) #object processing
library(limma)


args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#input file
dataset <- readRDS(args[1])

#cell batch label vector
batch_vector <- dataset$dataset
#run limma
t1 = Sys.time()

assay(dataset, "corrected") <- removeBatchEffect(x = assay(dataset, "logcounts"), batch = batch_vector)

t2 = Sys.time()
print(t2-t1)

saveRDS(dataset, file = args[2])

print("congratulations, LIMMA worked!!!")