#!/usr/bin/env Rscript

#libraries
#library(SingleCellExperiment) #object processing
#library(scater) #object processing
library(scran) #mnnCorrect

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {
  print("Provide a valid input file name --> RDS file")
    }
if (is.null(args[["output"]])) {
  print("Provide a valid outpsut file name --> RDS file")
          }
#input file
dataset <- readRDS(args[["input"]])

#create a list of the batches in the dataset
#len <- length(names(table(dataset$Batch)))
batch_vector <- as.character(dataset$Batch)
len <- length(table(as.character(batch_vector)))

batch_list <- lapply(1:len, function(x) {logcounts(dataset[, batch_vector == names(table(batch_vector))[x]])})

t1 = Sys.time()
#run mnnCorrect
corrected <- do.call('mnnCorrect', c(batch_list, c(k=30, sigma=0.1, cos.norm.in=TRUE, svd.dim=2)))
  
  #fix the names of corrected matrices                                  
  for (i in 2:len) {
    colnames(corrected$corrected[[i]]) <- colnames(dataset[ ,batch_vector == names(table(batch_vector)[i])])
    rownames(corrected$corrected[[i]]) <- rownames(dataset)
}
#append the corrected expression matrix to object  
assay(dataset, "corrected") <- do.call(cbind, corrected$corrected)
  
t2 = Sys.time()
print(t2 - t1)
 
saveRDS(dataset, file = args[["output"]])
print("congratulations, this worked!!!")
