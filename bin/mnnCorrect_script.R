#!/usr/bin/env Rscript

#libraries
#library(SingleCellExperiment) #object processing
#library(scater) #object processing
library(scran) #mnnCorrect

args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}  

#input file
dataset <- readRDS(args[1])


#create a list of the batches in the dataset
len <- length(names(table(dataset$dataset)))

batch_list <- lapply(1:len, function(x) {logcounts(dataset[, dataset$dataset == names(table(dataset$dataset)[x])])})

t1 = Sys.time()
#run mnnCorrect
corrected <- do.call('mnnCorrect', c(batch_list, c(k=30, sigma=0.1, cos.norm.in=TRUE, svd.dim=2)))
  
  #fix the names of corrected matrices                                  
  for (i in 2:len) {
    colnames(corrected$corrected[[i]]) <- colnames(dataset[ ,dataset$dataset == names(table(dataset$dataset)[i])])
    rownames(corrected$corrected[[i]]) <- rownames(dataset)
}
#append the corrected expression matrix to object  
assay(dataset, "corrected") <- do.call(cbind, corrected$corrected)
  
t2 = Sys.time()
print(t2 - t1)
 
saveRDS(dataset, file = args[2])
print("congratulations, this worked!!!")
