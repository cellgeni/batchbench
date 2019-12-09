#!/usr/bin/env Rscript

#library(scran) #fastMNN is depricated in scran, now on Batchelor package
library(batchelor)

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
batch_vector <- as.character(dataset$Batch)
len <- length(table(as.character(batch_vector)))

batch_list <- lapply(1:len, function(x) {logcounts(dataset[, batch_vector == names(table(batch_vector))[x]])})

t1 = Sys.time()

#run fastMNN
correction <- do.call('fastMNN', c(batch_list, c(k=30, d = 25,  cos.norm=T,  pc.input = F)))

#attach the batch corrected low_d embedding to reducedDims of the SCE object
dataset@reducedDims@listData[['corrected_embedding']] <- correction@reducedDims@listData[["corrected"]]

t2 = Sys.time()
print(t2 - t1)
 
saveRDS(dataset, file = args[["output"]])
print("congratulations, this worked!!!")

