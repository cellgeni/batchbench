#!/usr/bin/env Rscript
  
#libraries
library(scater) #object processing
library(harmony) #Harmony
library(magrittr) #Harmony

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {
  print("Provide a valid input file name --> RDS file")
  }
if (is.null(args[["output"]])) {
  print("Provide a valid outpsut file name --> RDS file")
    }

#INPUT!
dataset <- readRDS(args[["input"]])

#run PCA
dataset <- runPCA(dataset, method = "prcomp", exprs_values = "logcounts", ncomponents = 50)
pca <- dataset@reducedDims@listData[["PCA"]]
#cell batch label vector
batch_vector <-  as.character(dataset$Batch)

#run Harmony
dataset@reducedDims@listData[['corrected_embedding']] <- HarmonyMatrix(pca, batch_vector, theta=4)

print("congratulations, HARMONY worked!!!")

#OUTPUT!
saveRDS(dataset, args[["output"]])
