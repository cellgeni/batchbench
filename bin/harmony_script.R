
#!/usr/bin/env Rscript
  
#libraries
library(SingleCellExperiment) #object processing
library(scater) #object processing
library(harmony) #Harmony
library(magrittr) #Harmony

args <- commandArgs(trailingOnly = TRUE)

#INPUT!
dataset <- readRDS(args[1])

#run PCA
dataset <- runPCA(dataset, method = "prcomp", exprs_values = "logcounts", ncomponents = 10)
pca <- dataset@reducedDims@listData[["PCA"]]
#cell batch label vector
batch_vector <-  dataset$dataset

#run Harmony
dataset@reducedDims@listData[['corrected_embedding']] <- HarmonyMatrix(pca, batch_vector, theta=4)

print("congratulations, HARMONY worked!!!")

#OUTPUT!
saveRDS(dataset, args[2])