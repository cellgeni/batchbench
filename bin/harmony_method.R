#!/usr/bin/env Rscript
  
suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(harmony))

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {
  print("Provide a valid input file name --> RDS file")
  }
if (is.null(args[["assay_name"]])) {
  print("Assay to use")
  }
if (is.null(args[["n_pcs"]])) {
  print("Number of PCs for embedding")
  }
if (is.null(args[["output"]])) {
  print("Provide a valid outpsut file name --> RDS file")
    }

# input
dataset <- readRDS(args[["input"]])
# run PCA
dataset <- runPCA(dataset, exprs_values = args[["assay_name"]], ncomponents = args[["n_pcs"]])
pca <- dataset@reducedDims@listData[["PCA"]]
# cell batch label vector
batch_vector <-  as.character(dataset$Batch)
# run Harmony
dataset@reducedDims@listData[['corrected_embedding']] <- HarmonyMatrix(pca, batch_vector, theta=4, do_pca = F)
# Add rownames to corrected embedding 
rownames(dataset@reducedDims@listData[['corrected_embedding']]) <- colnames(dataset)
# save object with corrected embedding
saveRDS(dataset, args[["output"]])
print("congratulations, HARMONY worked!")
