#!/usr/bin/env Rscript

library(scater)

#value arguments
args <-R.utils::commandArgs(asValues=TRUE)
if (is.null(args[["input"]])) {
  print("Provide a valid input file name (Batch corrected object) --> RDS file")
}
if (is.null(args[["output"]])) {
  print("Provide a valid output file name (Per batch and cell type entropy) --> CSV file")
}

#input file
object <- readRDS(args[["input"]])

#runPCA
object <- runPCA(object, exprs_values = "logcounts", ncomponents = 50)

#umap function
umap_coordinates <- function(object){
  object <- runUMAP(object, ncomponents= 2, use_dimred = 'PCA', n_dimred = 30)
  umap <- object@reducedDims@listData[["UMAP"]]
  return(umap)
}
#save results
save_results <- function(x){
  write.table(x, file = args[["output"]], sep = "\t",  row.names = FALSE)
}

#(compute UMAP)
umap <- umap_coordinates(object)
#write result to table
save_results(umap)
