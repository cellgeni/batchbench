#!/usr/bin/env Rscript

#This script is to calculate the UMAP coordinates of the batch-corrected objects.

#libraries
library(scater)

args <-R.utils::commandArgs(asValues=TRUE)
if (is.null(args[["input"]])) {
  print("Provide a valid input file name (Batch corrected object) --> RDS file")
}
if (is.null(args[["output_umap"]])) {
  print("Provide a valid output file name (UMAP coordinates of batch-corrected object) --> txt file")
}

#umap function
umap_coordinates <- function(object, dr_type,  output_umap){
  object <- runUMAP(object, ncomponents = 2, use_dimred = dr_type, n_dimred = 30)
  umap <- object@reducedDims@listData[["UMAP"]]
  write.table(umap, file= output_umap, sep = "\t",  row.names= colnames(object), col.names = paste0("UMAP_", c(1:ncol(umap))))
}
umap_seurat_coordinates <- function(object, dr_type,  output_umap){
  object <- RunUMAP(object = object, reduction.use = dr_type, reduction.name = "UMAP", max.dim = 2)
  umap <- GetCellEmbeddings(object, reduction.type = "UMAP")
  write.table(umap, file= output_umap,sep = "\t", row.names= colnames(object@data), col.names = colnames(umap))
}

#input file
object <- readRDS(args[["input"]])

# 1) SCE objects 
if (class(object) == "SingleCellExperiment"){
  print("The input object is a Single Cell Experiment class object")

# 1.1) tools that correct low_D embeddings
  if ("corrected_embedding" %in% reducedDimNames(object)){
      print("A")
      umap_coordinates(object, dr_type = "corrected_embedding", output_umap = args[["output_umap"]])
  }

# 1.2) methods that correct expression matrix
  if ("corrected" %in% names(assays(object))){
    object <- runPCA(object, method = "prcomp", exprs_values = "logcounts", ncomponents = 30)
    umap_coordinates(object, dr_type = "PCA", output_umap = args[["output_umap"]])
  }
}

#2) Seurat objects
if (class(object) == "seurat"){
  library(Seurat)
  print("The input object is a Seurat object")
  umap_seurat_coordinates(object, dr_type = "cca", output_umap = args[["output_umap"]])
umap_seurat_coordinates(object, dr_type = "cca.aligned", output_umap = args[["output_umap"]])
}