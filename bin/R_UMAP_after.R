#!/usr/bin/env Rscript

#This script is to calculate the UMAP coordinates of the batch-corrected objects.

#libraries
library(scater)

args <-R.utils::commandArgs(asValues=TRUE)
if (is.null(args[["input"]])) {
  print("Provide a valid input file name (Batch corrected object) --> RDS file")
}
if (is.null(args[["output"]])) {
  print("Provide a valid output file name (UMAP coordinates of batch-corrected object) --> CSV file")
}

#umap function
umap_coordinates <- function(object, dr_type,  output_umap){
  object <- runUMAP(object, ncomponents = 2, use_dimred = dr_type, n_neighbors = 30)
  umap <- object@reducedDims@listData[["UMAP"]]
  write.csv(umap, file= output_umap, sep = "",  row.names= colnames(object), col.names = paste0("UMAP_", c(1:ncol(umap))))
}

#umap_seurat_coordinates <- function(object, dr_type,  output_umap){
#  object <- RunUMAP(object = object, reduction.use = dr_type, reduction.name = "UMAP", max.dim = 2)
#  umap <- GetCellEmbeddings(object, reduction.type = "UMAP")
#  write.table(umap, file= output_umap,sep = "\t", row.names= colnames(object@data), col.names = colnames(umap))
#}

#input file
object <- readRDS(args[["input"]])

# 1) SCE objects 
if (class(object) == "SingleCellExperiment"){
  print("The input object is a Single Cell Experiment class object")

# 1.1) tools that correct low_D embeddings
  if ("corrected_embedding" %in% reducedDimNames(object)){
    print("A")
    umap_coordinates(object, dr_type = "corrected_embedding", output_umap = args[["output"]])
    print("UMAP coordinates computed succesfully")
  }

# 1.2) methods that correct expression matrix
  if ("corrected" %in% names(assays(object))){
    object <- runPCA(object, exprs_values = "logcounts", ncomponents = 30)
    umap_coordinates(object, dr_type = "PCA", output_umap = args[["output"]])
    print("UMAP coordinates computed succesfully")
  }
}

#2) Seurat objects
# 2.1) Seurat v2 multiCCA object
if (class(object) == "seurat"){
    library(Seurat)
    print("The input object is a Seurat object")
    print("Object corrected by Seurat_v2_multiCCA")
    #convert Seurat v2 to SCE object
    object <- SingleCellExperiment(assays = list(logcounts = as.matrix(object@data)),
                     colData = as.data.frame(object@meta.data) ,
                     reducedDims =  list( "cca" =  attr(object@dr[["cca"]], "cell.embeddings"),
                                          "cca.aligned" = attr(object@dr[["cca.aligned"]], "cell.embeddings")))
    umap_coordinates(object, dr_type = "cca.aligned", output_umap = args[["output"]])
    print("UMAP coordinates computed succesfully")
    
}


# 2.2) Seurat v3 Anchors object
if (class(object) == "Seurat"){

    library(Seurat)
    print("The input object is a Seurat object")
    print("Object corrected by Seurat_v3_anchors")    
    
    object <- SingleCellExperiment(assays = list("corrected" = as.matrix(object@assays[["integrated"]]@data)),
                                   colData = as.data.frame(object@meta.data))
    #convert Seurat v3 to SCE object                             
    object <- runPCA(object, exprs_values = "corrected", ncomponents = 30)
    umap_coordinates(object, dr_type = "PCA", output_umap = args[["output"]])
    print("UMAP coordinates computed succesfully")
    
}
 

