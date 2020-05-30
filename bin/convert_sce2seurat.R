#!/usr/bin/env Rscript

#Convert SCE objects to Seurat.

#TODO
##Include assay name specification option
##One may decide to include the integrated assay into @data and the uncorrected assay into @counts
##MAYBE: Specify which low-D embedding to incorporate in the Seurat object. --> For Clustering analysis!
##Incorporate argparser

suppressPackageStartupMessages(require(Seurat))

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {stop("Provide valid path to input SCE object.")}
#if (is.null(args[["assay_name"]])) {stop("")}
if (is.null(args[["output"]])) {stop("Provide valid path to output Seurat object.")}

sce <- readRDS(args[["input"]])
#check "data_assay" is present in assay names of the SCE object
#if(!(args[["assay_name"]] %in% assayNames(sce))) stop("The provided assay name is not present in the input SCE object")
#convert sce to seurat
sce2seurat <- as.Seurat(sce, assay = "corrected", data = "corrected", counts = "logcounts")
#save object
saveRDS(sce2seurat, file = args[["output"]])
print("SCE successfully converted to Seurat object!")


