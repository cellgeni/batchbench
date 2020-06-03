#!/usr/bin/env Rscript

#Convert Seurat to SCE

#TODO
##Sceasy package doe snot support the conversion Seurat to SCE, although they say so in the documentation. 
##Report as GitHub ISSUE!

suppressPackageStartupMessages(require(Seurat))
suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(sceasy))

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {stop("Provide valid path to input Seurat object.")}
if (is.null(args[["output"]])) {stop("Provide valid path to output SCE object.")}

#read Seurat object
seurat_obj <- readRDS(args[["input"]])
#convert Seurat to SCE and Save
if (is.null(args[["assay_name"]])) { print("Seurat assay to extract not specified. Setting assay to 'RNA'")
  assay_name <- "RNA" }

#Convert Seurat to SCE
seurat2sce <- as.SingleCellExperiment(seurat_obj, assay = assay_name)

print("Seurat successfully converted to SCE object!")