#!/usr/bin/env Rscript

#libraries
library(SingleCellExperiment) #object processing
#library(scater)
library(Seurat)


#input
args <- commandArgs(trailingOnly = TRUE)

dataset <- readRDS(args[1])

#convert SCE object to seurat object
#dataset <- Convert(dataset, to = "seurat")
#if more than one data slot available(get the normcounts)
dataset <- Convert(dataset, to = "seurat", raw.data.slot = "logcounts", data.slot = "logcounts")

#subset batches of dataset object
len <- length(names(table(dataset@meta.data$dataset)))

batch_list <- lapply(1:len, function(x) {SubsetData(dataset, cells.use = (dataset@meta.data[["dataset"]] == names(table(dataset@meta.data$dataset)[x])))})

#fing variable genes and scale data
for (i in 1:len){
  batch_list[[i]] <- FindVariableGenes(batch_list[[i]])
  batch_list[[i]] <- ScaleData(batch_list[[i]])
}

#select genes to use in CC alignment
genes.use <- c()
for (i in 1:len) {
  genes.use <- c(genes.use, rownames(batch_list[[i]]@hvg.info)#[1:5000] #select as many HVG as wanted!
                 )
}
#select those HVG present in more than X subsets
#genes.use <- names(which(table(genes.use) > 2))
for (i in 1:len) {
  genes.use <- genes.use[genes.use %in% rownames(batch_list[[i]]@scale.data)]
}

#runMultiCCA
dataset <- RunMultiCCA(batch_list, genes.use = genes.use, num.ccs = 5)

#OPTIONAL: plot to determine the optimal number of CCs to align the subspaces
#MetageneBicorPlot(dataset_cca, grouping.var = "dataset", dims.eval = 1:10)

#Calculate the ratio of variance explained by PCA to CCA
dataset <- CalcVarExpRatio(object = dataset, reduction.type = "pca",
                              grouping.var = "dataset", dims.use = 1:5)

#OPTIONAL: Subset those cells with cutoff calue = accept.low
#dataset <- SubsetData(dataset_cca, subset.name = "var.ratio.pca",
#                         accept.low = 0.5)
#Align subspaces
dataset<- AlignSubspace(dataset,
                            reduction.type ="cca",
                            grouping.var = "dataset",
                            dims.align = 1:5)


print("Congratulations, multiCCA worked!")

#OUTPUT!
saveRDS(dataset, args[2])
