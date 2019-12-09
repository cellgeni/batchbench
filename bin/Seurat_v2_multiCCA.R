#!/usr/bin/env Rscript
library(Seurat)
library(SingleCellExperiment)

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {
  print("Provide a valid input file name --> RDS file")
}
if (is.null(args[["output"]])) {
  print("Provide a valid outpsut file name --> RDS file")
}

#INPUT!
dataset <- readRDS(args[["input"]])

#if more than one data slot available(get the normcounts)
dataset <- Convert(dataset, to = "seurat", raw.data.slot = "logcounts", data.slot = "logcounts")

#subset batches of dataset object
batch_vector <- as.character(dataset@meta.data$Batch)
len <- length(names(table(batch_vector)))

batch_list <- lapply(1:len, function(x) {SubsetData(dataset, cells.use = (dataset@meta.data[["Batch"]] == names(table(batch_vector)[x])))})

#find variable genes and scale data
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

#if there are 2 batches
if (len < 3){ 
  #runCCA
  dataset <- RunCCA(object = batch_list[[1]], object2 = batch_list[[2]], genes.use = genes.use, num.cc = 15)
  
  #Calculate the ratio of variance explained by PCA to CCA
  dataset <- CalcVarExpRatio(object = dataset, reduction.type = "pca",
                             grouping.var = "Batch", dims.use = 1:15)
  
  #OPTIONAL: Subset those cells with cutoff calue = accept.low
  #dataset <- SubsetData(dataset_cca, subset.name = "var.ratio.pca",
  #                         accept.low = 0.5)
  #Align subspaces
  dataset<- AlignSubspace(dataset,
                          reduction.type ="cca",
                          grouping.var = "Batch",
                          dims.align = 1:15)
  print("Congratulations, CCA worked!")
  
  saveRDS(dataset, args[["output"]])
  print("Object saved")
}

#if more that 2 batches
if (len >= 3){ 
  #runMultiCCA
  dataset <- RunMultiCCA(batch_list, genes.use = genes.use, num.ccs = 15)
  
  #Calculate the ratio of variance explained by PCA to CCA
  dataset <- CalcVarExpRatio(object = dataset, reduction.type = "pca",
                             grouping.var = "Batch", dims.use = 1:15)
  
  #OPTIONAL: Subset those cells with cutoff calue = accept.low
  #dataset <- SubsetData(dataset_cca, subset.name = "var.ratio.pca",
  #                         accept.low = 0.5)
  #Align subspaces
  dataset<- AlignSubspace(dataset,
                          reduction.type ="cca",
                          grouping.var = "Batch",
                          dims.align = 1:15)
  print("Congratulations, multiCCA worked!")
  saveRDS(dataset, args[["output"]])
  print("Object saved")
}
