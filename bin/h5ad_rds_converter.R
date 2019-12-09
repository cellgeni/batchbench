#!/usr/bin/env Rscript

library(reticulate)
library(SingleCellExperiment)
library(SummarizedExperiment)

anndata <- reticulate::import('anndata', convert = F)
pandas <- reticulate::import('pandas', convert = F)
reticulate::import('scipy', convert = F)

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {
  print(".h5ad file")
}
if (is.null(args[["output"]])) {
  print(".rds file")
}


anndata2sce <- function(obj)
{
  #read h5ad file
  dataset <- anndata$read_h5ad(obj)
  #extract cell names
  cell_names <- reticulate::py_to_r(pandas$Series(dataset$obs_names)$values)
  #extract gene names
  gene_names <- reticulate::py_to_r(pandas$Series(dataset$var_names)$values)
  #counts matrix
  counts_matrix <- t(as.matrix(reticulate::py_to_r(dataset$X)$todense()))
  rownames(counts_matrix) <- as.character(gene_names)
  colnames(counts_matrix) <- as.character(cell_names)
  #cell metadata
  cell_metadata <- reticulate::py_to_r(dataset$obs)
  #embeddings
  emb_list <- list()
  emb_names = reticulate::py_to_r(dataset$obsm_keys())
  if (length(emb_names) > 0) {
    for (i in 1:length(emb_names)) {
      emb_list[[i]] <- reticulate::py_to_r(dataset$obsm[emb_names[[i]]])
    }
    names(emb_list) <- emb_names
  }
  #build sce
  sce <- SingleCellExperiment(assays = list(logcounts = counts_matrix),
                              colData = data.frame(cell_metadata),
                              rowData = list('feature_names' =  as.character(gene_names)),
                              reducedDims = emb_list
                              )
  return(sce)
}


#exec function
sce <- anndata2sce(obj = args[['input']])
#save file
saveRDS(sce, args[["output"]])

