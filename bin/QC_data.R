#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to rds input file'
  ),
  make_option(
    c("-a", "--assay_name"),
    action = "store",
    default = "logcounts",
    type = 'character',
    help = 'Assay name'
  ),
  make_option(
    c("-b", "--batch_key"),
    action = "store",
    default = "Batch",
    type = 'character',
    help = 'Batch key in cell metadata'
  ),
  make_option(
    c("-c", "--celltype_key"),
    action = "store",
    default = "cell_type1",
    type = 'character',
    help = 'Cell type key in cell metadata'
  ),
  make_option(
    c("-m", "--min_cells"),
    action = "store",
    default = 3,
    type = 'integer',
    help = 'Minimum number of cells for a gene to be expressed in.'
  ),
  make_option(
    c("-g", "--min_genes"),
    action = "store",
    default = 200,
    type = 'integer',
    help = 'Minimum number of genes expressed per cell for a cell to be considered.'
  ),
  make_option(
    c("-t", "--bt_thres"),
    action = "store",
    default = 0.05,
    type = 'numeric',
    help = 'Minimum proportion of total cells for a batch to be considered'
  ),
  make_option(
    c("-e", "--ct_thres"),
    action = "store",
    default = 0.01,
    type = 'numeric',
    help = 'Minimum proportion of total cells for a cell type to be considered'
  ),
  make_option(
    c("-o", "--output_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to rds output file'
  )
)
opt <- parse_args(OptionParser(option_list=option_list))

# args
assay_name <- opt$assay_name
batch_key <- opt$batch_key
celltype_key <- opt$celltype_key
min_genes <- opt$min_genes
min_cells <- opt$min_cells
bt_thres <- opt$bt_thres
ct_thres <- opt$ct_thres

# Read input object
dataset <- readRDS(opt$input_object)
# Vectors with batch and cell type annotation
batch_vector <- as.character(dataset[[batch_key]])
cell_type_vector <- as.character(dataset[[celltype_key]])

# print table with batch and cell types filtering information
print_filtering <- function(dataset, filter_vec, threshold, meta_name){
  
  cell_count <- table(filter_vec)[order(table(filter_vec), decreasing = T)]
  print(paste("**", meta_name ,  "containing less than:",  threshold,  "of total cells are removed."))
  print(paste("** TABLE:", meta_name, "filtered based on our threshold:"))
  # dataframe informing about the filtering about to be done
  filtering_df <- data.frame("X_X" = names(cell_count),
                             "n_cells" = as.vector(cell_count),
                             "percent_cells" = as.vector(cell_count)/ncol(dataset),
                             "Excluded"= as.vector(cell_count)/ncol(dataset) < threshold
  )
  colnames(filtering_df)[1] <- meta_name
  print(filtering_df)
  
  removal_names <- filtering_df[meta_name][filtering_df["Excluded"] == T]
  return(removal_names)
}

# QC pre processing funciton to remove: cells annotated as NA in batch or cell type, genes expressed in less than min_cells,and cells with less than min_genes expressed.
pre_processing <- function(dataset, min_genes, min_cells){
  n_NAs_b <- length(is.na(as.character(dataset[[batch_key]]))[is.na(as.character(dataset[[batch_key]])) == T])
  dataset <- dataset[,!is.na(as.character(dataset[[batch_key]]))]
  print(paste0("Cells annotated as NA for dataset$Bacth = ", n_NAs_b, ". Removed"))
  
  n_NAs_ct <- length(is.na(as.character(dataset[[celltype_key]]))[is.na(as.character(dataset[[celltype_key]])) == T])
  dataset <- dataset[,!is.na(as.character(dataset[[celltype_key]]))]
  print(paste0("Cells annotated as NA for dataset$cell_type1 = ", n_NAs_ct, ". Removed"))
  
  dataset <- dataset[apply(logcounts(dataset), 1, function(x) sum(x > 0) >= min_cells),]
  dataset <- dataset[, apply(logcounts(dataset), 2, function(x) sum(x > 0) >= min_genes)]
  
  print(paste("** CELLS with less than", min_genes, "genes expressed filtered out."))
  print(paste("** GENES expressed in less than", min_cells, "cells filtered out."))
  
  print("** Dimensions after preprocessing:")
  print(dim(dataset))
  return(dataset)
}
# Remove those batch / cell types present less than X_threshold times
remove_elements <- function (dataset, filter_vec, threshold, meta_name){
  removal_names <- print_filtering(dataset, filter_vec, threshold, meta_name)
  
  for(i in 1:length(removal_names)){
    if(length(removal_names) == 0){next()}
    else{
      dataset <- dataset[, dataset[[batch_key]] != removal_names[[i]]]
      dataset <- dataset[, dataset[[celltype_key]] != removal_names[[i]]]
    }
  }
  
  return(dataset)
}
# Remove items with zero variance
rm_zero_var <- function(dataset, assay, axis){
  if(!(axis %in% c(1, 2))){stop("Axis provided must be 1 or 2!")}
  ax_zero_var <- as.numeric(which(apply(assay(dataset, assay), axis, var) == 0))
  print(paste0("Removing ", length(ax_zero_var), " items with zero variance in axis = ", axis))
  if(length(ax_zero_var) > 0){
    if(axis == 1) dataset <- dataset[-ax_zero_var, ]
    if(axis == 2) dataset <- dataset[, -ax_zero_var]
  } 
  dataset
}

print("** Dimensions pre-filtering:")
print(dim(dataset))

# apply filtering functions!
dataset <- pre_processing(dataset, min_genes = min_genes, min_cells = min_cells)
# Remove batches
dataset <- remove_elements(dataset, filter_vec = dataset[[batch_key]], threshold = bt_thres, meta_name = "Batch/es")
# Remove cell types
dataset <- remove_elements(dataset, filter_vec = dataset[[celltype_key]], threshold = ct_thres, meta_name = "Cell type/s")
# Remove cells with zero variance
dataset <- rm_zero_var(dataset, assay = assay_name, axis = 1)
# Remove genes with zero variance
dataset <- rm_zero_var(dataset, assay = assay_name, axis = 2)

print("** Dimensions post-filtering:")
print(dim(dataset))

# save output
saveRDS(dataset, opt$output_object)
print("Dataset successfully QCded")
