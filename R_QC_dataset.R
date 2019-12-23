#!/usr/bin/env Rscript

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {print("Provide a valid input file name --> RDS file")}
if (is.null(args[["bt_thres"]])) {
  print("Provide a valid bt_thres value --> float. Minimum proportion of cells required for a BATCH to be present")
  }
if (is.null(args[["ct_thres"]])) {
  print("Provide a valid ct_thres value --> float. Minimum proportion of cells required for a CELL TYPE to be present")
}
if (is.null(args[["min_genes"]])) {
  print("Provide a valid min_genes value --> integer. Minimum number of genes to be expressed in a cell")
}
if (is.null(args[["min_cells"]])) {
  print("Provide a valid min_cells value --> integer. Minimum number of cells for a gene to be expressed in")
}
if (is.null(args[["output"]])) {print("Provide a valid outpsut file name --> RDS file")}

#inputs
dataset <- readRDS(args[["input"]])
bt_thres <- as.numeric(args[["bt_thres"]])
ct_thres <- as.numeric(args[["ct_thres"]])
min_genes <- as.integer(args[["min_genes"]])
min_cells <- as.integer(args[["min_cells"]])
#parameters
batch_vector <- as.character(dataset$Batch)
cell_type_vector <- as.character(dataset$cell_type1)

#print table with batches/cell types filtering information
print_filtering<- function(dataset, filter_vec, threshold, meta_name){
  
  cell_count <- table(filter_vec)[order(table(filter_vec), decreasing = T)]
  
  print(paste("**", meta_name ,  "containing less than:",  threshold,  "of total cells are removed."))
  print(paste("** TABLE:", meta_name, "filtered based on our threshold:")) 
  #dataframe informing about the filtering about to be done
  exclude_df <- data.frame("X_X" = names(cell_count), 
             "n_cells" = as.vector(cell_count), 
             "percent_cells" = as.vector(cell_count)/ncol(dataset), 
             "Excluded"= as.vector(cell_count)/ncol(dataset) < threshold
             )
  names(exclude_df) <- c(meta_name, "n_cells", "percent_cells", "Excluded")
  print(exclude_df)
  
  removal_names <- exclude_df[meta_name][exclude_df["Excluded"] == T]
  return(removal_names)
}
pre_processing <- function(dataset, min_genes, min_cells){
  #QC pre processing funciton to remove: cells annotated as NA in batch or cell type, 
	#genes expressed in less than min_cells,and cells with less than min_genes expressed.
	n_NAs_b <- length(is.na(as.character(dataset$Batch))[is.na(as.character(dataset$Batch)) == T])
 	dataset <- dataset[,!is.na(as.character(dataset$Batch))]
	print(paste0("Cells annotated as NA for dataset$Batch = ", n_NAs_b, ". Removed"))
	
	n_NAs_ct <- length(is.na(as.character(dataset$cell_type1))[is.na(as.character(dataset$cell_type1)) == T])
 	dataset <- dataset[,!is.na(as.character(dataset$cell_type1))]
	print(paste0("Cells annotated as NA for dataset$cell_type1 = ", n_NAs_ct, ". Removed"))
	
  	dataset <- dataset[apply(logcounts(dataset), 1, function(x) sum(x > 0) >= min_cells),] 
        dataset <- dataset[, apply(logcounts(dataset), 2, function(x) sum(x > 0) >= min_genes)]
 	#Note: we split the apply by axis to find equivalency with the Py_QC.py script  
  
  	print(paste("** CELLS with less than", min_genes, "genes expressed filtered out."))
  	print(paste("** GENES expressed in less than", min_cells, "cells filtered out."))
  
  	print("** Dimensions after preprocessing:")
  	print(dim(dataset))
  	return(dataset)
}
#remove those batch / cell types present less than X_threshold times
remove_elements <- function (dataset, filter_vec, threshold, meta_name){
  removal_names <- print_filtering(dataset, filter_vec, threshold, meta_name)
  
    for(i in 1:length(removal_names)){
      if(length(removal_names) == 0){next()}
      else{
        dataset <- dataset[, dataset$Batch != removal_names[[i]]]
        dataset <- dataset[, dataset$cell_type1 != removal_names[[i]]]
      }
    }
  
  return(dataset)
}

print("** Dimensions pre-filtering:")
print(dim(dataset))
#apply functions!
dataset <- pre_processing(dataset, min_genes = min_genes, min_cells = min_cells)
dataset <- remove_elements(dataset, filter_vec = as.character(dataset$Batch), threshold = bt_thres, meta_name = "Batch/es")
dataset <- remove_elements(dataset, filter_vec = as.character(dataset$cell_type1), threshold = ct_thres, meta_name = "Cell type/s")
print("** Dimensions post-filtering:")
print(dim(dataset))

#save output
saveRDS(dataset, args[["output"]])
