#!/usr/bin/env Rscript

#libraries

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args[["input"]])) {
  print("Provide a valid input file name --> RDS file")
  }
if (is.null(args[["batch_threshold"]])) {
  print("Provide a valid batch_threshold value --> integer")
  }
if (is.null(args[["cell_type_threshold"]])) {
  print("Provide a valid cell_type_threshold value --> integer")
  }
if (is.null(args[["output"]])) {
  print("Provide a valid outpsut file name --> RDS file")
    }

#input file
object <- readRDS(args[["input"]])
#parameters
batch_vector <- object$Batch
batch_threshold = as.integer(args[["batch_threshold"]])
cell_type_vector <- object$cell_type1
cell_type_threshold = as.integer(args[["cell_type_threshold"]])

print("Dimensions pre-filtering:")
print(dim(object))

#remove those batch / cell types present less than X_threshold times
remove_elements <- function (object, filter_vec, threshold, info){
  removal_names <- names(table(filter_vec)[table(filter_vec) <= threshold])
  if (length(removal_names) == 0) print(paste("There is no", info, "with less than", threshold, "cells"))
  if (length(removal_names) != 0){
    print(paste("There are", length(removal_names), info,  "with less than", threshold, "cells :", paste(removal_names, collapse = ", "), "removed!"))
    for(i in 1:length(removal_names)){
      object <- object[, object$Batch != removal_names[[i]]]
      object <- object[, object$cell_type1 != removal_names[[i]]]
    }
  }
  return(object)
}

#apply function
object <- remove_elements(object, filter_vec = batch_vector, threshold = batch_threshold, info = "Batch/es")
object <- remove_elements(object, filter_vec = cell_type_vector, threshold = cell_type_threshold, info = "Cell type/s")

print("Dimensions post-filtering:")
print(dim(object))

#save output
saveRDS(object, args[["output"]])