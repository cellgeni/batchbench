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
    help = 'Counts assay to add to the h5ad object'
  ),
  make_option(
    c("-c", "--corrected_assay"),
    action = "store",
    default = "corrected",
    type = 'character',
    help = 'Corrected counts assay name'
  ),
  make_option(
    c("-b", "--batch_key"),
    action = "store",
    default = "Batch",
    type = 'character',
    help = 'Batch key in cell metadata'
  ),
  make_option(
    c("-t", "--celltype_key"),
    action = "store",
    default = "cell_type1",
    type = 'character',
    help = 'Cell type key in cell metadata'
  ),
  make_option(
    c("-l", "--logFC.thres"),
    action = "store",
    default = 0.5,
    type = 'numeric',
    help = 'logFC threshold for marker gene detection.'
  ),
  make_option(
    c("-o", "--output_markers"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to csv markers file'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(Seurat))

# args
assay_name <- opt$assay_name
batch_key <- opt$batch_key
celltype_key <- opt$celltype_key
corrected_assay <- opt$corrected_assay
logFC.thres <- opt$logFC.thres

# Find markers function 
find_markers_seurat <- function(dataset, batch){
  Idents(dataset) <- dataset[[celltype_key]]
  markers_list <- list() #list of dataframes with marker gene info
  for (i in 1:length(levels(dataset))){
    tryCatch({
      print(paste0("** Cell type ", i, "--", levels(dataset)[i]))
      #if computing over the merged dataset, 
      if(batch == F){
        markers_list[["merged_dataset"]][[levels(dataset)[i]]] <- FindMarkers(object = dataset, slot = "data", ident.1 = levels(dataset)[i], ident.2 = NULL, min.pct = 0.5, logfc.threshold =logFC.thres )
      }
      if (batch == T){
        markers_list[[levels(dataset)[i]]] <- FindMarkers(object = dataset,  slot = "data", ident.1 = levels(dataset)[i], ident.2 = NULL, min.pct = 0.5,  logfc.threshold = logFC.thres)
      }
      # continue if error
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  return(markers_list)
}
#. Find_markers_by_batch
find_markers_by_batch <- function(dataset){
  #table_batches <- table(as.character(dataset[[batch_key]]))
  #batch_list <- lapply(1:length(table_batches), function(x) dataset[, dataset[[batch_key]]==names(table_batches)[x]])
  #names(batch_list) <- names(table_batches)
  batch_list <- SplitObject(dataset, split.by = batch_key)
  by_batch_markers_list <- list()
  #loop through the batches
  for(i in 1:length(batch_list)){
    print(paste0("* Batch number ", i, ": ", names(batch_list)[i]))
    #execute function
    markers_list <- find_markers_seurat(batch_list[[i]], batch = T)
    by_batch_markers_list[[paste0("Batch_", names(batch_list)[i])]] <- markers_list
  }
  return(by_batch_markers_list)
}

#EXECUTE MARKERS
#1. read input object
dataset <- readRDS(opt$input_object)
#2. Find Markers merged dataset
markers_list <- find_markers_seurat(dataset, batch = F)
#3. Find Markers by batch
by_batch_markers_list <- find_markers_by_batch(dataset)
## Combine lists
total_markers_list <- c(markers_list, by_batch_markers_list)
##Save output
saveRDS(total_markers_list, file=opt$output_markers)
print("Marker genes successfully computed for merged dataset and individual batches.")
