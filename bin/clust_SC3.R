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
    c("-c", "--celltype_key"),
    action = "store",
    default = "cell_type1",
    type = 'character',
    help = 'Cell type key in cell metadata'
  ),
  make_option(
    c("-o", "--output_clusters"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to h5ad output file'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(SC3))

#build sce with corrected counts matrix
exp_sce <- function(dataset){
  print("inside!")
  counts_mat <- dataset@assays[["corrected"]]
  col_data <- dataset@colData
  sce <- SingleCellExperiment(assays = list(logcounts = counts_mat), 
                              colData = col_data)
  return(sce)
}

#common modifications applied to sce object
common_modif_sc3 <- function(sce){
  rowData(sce)$feature_symbol <- rownames(sce)
  counts(sce) <- logcounts(sce)
  return(sce)
}
#run sc3 function
run_sc3 <- function(sce){
  k <- length(table(sce[[celltype_key]])[table(sce[[celltype_key]]) > 0])
  k_vec <- c(k)
  print(paste0("k values on clusteing:", k_vec))
  #biology = T enables calculation of DE genes, and marker genes
  sce_sc3 <- SC3::sc3(sce, ks = k_vec, gene_filter = F, biology = T)
  return(sce_sc3)
}

# args
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay
celltype_key <- opt$celltype_key
k_num = opt$k_num

#1. Obtain SCE object from input_path
sce <- distribute_input(input_path = opt$input_object)
#2. Add modifications to make SC3 work
sce <- common_modif_sc3(sce)
#3. Run SC3
sce_sc3 <- run_sc3(sce)
#4. Extract metadata where the cluster annotation is
sce_sc3_metadata <- sce_sc3@colData
print(head(sce_sc3_metadata))
#5. Save metadata
write.csv(sce_sc3_metadata, file = opt$clusters, row.names = T, col.names = F)
print("SC3 Cluster annotation saved")
