#!/usr/bin/env Rscript

# Run SC3 clustering. For integrated expression matrices. 

# TODO: Set up the GENE FILTERIN

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
    c("-m", "--method"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Batch correction method the input file comes from'
  ),
  make_option(
    c("-t", "--celltype_key"),
    action = "store",
    default = "cell_type1",
    type = 'character',
    help = 'Cell type key in cell metadata'
  ),
  make_option(
    c("-b", "--biology"),
    action = "store",
    default = TRUE,
    type = 'logical',
    help = 'Wether to calculate biological features based on the identified cell clusters (DE genes, marker genes etc).'
  ),
  make_option(
    c("-o", "--output_clusters"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output dataframe with cluster annotation'
  ), 
  make_option(
    c("-r", "--output_rowdata"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output dataframe with gene-based analysis'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(require(SC3))


# FUNCTIONS #
# build custom SCE object
exp_sce <- function(dataset, assay_name){
  counts_mat <- assay(dataset, assay_name)
  col_data <- dataset@colData
  sce <- SingleCellExperiment(assays = list(logcounts = counts_mat),
                              colData = col_data)
  sce
}

# common modifications applied to SCE object
common_modif_sc3 <- function(sce){
  rowData(sce)$feature_symbol <- rownames(sce)
  counts(sce) <- logcounts(sce)
  sce
}
#run SC3 function
run_SC3 <- function(sce, celltype_key, biology){
  
  k <- length(table(sce[[celltype_key]])[table(sce[[celltype_key]]) > 0])
  k_vec <- c(k)
  print(paste0("k values on clustering:", k_vec))
  # biology = T enables calculation of DE genes, and marker genes
  sce_sc3 <- SC3::sc3(sce, ks = k_vec, gene_filter = FALSE, biology = biology)
  return(sce_sc3)
}

# args
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay
method <- opt$method
if(is.null(method) || is.na(method)){ stop("Please provide the batch correction method the input file comes from") }
celltype_key <- opt$celltype_key
biology <- opt$biology

# read input file
dataset <- readRDS(opt$input_object)
# 1. Build SC3 custom SCE object
if(method %in% c("logcounts", "Logcounts")){
  sce <- exp_sce(dataset, assay_name = assay_name)
  }else{
  sce <- exp_sce(dataset, assay_name = corrected_assay)
}
# 2. Add modifications to make SCE object required by SC3
sce <- common_modif_sc3(sce)
# 3. Run SC3
sce_sc3 <- run_SC3(sce, celltype_key, biology = biology)
# 4. Extract cluster annotation
sc3_clust_annot <- sce_sc3@colData[, grep("sc3_", colnames(sce_sc3@colData))]
print(head(sc3_clust_annot))
# 5. Save cluster annotation
write.csv(sc3_clust_annot, file = opt$output_clusters, row.names = T)
print("SC3 Cluster annotation saved")

if(biology == TRUE) {
  # Extract gene 
  gene_annot <- rowData(sce_sc3)
  # Save 
  write.csv(gene_annot, file = opt$output_rowdata, row.names = T)
  print("SC3 Gene annotation saved")
}
