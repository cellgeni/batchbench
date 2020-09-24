#!/usr/bin/env Rscript

# RaceID clustering

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to SCE object input file'
  ),
  make_option(
    c("-f", "--input_features"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Features names subsetted by coefficient of variation'
  ),
  make_option(
    c("-a", "--assay_name"),
    action = "store",
    default = "logcounts",
    type = 'character',
    help = 'Uncorrected assay name.'
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
    help = 'Bacth correction method the input comes from'
  ),
  make_option(
    c("-d", "--dist_metric"),
    action = "store",
    default = "pearson",
    type = 'character',
    help = 'Distance metric'
  ),
  make_option(
    c("-o", "--output_clusters"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to csv output file'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

# args
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay
method <- opt$method
dist_metric <- opt$dist_metric

# Functions
## Remove items with zero variance
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
## Remove negative values vrom corrected expression matrix
rm_neg_values <- function(dataset, assay){
  assay(dataset, assay)[assay(dataset, assay) < 0 ] <- 0 
  dataset
}

suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(RaceID))


# Read input file
dataset <- readRDS(opt$input_object)
features <- as.character(read.csv(opt$input_features, row.names =1, header = T)$x)
# Subset input object by features
dataset <- dataset[features, ]
# Select assay to cluster depending on method
if( method %in% c("Logcounts", "logcounts")){ clust_assay <- assay_name 
  }else{ clust_assay <- corrected_assay }
## 1. Remove items with zero variance (this introduces NA values in correlation matrix and clustering is aborted!)
# Remove cells with zero variance
dataset <- rm_zero_var(dataset, assay = clust_assay, axis = 1)
# Remove genes with zero variance
dataset <- rm_zero_var(dataset, assay = clust_assay, axis = 2)
# 2. Remove negative expression values if present
dataset <- rm_neg_values(dataset, assay= clust_assay)
# 3. Coerce corrected matrix to 'dgCMatrix' class
mat <- as(assay(dataset, clust_assay), "dgCMatrix")
# 3. Build RaceID object
sc <- SCseq(mat)
# 3. Add features to object (to allow successful execution)
# 3.1 Set expression data (already filtered and normalized) as filtered data (@ndata)
sc@ndata <- sc@expdata
# 3.2 Add total transcript counts
sc@counts <- colSums(as.matrix(sc@expdata))
# 3.3 Add feature names (for distance calculation)
sc@genes <- rownames(sc@ndata)
# 4. Compute distance matrix
sc <- compdist(sc, metric=dist_metric, FSelect = FALSE)
# 5. Cluster. FUNcluster arg must be one of: "kmedoids", "kmeans", "hclust".
# Run all possible clustering types (FUNcluster arg)
sc_kmedoids <- clustexp(sc, clustnr = 30, FUNcluster = "kmedoids", verbose = FALSE)
sc_kmeans <- clustexp(sc, clustnr = 30, FUNcluster = "kmeans", verbose = FALSE)
sc_hclust <- clustexp(sc, clustnr = 30, FUNcluster = "hclust", verbose = FALSE)
# 6. Save cluster annotation
write.csv(data.frame("kmedoids.RaceID" = sc_kmedoids@cluster[["kpart"]]), file = paste0("kmedoids_", opt$output_clusters))
write.csv(data.frame("kmeans.RaceID" = sc_kmeans@cluster[["kpart"]]), file = paste0("kmeans_", opt$output_clusters))
write.csv(data.frame("hclust.RaceID" = sc_hclust@cluster[["kpart"]]), file = paste0("hclust_", opt$output_clusters))

print("RaceID clustering annotations saved")
