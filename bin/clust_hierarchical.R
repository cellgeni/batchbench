#!/usr/bin/env Rscript

# Hierarchical Clustering on batch corrected matrices
# ref: https://scrnaseq-course.cog.sanger.ac.uk/website/biological-analysis.html 
# NOTE: we only subset by features in for methods correcting count matrix (bbknn, harmony or fastMNN outpu lack feature space or is reduced to PCA space)

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
    help = 'Bacth correction method the input comes from'
  ),
  make_option(
    c("-e", "--corrected_emb"),
    action = "store",
    default = "corrected_emb",
    type = 'character',
    help = 'Name of the fastMNN corrected low-dimensional embedding.'
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

# FUNCTIONS #

# distance matrix
dist_matrix <- function(dat){
  dist_graph <- as.dist((1 - cor(dat, method = "pearson")))
  dist_graph
}
# hierarchical clustering
hierarch_clust <- function(dist_mat){
  hc <- hclust(dist_mat, method = "ward.D2")
  hc
}
# function to select k
select_k <- function(dat, hc){
  # desired number of clusters containing only 1 cell
  num.singleton <- 0
  kk <- 1
  for (i in 2:dim(dat)[2]) {
    clusters <- cutree(hc, k = i)
    clustersizes <- as.data.frame(table(clusters))
    singleton.clusters <- which(clustersizes$Freq < 2)
    if (length(singleton.clusters) <= num.singleton) { kk <- i } else { break; }
  }
  kk
}
# extract cluster df
clust_df <- function(hc, kk, cell_names){
  clusters <- cutree(hc, k = kk)
  clust_df <- data.frame("hclust" = clusters, row.names = cell_names)
  clust_df
}

##  EXECUTE  ##

# args
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay
corrected_emb <- opt$corrected_emb
method <- opt$method
if(is.null(method) || is.na(method)){ stop("Please provide the batch correction method the input file comes from") }

suppressPackageStartupMessages(require(SingleCellExperiment))
# Function 
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
# read input files
dataset <- readRDS(opt$input_object)
features <- as.character(read.csv(opt$input_features, row.names =1, header = T)$x)

# for graph correcting methods (BBKNN)
if (method %in% c("bbknn", "BBKNN")){
  # distance matrix
  dist_mat <- assay(dataset, "distances")
  # hierarchical clustering
  hc <- hierarch_clust(dist_mat = as.dist(dist_mat))
  # select optimal k value
  #kk <- select_k(dat = dist_mat, hc = hc)
  # set k value to number of cell type sin original dataset 
  kk <- length(table(as.character(dataset$cell_type1)))
  # extract clusters
  clusters <- clust_df(hc = hc, kk = kk, cell_names = colnames(dataset))
  # save clusters
  write.csv(clusters, opt$output_clusters)
  print(paste0("Hclust clustering successfuly computed for method: ", method))
}


# for embedding correcting methods (fastMNN, harmony)
if (method %in% c("fastMNN", "FastMNN", "Harmony", "harmony")){
  low_d_emb <- reducedDim(dataset, corrected_emb)
  # distance matrix (tranpose)
  dist_mat <- dist_matrix(dat = t(low_d_emb))
  # hierarchical clustering
  hc <- hierarch_clust(dist_mat = dist_mat)
  # select optimal k value
  #kk <- select_k(dat = t(low_d_emb), hc = hc)
  # set k value to number of cell type sin original dataset 
  kk <- length(table(as.character(dataset$cell_type1)))
  # extract clusters
  clusters <- clust_df(hc = hc, kk = kk, cell_names = colnames(dataset))
  # save clusters
  write.csv(clusters, opt$output_clusters)
  print(paste0("Hclust clustering successfuly computed for method: ", method))
}

# for counts matrix correcting methods except Seurat (logcounts, mnnCorrect, limma, ComBat, Scanorama)
if (method %in% c("logcounts", "Logcounts", "mnnCorrect", "mnncorrect", "limma", "Limma", "ComBat", "combat", "Scanorama", "scanorama", "Seurat3", "seurat3", "Seurat")){
  # subset input object by features
  dataset <- dataset[features, ]
  # extract counts matrix
  if(method %in% c("logcounts", "Logcounts")){ clust_assay <- assay_name } else { clust_assay <- corrected_assay }
 ## 1. Remove items with zero variance (this introduces NA values in correlation matrix and clustering is aborted!)
  # Remove cells with zero variance
  dataset <- rm_zero_var(dataset, assay = clust_assay, axis = 1)
  # Remove genes with zero variance
  dataset <- rm_zero_var(dataset, assay = clust_assay, axis = 2) 
  # distance matrix
  counts_mat <- as.matrix(assay(dataset, clust_assay))
  dist_mat <- dist_matrix(dat = counts_mat)
  # hierarchical clustering
  hc <- hierarch_clust(dist_mat = dist_mat)
  # select optimal k value
  #kk <- select_k(dat = counts_mat, hc = hc)
  # set k value to number of cell type sin original dataset 
  kk <- length(table(as.character(dataset$cell_type1)))
  # extract clusters usters <- clust_df(hc = hc, kk = kk, cell_names = colnames(dataset))
  clusters <- clust_df(hc = hc, kk = kk, cell_names = colnames(dataset))
  write.csv(clusters, opt$output_clusters)
  print(paste0("Hclust clustering successfuly computed for method: ", method))
}
