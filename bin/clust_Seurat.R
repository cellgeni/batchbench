#!/usr/bin/env Rscript

# Run Seurat Louvain and Leiden clustering algorithms

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
    c("-e", "--corrected_emb"),
    action = "store",
    default = "corrected_emb",
    type = 'character',
    help = 'Name of the fastMNN corrected low-dimensional embedding.'
  ),
  make_option(
    c("-m", "--method"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Bacth correction method the input comes from'
  ),
  make_option(
    c("-t", "--celltype_key"),
    action = "store",
    default = "cell_type1",
    type = 'character',
    help = 'Cell type key in cell metadata'
  ),
  make_option(
    c("-p", "--n_pcs"),
    action = "store",
    default = 25,
    type = 'integer',
    help = 'Number of Principal Components to perform Dimensionality Reduction'
  ),
  make_option(
    c("-k", "--k_num"),
    action = "store",
    default = 30,
    type = 'integer',
    help = 'Number of nearest neighbors for the NN graph construction'
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

## FUNCTIONS ##
#1. Scale data
scale_data <- function(seurat_obj, assay_name, use_genes){
  seurat_obj <- ScaleData(seurat_obj, assay = assay_name, features = use_genes, verbose = F)
  seurat_obj
}
#2. Dimensionality reduction
dim_red <- function(seurat_obj, assay_name, use_genes, n_pcs){
  seurat_obj <- RunPCA(seurat_obj, assay = assay_name, features = use_genes, npcs = n_pcs, 
                       reduction.name =  paste0("pca_", assay_name), 
                       verbose = F)
  seurat_obj
}
#3. NN Graph 
compute_nn_graph <- function(seurat_obj, corrected_emb, k_num){
  seurat_obj <- FindNeighbors(seurat_obj, reduction = corrected_emb,
                              k.param = k_num,
                              graph.name = paste0("NN_", names(seurat_obj@reductions[1])),
                              compute.SNN = F, verbose = F)
  return(seurat_obj)
}
#4. Clustering
seurat_clustering <- function(seurat_obj, clust_alg){
  seurat_obj <- FindClusters(seurat_obj, graph.name = names(seurat_obj@graphs)[1],  resolution = 0.8, algorithm = clust_alg , verbose = F)
  seurat_obj@meta.data[["seurat_clusters"]]
}
#5. Merge Louvain & Leiden annotations
merge_annots <- function(louvain, leiden, cell_names){
  clust_dat <- data.frame("Louvain" = louvain,  "Leiden" = leiden, row.names = cell_names)
  clust_dat
}

suppressPackageStartupMessages(library(Seurat))
# read input file
dataset <- readRDS(opt$input_object)
# args
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay
corrected_emb <- opt$corrected_emb
celltype_key <- opt$celltype_key
k_num <-  opt$k_num
n_pcs <- opt$n_pcs
method <- opt$method
if(is.null(method) || is.na(method)){ stop("Please provide the batch correction method the input file comes from") }

## EXECUTE ##
# Case 1: for graph correcting methods (BBKNN)
if (method %in% c("bbknn", "BBKNN")){
  # add graph
  graph_name <- names(dataset@graphs)
  dataset@graphs[[graph_name]] <- as.Graph(dataset@graphs[[graph_name]])
  # find clusters Louvain 
  louvain_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 1)
  # find clusters Leiden
  leiden_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 4)
  # merge Louvain & Leiden annots
  clust_dat <- merge_annots(louvain = louvain_clusters, leiden = leiden_clusters, cell_names = colnames(dataset))
  # save data
  write.csv(clust_dat, file = opt$output_clusters, row.names = T)
  print(paste0("Louvain and Leiden markers successfuly computed for method: ", method))
}

# Case 2: for embedding correcting methods (fastMNN, harmony)
if (method %in% c("fastMNN", "FastMNN", "Harmony", "harmony")){
  # compute NN graph over batch corrected low-d embedding
  dataset <- compute_nn_graph(seurat_obj = dataset, corrected_emb = corrected_emb, k_num = k_num)
  # find clusters Louvain 
  louvain_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 1)
  # find clusters Leiden
  leiden_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 1)
  # merge Louvain & Leiden annots
  clust_dat <- merge_annots(louvain = louvain_clusters, leiden = leiden_clusters, cell_names = colnames(dataset))
  # save data
  write.csv(clust_dat, file = opt$output_clusters, row.names = T)
  print(paste0("Louvain and Leiden markers successfuly computed for method: ", method))
}

# Case 2: for embedding correcting methods (fastMNN, harmony)
if (method %in% c("fastMNN", "FastMNN", "Harmony", "harmony")){
  # compute NN graph over batch corrected low-d embedding
  dataset <- compute_nn_graph(seurat_obj = dataset, corrected_emb = corrected_emb, k_num = k_num)
  # find clusters Louvain 
  louvain_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 1)
  # find clusters Leiden
  leiden_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 1)
  # merge Louvain & Leiden annots
  clust_dat <- merge_annots(louvain = louvain_clusters, leiden = leiden_clusters, cell_names = colnames(dataset))
  # save data
  write.csv(clust_dat, file = opt$output_clusters, row.names = T)
  print(paste0("Louvain and Leiden markers successfuly computed for method: ", method))
}

# Case 3: for counts matrix correcting methods except Seurat (logcounts, mnnCorrect, limma, ComBat, Scanorama)
if (method %in% c("logcounts", "Logcounts", "mnnCorrect", "mnncorrect", "limma", "Limma", "ComBat", "combat", "Scanorama", "scanorama")){
  # logcounts
  if(method %in% c("logcounts", "Logcounts")){assay_name <- assay_name }else{ assay_name <- corrected_assay }
  all.genes <- rownames(dataset)
  # scale data
  dataset <- scale_data(dataset, assay_name = assay_name, use_genes = all.genes)
  # dim red
  dataset <- dim_red(dataset, assay_name = assay_name, use_genes = all.genes, n_pcs = n_pcs)
  # NN Graph
  dataset <- compute_nn_graph(dataset, corrected_emb = names(dataset@reductions), k_num = k_num)
  # find clusters Louvain 
  louvain_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 1)
  # find clusters Leiden
  leiden_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 1)
  # merge Louvain & Leiden annots
  clust_dat <- merge_annots(louvain = louvain_clusters, leiden = leiden_clusters, cell_names = colnames(dataset))
  # save data
  write.csv(clust_dat, file = opt$output_clusters, row.names = T)
  print(paste0("Louvain and Leiden markers successfuly computed for method: ", method))
}

# Case 4: for Seurat object
if (method %in% c("Seurat3", "seurat3", "Seurat")){
  all.genes <- rownames(dataset)
  # scale data
  dataset <- scale_data(dataset, assay_name = assay_name, use_genes = all.genes)
  # dim red
  dataset <- dim_red(dataset, assay_name = assay_name, use_genes = all.genes, n_pcs = n_pcs)
  # NN Graph
  dataset <- compute_nn_graph(dataset, corrected_emb = names(dataset@reductions), k_num = k_num)
  # find clusters Louvain 
  louvain_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 1)
  # find clusters Leiden
  leiden_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 1)
  # merge Louvain & Leiden annots
  clust_dat <- merge_annots(louvain = louvain_clusters, leiden = leiden_clusters, cell_names = colnames(dataset))
  # save data
  write.csv(clust_dat, file = opt$output_clusters, row.names = T)
  print(paste0("Louvain and Leiden markers successfuly computed for method: ", method))
}
