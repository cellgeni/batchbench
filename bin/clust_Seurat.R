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
    c("-k", "--k_num"),
    action = "store",
    default = 30,
    type = 'integer',
    help = 'Number of nearest neighbors for the NN graph construction'
  ),
  make_option(
    c("-p", "--n_pcs"),
    action = "store",
    default = 25,
    type = 'integer',
    help = 'Number of Principal Components to perform Dimensionality Reduction'
  ),
  make_option(
    c("-o", "--louvain_clusters"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to h5ad output file'
  ),
  make_option(
    c("-q", "--leiden_clusters"),
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
  seurat_obj <- FindClusters(seurat_obj, graph.name = tail(names(dataset@graphs), 1),  resolution = 0.8, algorithm = clust_alg , verbose = F)
  seurat_obj@meta.data[["seurat_clusters"]]
}

suppressPackageStartupMessages(library(Seurat))

# args
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay
corrected_emb <- opt$corrected_emb
celltype_key <- opt$celltype_key
k_num <-  opt$k_num
n_pcs <- opt$n_pcs
method <- opt$method
if(is.null(method) || is.na(method)){ stop("Please provide the batch correction method the input file comes from") }

# read input file
dataset <- readRDS(opt$input_object)
features <- as.character(read.csv(opt$input_features, row.names =1, header = T)$x)

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
  # Write results separately
  write.csv(louvain_clusters, file = opt$louvain_clusters, row.names = colnames(dataset))
  print(paste0("Louvain clustering successfuly computed for method: ", method))
  write.csv(leiden_clusters, file = opt$leiden_clusters, row.names = colnames(dataset))
  print(paste0("Leiden clustering successfuly computed for method: ", method))
}

# Case 2: for embedding correcting methods (fastMNN, harmony)
if (method %in% c("fastMNN", "FastMNN", "Harmony", "harmony")){
  # compute NN graph over batch corrected low-d embedding
  dataset <- compute_nn_graph(seurat_obj = dataset, corrected_emb = corrected_emb, k_num = k_num)
  # find clusters Louvain 
  louvain_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 1)
  # find clusters Leiden
  leiden_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 4)
  # Write results separately
  write.csv(louvain_clusters, file = opt$louvain_clusters, row.names = colnames(dataset))
  print(paste0("Louvain clustering successfuly computed for method: ", method))
  write.csv(leiden_clusters, file = opt$leiden_clusters, row.names = colnames(dataset))
  print(paste0("Leiden clustering successfuly computed for method: ", method))
}

# Case 3: for counts matrix correcting methods except Seurat (logcounts, mnnCorrect, limma, ComBat, Scanorama)
if (method %in% c("logcounts", "Logcounts", "mnnCorrect", "mnncorrect", "limma", "Limma", "ComBat", "combat", "Scanorama", "scanorama")){
  # logcounts
  if(method %in% c("logcounts", "Logcounts")){select_assay <- assay_name }else{ select_assay <- corrected_assay }
  # subset input object by features
  dataset <- dataset[features, ]
  
  all.genes <- rownames(dataset)
  # scale data
  dataset <- scale_data(dataset, assay_name = select_assay, use_genes = all.genes)
  # dim red
  dataset <- dim_red(dataset, assay_name = select_assay, use_genes = all.genes, n_pcs = n_pcs)
  # NN Graph
  dataset <- compute_nn_graph(dataset, corrected_emb = tail(names(dataset@reductions) ,1),  k_num = k_num)
  # find clusters Louvain 
  louvain_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 1)
  # find clusters Leiden
  leiden_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 4)
  # Write results separately
  write.csv(louvain_clusters, file = opt$louvain_clusters, row.names = colnames(dataset))
  print(paste0("Louvain clustering successfuly computed for method: ", method))
  write.csv(leiden_clusters, file = opt$leiden_clusters, row.names = colnames(dataset))
  print(paste0("Leiden clustering successfuly computed for method: ", method))
}

# Case 4: for Seurat object
if (method %in% c("Seurat3", "seurat3", "Seurat")){
  all.genes <- rownames(dataset)
  # scale data
  dataset <- scale_data(dataset, assay_name = corrected_assay, use_genes = all.genes)
  # dim red
  dataset <- dim_red(dataset, assay_name = corrected_assay, use_genes = all.genes, n_pcs = n_pcs)
  # NN Graph
  dataset <- compute_nn_graph(dataset, corrected_emb = tail(names(dataset@reductions),1), k_num = k_num)
  # find clusters Louvain 
  louvain_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 1)
  # find clusters Leiden
  leiden_clusters <- seurat_clustering(seurat_obj = dataset, clust_alg = 4)
  # Write results separately
  write.csv(louvain_clusters, file = opt$louvain_clusters, row.names = colnames(dataset))
  print(paste0("Louvain clustering successfuly computed for method: ", method))
  write.csv(leiden_clusters, file = opt$leiden_clusters, row.names = colnames(dataset))
  print(paste0("Leiden clustering successfuly computed for method: ", method))
}
