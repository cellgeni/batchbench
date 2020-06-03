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

suppressPackageStartupMessages(library(Seurat))

# Scale data
scale_data <- function(seurat_obj, assays_interest, use_genes){
  all.genes <- rownames(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj, assay = corrected_assay, features = all.genes, verbose = F)
  return(seurat_obj)
}

# Dimensionality reduction
dim_red <- function(seurat_obj, assays_interest, use_genes ){
  
    seurat_obj <- RunPCA(seurat_obj, assay = corrected_assay, features = all.genes, npcs = n_pcs, 
                         reduction.name =  paste0("pca_", corrected_assay), 
                         verbose = F)
  return(seurat_obj)
}

# Graph 
neares_neighbour <- function(seurat_obj){
    base_seurat <- FindNeighbors(seurat_obj, reduction = names(seurat_obj@reductions)[1], 
                                 k.param = k_num,  
                                 dims = k_pcs, 
                                 graph.name = paste0("NN_", names(seurat_obj@reductions[1])),
                                 compute.SNN = F, verbose = F)
  return(seurat_obj)
}

# Clustering
seurat_clustering <- function(seurat_obj, algorithm){
    base_seurat <- FindClusters(seurat_obj, graph.name = names(seurat_obj@graphs)[1],  resolution = 0.8, algorithm = clust_alg , verbose = F) 
  return(seurat_obj)
}

#7. EXTRACT CLUSTERING ANNOTATION
extract_clustering_annot <- function(seurat_obj = base_seurat){
  clustering_annotation <- seurat_obj@meta.data[, startsWith(colnames(seurat_obj@meta.data), "NN_")]
  return(clustering_annotation)
}
#8. SAVE CLUSTERING ANNOTATION
save_clustering_annot <- function(clustering_annotation, save_path, algorithm, dataset_name){
  save_name <- paste0(algorithm,".", dataset_name, ".csv")
  write.csv(clustering_annotation, file = paste0(save_path, save_name))
}

base_seurat <- readRDS(opt$input_object)

# args
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay
celltype_key <- opt$celltype_key
k_num <-  opt$k_num
n_pcs <- opt$n_pcs

# Scale data of corrected expression matrices
all.genes <- rownames(base_seurat)
base_seurat <- scale_data(seurat_obj = base_seurat, assays_interest = assays_interest, use_genes = all.genes)
 
# Dimensionality reduction of corrected expression matrices
base_seurat <- dim_red(seurat_obj = base_seurat, assays_interest = assays_interest, use_genes = all.genes)

#5. Neares neighbour calculation of Dim reds
base_seurat <- neares_neighbour(seurat_obj = base_seurat)

# Louvain clustering
base_seurat <- seurat_clustering(seurat_obj = base_seurat, algorithm = 1)
# Extract Louvain annotation
clustering_annotation <- extract_clusteringa_annotation(seurat_obj = base_seurat)
# Save Louvain annotation
save_clustering_annot(clustering_annotation = clustering_annotation, save_path = opt$output_clusters, 
                      algorithm = "louvain", dataset_name = dataset_name)

# Leiden clustering
base_seurat <- seurat_clustering(seurat_obj = base_seurat, algorithm = 4)
# Extract Leiden annotation
clustering_annotation <- extract_clusteringa_annotation(seurat_obj = base_seurat)
# Save Leiden annotation
save_clustering_annot(clustering_annotation = clustering_annotation, save_path =  opt$output_clusters, 
                      algorithm = "leiden", dataset_name = dataset_name)



