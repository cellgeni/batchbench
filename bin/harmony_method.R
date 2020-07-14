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
    c("-b", "--batch_key"),
    action = "store",
    default = "Batch",
    type = 'character',
    help = 'Batch key in cell metadata'
  ),
  make_option(
    c("-a", "--assay_name"),
    action = "store",
    default = "logcounts",
    type = 'character',
    help = 'Counts assay to add to the h5ad object'
  ),
  make_option(
    c("-e", "--corrected_emb"),
    action = "store",
    default = "corrected_emb",
    type = 'character',
    help = 'Name of the fastMNN corrected low-dimensional embedding.'
  ),
  make_option(
    c("-n", "--n_pcs"),
    action = "store",
    default = 25,
    type = 'integer',
    help = 'Number of PCs to perform dimensionality reduction'
  ),
  make_option(
    c("-t", "--theta"),
    action = "store",
    default = 2,
    type = 'integer',
    help = 'Diversity clustering penalty parameter. Default theta=2. theta=0 does not encourage any diversity. Larger values of theta result in more diverse clusters.'
  ),
  make_option(
    c("-o", "--output_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to h5ad output file'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(scater))
suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(harmony))
# args
assay_name <- opt$assay_name
batch_key <- opt$batch_key
corrected_emb <- opt$corrected_emb
n_pcs <- opt$n_pcs
theta <- opt$theta

# read input object
dataset <- readRDS(opt$input_object)
# run PCA
dataset <- runPCA(dataset, exprs_values = assay_name, ncomponents = n_pcs)
pca <- dataset@reducedDims@listData[["PCA"]]
# cell batch label vector
batch_vector <- as.character(dataset[[batch_key]])
# run Harmony
dataset@reducedDims@listData[[corrected_emb]] <- HarmonyMatrix(data_mat=pca, meta_data=colData(dataset), vars_use = batch_key, theta=theta, do_pca = F)
# Add rownames to corrected embedding 
rownames(dataset@reducedDims@listData[[corrected_emb]]) <- colnames(dataset)
# save object with corrected embedding
saveRDS(dataset, opt$output_object)
print("Harmony worked!")
