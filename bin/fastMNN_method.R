#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
print("A")
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
        c("-k", "--k_num"),
        action = "store",
        default = "30",
        type = 'integer',
        help = 'number of nearest neighbors to consider when identifying MNNs.'
    ),
    make_option(
        c("-p", "--n_pcs"),
        action = "store",
        default = 25,
        type = 'integer',
        help = 'Number of PCs to perform dimensionality reduction'
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
        c("-n", "--cos_norm"),
        action = "store",
        default = TRUE,
        type = 'logical',
        help = ' should cosine normalization be performed on the input data prior to PCA.' 
    ),
    make_option(
        c("-o", "--output_object"),
        action = "store",
        default = NA,
        type = 'character',
        help = 'Path to rds output file'
    )
)


opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(batchelor))

# args
batch_key <- opt$batch_key
assay_name <- opt$assay_name
k_num <- opt$k_num
n_pcs <- opt$n_pcs
cos_norm <- opt$cos_norm
corrected_emb <- opt$corrected_emb
# input file
dataset <- readRDS(opt$input_object)
# batch cell label vector
batch_vector <- as.character(dataset[[opt$batch_key]])
batch_names <- unique(batch_vector)
N_batches <- length(batch_names)
print(batch_names)
print(N_batches)
# split object by batches
batch_list <- lapply(1:N_batches, function(x) {assay(dataset[, batch_vector == batch_names[x]], assay_name)})
# run fastMNN
correction <- do.call('fastMNN', c(batch_list, c(k = k_num, d = n_pcs,  cos.norm = cos_norm,  pc.input = F)))
# attach the batch corrected low_d embedding to reducedDims of the SCE object
dataset@reducedDims@listData[[corrected_emb]] <- correction@reducedDims@listData[["corrected"]]
# add colnames to embedding
colnames(dataset@reducedDims@listData[[corrected_emb]]) <- paste0("PC_", c(1:ncol(dataset@reducedDims@listData[[corrected_emb]])))
# order cells of embedding as in expression matrix 
dataset@reducedDims@listData[[corrected_emb]] <- dataset@reducedDims@listData[[corrected_emb]][colnames(assay(dataset, assay_name)), ]
# save corrected object 
saveRDS(dataset, file = opt$output_object) 
print("congratulations, fastMNN worked!")
