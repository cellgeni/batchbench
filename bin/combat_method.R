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
        help = 'Minimum number of cells for a gene to be expressed in.'
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

suppressPackageStartupMessages(require(SingleCellExperiment))
suppressPackageStartupMessages(require(sva))

#args 
assay_name <- opt$assay_name
corrected_assay <- opt$corrected_assay 
batch_key <- opt$batch_key
# input file
dataset <- readRDS(opt$input_object)
batch_vector <- as.character(dataset[[batch_key]])
# run ComBat
mod_data <- as.data.frame(t(as.matrix(assay(dataset,assay_name))))
mod0 = model.matrix(~ 1, data = mod_data)

assay(dataset, corrected_assay) <- ComBat(
  dat = t(mod_data),
  batch = as.character(dataset$Batch),
  mod = mod0,
  par.prior = TRUE,
  prior.plots = FALSE
)
# save corrected object
saveRDS(dataset, file = opt$output_object) 
print("ComBat worked!")
