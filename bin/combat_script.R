
#!/usr/bin/env Rscript
  
#libraries
library(SingleCellExperiment) #object processing
library(sva) #ComBat

args <- commandArgs(trailingOnly = TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

#input file
dataset <- readRDS(args[1])

#data
mod_data <- as.data.frame(t(as.matrix(logcounts(dataset))))
# Basic batch removal
mod0 = model.matrix(~ 1, data = mod_data)
#run ComBat
t1 = Sys.time()

assay(dataset, "corrected") <- ComBat(
  dat = t(mod_data),
  batch = dataset$dataset,
  mod = mod0,
  par.prior = TRUE,
  prior.plots = FALSE
)

t2 = Sys.time()
print(t2-t1)

#Note: an error appears when telling:Error in while (change > conv) { : missing value where TRUE/FALSE needed
#This is fixed, removing those genes with 0 variance across cells!
#This may be also because of -Inf values in logcounts

print("ComBat worked!")

saveRDS(dataset, file = args[2])

print("congratulations, this worked!!!")