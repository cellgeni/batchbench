FROM quay.io/cellgeni/cellgeni-jupyter:16.07-01

USER root

# Install other CRAN
RUN Rscript -e 'install.packages(c("Seurat", "rJava", "umap","ggplot2", "ggfortify", "Rmagic", "BiocManager", "devtools","lsa","uwot" ), dependencies = TRUE)'

# Install Bioconductor packages
RUN Rscript -e 'BiocManager::install(c("limma", "SingleCellExperiment", "Rhdf5lib", "scater", "scran", "RUVSeq", "sva", "MultiAssayExperiment", "SummarizedExperiment", "batchelor"))'

# install github packages
RUN Rscript -e 'devtools::install_github(c("immunogenomics/harmony", "LTLA/beachmat", "MarioniLab/DropletUtils", "cellgeni/sceasy", "mojaveazure/loomR"))'

# Python packages
pip install umap-learn

