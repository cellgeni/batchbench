FROM r-base:4.0.0
#
RUN apt-get update --fix-missing && apt-get install -y procps wget bzip2 ca-certificates libglib2.0-0 libxext6 libsm6 libxrender1 git mercurial subversion build-essential libcurl4-gnutls-dev libssl-dev libxml2-dev libhdf5-dev


RUN mkdir -p /root/.conda

ENV MINICONDA_VERSION=4.7.10
ENV PATH /opt/conda/bin:$PATH
RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

# python packages
ENV PYTHON_VERSION=3.7.4
RUN conda install python=${PYTHON_VERSION}
RUN pip install scipy scanpy bbknn scanorama leidenalg

# Install other CRAN
RUN Rscript -e 'install.packages(c("Seurat", "rJava", "umap","ggplot2", "ggfortify", "Rmagic", "BiocManager", "devtools","lsa","uwot", "SC3", "optparse"), dependencies = TRUE)'

# Install Bioconductor packages
RUN Rscript -e 'BiocManager::install(c("limma", "SummarizedExperiment", "SingleCellExperiment", "DropletUtils", "LoomExperiment", "Rhdf5lib", "scater", "scran", "RUVSeq", "sva", "MultiAssayExperiment", "batchelor"))'

# install github packages
RUN Rscript -e 'devtools::install_github(c("immunogenomics/harmony", "LTLA/beachmat", "cellgeni/sceasy", "mojaveazure/loomR"))'

ENV OPENBLAS_CORETYPE=nehalem
