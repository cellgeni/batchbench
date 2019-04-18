FROM quay.io/cellgeni/notebooks-base:master

USER root

# Install other CRAN
RUN Rscript -e 'install.packages("lsa")'

# install github packages
# see here for with_libpaths description:
# https://stackoverflow.com/questions/24646065/how-to-specify-lib-directory-when-installing-development-version-r-packages-from
# (do not install anything in the home directory, it will be wiped out when a volume is mounted to the docker container)
#RUN Rscript -e 'withr::with_libpaths(new = "/usr/lib/R/site-library/", devtools::install_github(c("immunogenomics/harmony")))'

# Install Bioconductor packages
#RUN Rscript -e 'BiocManager::install(c("bla-bla1", "bla-bla2"), version = "3.8")'

# Install python packages
RUN pip install umap-learn

