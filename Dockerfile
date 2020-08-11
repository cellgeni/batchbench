FROM quay.io/cellgeni/batchbench:v0.26 
# Install other CRAN
RUN Rscript -e 'install.packages(c("RaceID"), dependencies = TRUE)'
