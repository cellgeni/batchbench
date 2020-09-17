FROM quay.io/cellgeni/batchbench:v0.26
# Re-install RaceID
RUN Rscript -e 'install.packages(c("RaceID"), dependencies = TRUE)'
# Re-install SC3 from GitHub
RUN Rscript -e 'install.packages("devtools")'
RUN Rscript -e 'devtools::install_github("hemberg-lab/SC3")'
