FROM quay.io/cellgeni/batchbench:v0.26
# Re-install RaceID
Rscript -e 'install.packages(c("RaceID"), dependencies = TRUE)'
#Install SC3 via GitHub to solve issues
install.packages("devtools")
devtools::install_github("hemberg-lab/SC3") 
