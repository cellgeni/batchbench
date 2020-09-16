FROM quay.io/cellgeni/batchbench:v0.27
#Install SC3 via GitHub to solve issues
install.packages("devtools")
devtools::install_github("hemberg-lab/SC3") 
