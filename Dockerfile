FROM quay.io/cellgeni/batchbench:v0.14
# Install SC3 via GitHub to solve issues
Rscript -e 'devtools::install_github("hemberg-lab/SC3")'
# 
Rscript -e 'install.packages(c("RaceID"))'
# 
pip install leidenalg 
