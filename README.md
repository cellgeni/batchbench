# batchbench

BatchBench is a Nextflow workflow for running the following scRNA-Seq data batch effect correction methods:
* mnnCorrect 
* limma
* ComBat
* Seurat 3
* Scanorama
* Harmony
* FastMNN
* BBKNN

## Run BatchBench

The pieline is run through the bash script `RESUME_Batchbench`. Here the profiles, data directory and output directory are specified. Nextflow can generate an execution report, a timeline and a flowchart which output directory is specified here. 
Additionally, the `-resume` command line option allows for the continuation of a workflow execution
```
source=/home/ubuntu/BatchBench.nf

$HOME/bin/nextflow run $source \
        -profile docker,local\
        --datadir /home/ubuntu/BatchBench/data\
        --metadata /home/ubuntu/BatchBench/data/dataset_list.txt\
        --outdir /home/ubuntu/BatchBench/results\
        -with-report /home/ubuntu/BatchBench/reports/report.html\
        -with-trace /home/ubuntu/BatchBench/reports/trace.txt\
        -with-timeline /home/ubuntu/BatchBench/reports/timeline.html\
        -with-dag /home/ubuntu/BatchBench/reports/flowchart.png\
        -resume
```

Equivalently, the pipeline can be run from the command line as:
```
nextflow run BatchBench.nf -profile docker,local --datadir /home/ubuntu/BatchBench/data --metadata /home/ubuntu/BatchBench/data/dataset_list.txt --outdir /home/ubuntu/BatchBench/results -resume
```
## Input
As a input to the workflow equivalent __SingleCellExperiment rds__ object and __AnnData h5ad__ object should be present in the data directory, and the data set name should be specified in the `metadata.txt` file. 

The workflow expects the log normalized counts to be stored in the `logcounts` assay in the rds object, and in `.X` matrix in the AnnData object.

#### Object Conversion:
Usually, data will be present in one format. To convert from one format to the other we provide `converter scripts`:
- Rds to h5ad:
```
rds_h5ad_converter.R\
    --input path_to_rds\
    --output path_to_new_h5ad
```
- H5ad to rds:
```
h5ad_rds_converter.R\
    --input path_to_h5ad\
    --output path_to_new_rds
```
#### Metadata specifications
For the workflow to run, both rds and h5ad objects must specify their batch annotation and cell type annotation in the object metadata, as __Batch__ and __cell_type1__ respectively. 

- For rds object:
```
SingleCellExperiment@colData@listData[["Batch"]]
SingleCellExperiment@colData@listData[["cell_type1"]]
```
- For h5ad object:
```
AnnData.obs[["Batch"]]
AnnData.obs[["cell_type1"]]
```

## BatchBench
The workflow initially checks the presence of both rds and h5ad objects for each of the datasets listed in `dataset_list.txt`.
### 1. QC
This is a mild quality control step where: 
- Cells with `Batch` or `cell_type1` labelled as NA are removed. 
- Cells with less than `min_genes` genes are removed. 
- Genex expressed in less that `min_cells` cells are removed. 
- Batches representing less than `bt_thres` per cent of the total number of cells are removed. 
- Cell types representing less than `ct_thres` per cent of the total number of cells are removed.

The purpose of `bt_thres` and `ct_thres` is avoiding very small batche or cell populations affecting the downstream entropy calculation. 

### 2. Batch correction
Each of the batch corrected object is saved as specified in the `publisDir` parameter. This way one can check how the different methods correct the input data. 

### 3. Outputs
#### 3.1 Compute Shannon Entropy
Shannon entropy is an indicator of the degree of mixing of the data. 
We calculate: 
- __batch entropy__: to evaluate the degree of alignment of the batches after the correction. 
- __cell type entropy__: to examine if different cell populations are being mixed. 
    
An ideal batch correction should output a high batch entropy , and low cell type entropy. Indicating alignment of batches while maintaining different cell types separate. 
 
Each of the batch corrected objects both entropies are computed, stored in a CSV file and saved as specified in the `publisDir` parameter. 

#### 3.2 Compute UMAP

We may want to visualize the batch correction performed by each of the methods this is why __UMAP__ (Uniform Manifold Approximation and Projection) coordinates are computed. 

To ensure reproducibility in the coordinates calculation, all corrected objects are converted to h5ad and then UMAP is computed thorugh scanpy function `sc.tl.umap`.

As output, a CSV file with the 2 first UMAP components is saved as specified in the `publisDir` parameter. 

#### 3.3 Nextflow reports

If specified when running, Nextflow can generate various reports about the execution: 
- HTML execution report
![alt text](https://www.nextflow.io/docs/latest/_images/report-summary-min.png) 
![alt text](https://www.nextflow.io/docs/latest/_images/report-resource-cpu.png)

- Trace report

Contains some useful information about each process executed in your pipeline script, including: submission time, start time, completion time, cpu and memory used.


- Timeline report
![alt text](https://www.nextflow.io/docs/latest/_images/timeline-min.png)

- DAG visualization
![alt text](https://www.nextflow.io/docs/latest/_images/dag.png)





[![Docker Repository on Quay](https://quay.io/repository/cellgeni/batchbench/status "Docker Repository on Quay")](https://quay.io/repository/cellgeni/batchbench)
