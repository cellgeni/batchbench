# BatchBench: flexible comparison of batch correction methods for single-cell RNA-seq

BatchBench is a modular and flexible pipeline for comparing batch correction methods for single-cell RNA-seq data. It includes eight popular batch correction methods as well as a set of downstream analysis for assessing their performance including __calculation of mixing entropy , __clustering analysis__, __marker genes__, and __UMAP embedding generation__. 

The batch effect removal tools considered in Batchbench are (link to publication): 

* [mnnCorrect](https://www.nature.com/articles/nbt.4091)
* [limma](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6687398/)
* [ComBat](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3307112/)
* [Seurat 3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6687398/) 
* [Scanorama](https://www.nature.com/articles/s41587-019-0113-3)
* [Harmony](https://www.nature.com/articles/s41592-019-0619-0)
* FastMNN
* [BBKNN](https://academic.oup.com/bioinformatics/article/36/3/964/5545955)

## Run BatchBench

The pieline is run through the bash script `run_batchbench.sh`. 

In this script you can specify the source directory for running the pipeline (__source_dir__), the results output directory (__output_dir__), and the nextflow report directory (__report_dir__). Additionally, the `-resume` option allows for the continuation of a workflow execution. 

## Input
Input to Batchbench must be an __rds__ R object of class __SingleCellExperiment__. 

To enable full funcitonallity the `nextflow.config` file parameters must be modified to you own case: 
1. __data_dir__: Directory containing the data. 
2. __dataset_list__: Path to a text file containing a list of the datasets to be included (including the file base name, without preffix).
3. __batch_key__: Object metadata column containing the batch identity of each cell. 
4. __celltype_key__: Object metadata column containing the cell type identity of each cell. 
5. __assay_name__: Name of the assay to perform batch correction. Default value is `logcounts`.


## Workflow structure

### 1. data QC
BatchBench initially performs a filtering of the input dataset. This can be customized according to the parameters in the __QC_rds__ process of the `nextflow.config` file.

- Cell filtering: Cells including less than `min_genes` features are excluded. Cells with `batch_key` or `celltype_key` labelled as NA are removed. 
- Feature filtering: features present in less that `min_cells` cells are excluded. 
- Batch filtering: Batches representing less than `bt_thres` of the total number of cells are removed. 
- Cell type filtering: Cell types representing less than `ct_thres` per cent of the total number of cells are removed.

In order to disable the data QC step, the previous arguments can be stet to `0` in the `nextflow.config` file.

### 2. Batch correction methods

Each of the batch corrected methods can be enabled or disabled from the `nextflow.config` file.

The batch corrected output objects are saved in the `Corrected_objects` subdirectory within the results directory.

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
