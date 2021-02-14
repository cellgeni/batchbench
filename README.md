# BatchBench: flexible comparison of batch correction methods for single-cell RNA-seq

BatchBench is a modular and flexible pipeline for comparing batch correction methods for single-cell RNA-seq data. It performs __batch correction__ considering eight popular methods for scRNA-seq data as well as a set of downstream analysis for assessing their performance including __calculation of mixing entropy__ , __clustering analysis__, __marker genes__, and __UMAP embedding generation__. 

__Publication__: 

Ruben Chazarra-Gil, Stijn van Dongen, Vladimir Yu Kiselev, Martin Hemberg, [Flexible comparison of batch correction methods for single-cell RNA-seq using BatchBench]( https://doi.org/10.1093/nar/gkab004), Nucleic Acids Research, 2021;, gkab004, https://doi.org/10.1093/nar/gkab004 

## Run BatchBench

The pieline is run through the bash script `run_batchbench.sh`. 

In this script you can specify the source directory for running the pipeline (__source_dir__), and the nextflow report directory (__report_dir__) where nextflow execution reporst will be saved. Additionally, the `-resume` option allows for the continuation of a workflow execution. 

## Input
Input to Batchbench must be an __rds__ R object of class __SingleCellExperiment__. 

To enable full funcitonallity the `nextflow.config` file parameters must be modified to you own case: 
1. __data_dir__: Directory containing the data. 
2. __dataset_list__: Path to a text file containing a list of the datasets to be included (including the file base name, without suffix).
3. __batch_key__: Object metadata column containing the batch identity of each cell. 
4. __celltype_key__: Object metadata column containing the cell type identity of each cell. __Note__: in abscence of a cell type identity of cells, a numeric vector containing the cell cluster annotation can also be used. This cluster annotaion can be obtained from section ## 3. Clustering.
5. __assay_name__: Name of the assay to perform batch correction. Default value is `logcounts`.


## Workflow structure

### 1. data QC
BatchBench initially performs a filtering of the input dataset. This can be customized according to the parameters in the __QC_rds__ process of the `nextflow.config` file. This step performs:

- Cell filtering: Cells including less than `min_genes` features are excluded. Additionally, cells with `batch_key` or `celltype_key` labelled as NA are removed. 
- Feature filtering: features present in less that `min_cells` cells are excluded. 
- Batch filtering: Batches representing less than `bt_thres` of the total number of cells are removed. 
- Cell type filtering: Cell types representing less than `ct_thres` per cent of the total number of cells are removed.

In order to disable the data QC step, the previous arguments can be stet to `0` in the `nextflow.config` file.

### 2. Batch correction 

The batch effect removal tools considered in Batchbench are: 

* [mnnCorrect](https://www.nature.com/articles/nbt.4091)
* [limma](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6687398/)
* [ComBat](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3307112/)
* [Seurat 3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6687398/) 
* [Scanorama](https://www.nature.com/articles/s41587-019-0113-3)
* [Harmony](https://www.nature.com/articles/s41592-019-0619-0)
* FastMNN
* [BBKNN](https://academic.oup.com/bioinformatics/article/36/3/964/5545955)

Each of the batch corrected methods can be enabled or disabled from the `nextflow.config` file.

The batch corrected output objects are saved in the `Corrected_objects/` subdirectory within the results directory.

### 3. Mixing Entropy

We compute Shannon entropy as an indicator of the degree of mixing of the data before and after batch correction. The batch corrected output is converted into a cell k-Nearest Neighbour graph, and then for each cell entropy is calculated considering its Batch and Cell type identity. This results in: 
- __Batch Entropy__: quantifies the degree of merging of the batches before and after the correction. 
- __Cell type Entropy__: quantifies the degree of mixing of different cell populations before and after the correction.  
    
An ideal batch correction should result a high batch entropy and low cell type entropy, indicating alignment of batches while maintaining different cell types separate. 

The output is a csv file with the Batch Entropy and Cell type Entropy for each cell: 
```
"","batch_entropy","celltype_entropy"
"human1_lib1.final_cell_0001",0.819448371872804,0.318144667641655
"human1_lib1.final_cell_0002",0.932960189010999,0.333112501389628
"human1_lib1.final_cell_0003",0.819448371872804,0.318144667641655
"human1_lib1.final_cell_0004",0.781266804259214,0.312771784077802
 ```
Entropy files for each of the batch correction methods selected, as well as for the uncorrected data are published in the `Entropy/` subdirectory within the results directory.

Parameters for the mixing entropy computation can be edited in the __entropy__ entry of the `nextflow.config` file. 

### 4. Clustering analysis

To assess the effect of batch correction on clustering analysis, we apply five unsupervised clustering methods to the batch corrected data:  
1. [__Leiden__](https://www.nature.com/articles/s41598-019-41695-z)
2. [__Louvain__](https://iopscience.iop.org/article/10.1088/1742-5468/2008/10/P10008) 
3. [__SC3__](https://www.nature.com/articles/nmeth.4236) 
4. [__RaceID__](https://www.nature.com/articles/nature14966) 
5. [__Hierarchical clustering__]

The cell cluster annotation for each of the clustering algorithms for each of the batch correction methods is stored in a csv file located in the `Clustering/` subdirectory within the results directory: 
```
"","cluster"
"human1_lib1.final_cell_0001","7"
"human1_lib1.final_cell_0002","7"
"human1_lib1.final_cell_0003","7"
"human1_lib1.final_cell_0004","7"
```

#### 4.1 Feature subsetting for clustering analysis

We perform our clustering analysis establishing five fractions of features: 0.05, 0.1, 0.2, 0.5, and 1.0 (all of the features) by ranking genes descendingly by their coefficient of variation. 

This approach considerably increases the running time of the pipeline, since the clustering step is run 5 times. In order to avoid feature selection prior to clustering, one can modify the `PROP_GENES` channel in line 68 of `main.nf` file to only include 100% of input features like: `PROP_GENES = Channel.from(1)`. Equivalently, if you want to consider different feature fractions from those default, one can do: `PROP_GENES = Channel.from(0.1, 0.2, 0.3)`. 

__Note___ this feature selection won't applied to BBKNN, Harmony and FastMNN batch correction methods since they do not operate in the gene space.

### 5. Marker gene analysis

To further assess the effect of batch correction over downstream analysis, we perform a marker gene analysis. We compute marker genes for each of the cell types in the dataset considering: the whole dataset and each of the batches composing the dataset. 

Our underlying hypothesis is that if the batch correction process has approached the batches, the marker genes found in each of these should be more similar to the ones found in the whole dataset. In our case, this similarity was calculated with the Jaccard Index metric. 

The output is an rds object with the structure: `merged_dataset`, `Batch_A`, `Batch_B`, ..., `Batch_N`, where each of these is a list of marker gene dataframes for each of the cell types of the dataset (as defined by `celltype_key`. __Note:__ if a cell cluster annotation vector is provided in `celltype_key` the marker genes will be computed equivalently, on a per-cluster basis.

This output can be found under the `Markers/` subdirectory within the results directory.


#### 6. UMAP Embedding

In order to visualize the batch correction we generate a  __UMAP embedding__ with the output of each of the methods. The batch corrected outputs are Python and R objects, this is why, to ensure reproducibility in the UMAP coordinates calculation, all objects are converted to __h5ad__ format (Python) and then UMAP is computed thorugh scanpy function `sc.tl.umap`.

The output is CSV file with the UMAP coordinates of each cell, and can be found in the `UMAP/` subdirectory within the results directory. 
```
,UMAP_1,UMAP_2
human1_lib1.final_cell_0001,-2.1425686,3.838006
human1_lib1.final_cell_0002,-2.4662907,3.6042838
human1_lib1.final_cell_0003,-1.5784191,3.1411119
human1_lib1.final_cell_0004,-2.2281802,3.3420012
```

### 7. Object conversion 

As a result of the multiple object type conversions ocurred during the pipeline, dataset objects for each of the batch correction outputs can be found format: 
1. __SingleCellExperiment__ (rds object, R)
2. __Seurat__ (rds object, R)
3. __AnnData__ (h5ad object, Python)

These objects can be found in the `Converted_objects/` subdirectory within the results directory. 

### 8. Nextflow reports

Nextflow generates a series of execution reports which contain useful information: 

- __HTML report__: a single document which includes many useful metrics about a workflow execution including __resource usage__ or __execution time__ for each process. 
- __Trace report__: an execution tracing file that contains some useful information about each process executed in your pipeline script, including: submission time, start time, completion time, cpu and memory used. 
- __Timeline report__:  HTML timeline for all processes executed in your pipeline

Thes files are located in the __reports/__ directory.  

