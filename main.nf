#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

Channel.fromPath(params.dataset_list)
  .flatMap{ it.readLines() }
  .view()
  .set { DATALIST }

// fetch data
process get_datasets {
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 2.GB + 10.GB * (task.attempt - 1) }
 	tag "$datasetname"

    	input:
    	val datasetname from DATALIST
    	output:
    	//set val(datasetname), file('*.h5ad') into ch_QC_h5ad
    	set val(datasetname), file('*rds') into QC_RDS_CH
    
		
    	shell:
    	'''
    	#pyfile="!{params.datadir}/!{datasetname}.h5ad" #input file can be provided in rds format only
    	Rfile="!{params.data_dir}/!{datasetname}.rds"
    	if [[ ! -e $Rfile ]]; then
    	  echo "Please check existence of $Rfile"
    	  false
    	fi
    	#ln -s $pyfile .
    	ln -s $Rfile .
    	'''
}


// preform QC based on min_genes expressed per cell, min_cells with expression per gene, remove batch and cell types representing less than bt_thres and ct_thres proportion of totalcells 
process QC_rds{
    	publishDir "${params.output_dir}/${datasetname}" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 10.GB * (task.attempt - 1) }
	tag "QC $datasetname"

        input:
        set val(datasetname), file(datain) from QC_RDS_CH
        
	output:
	set val(datasetname), file('QC.*.rds') into SCE2H5AD_INPUT, HARMONY_METHOD, LIMMA_METHOD, COMBAT_METHOD, SEURAT3_METHOD, MNNCORRECT_METHOD, FASTMNN_METHOD
        set val(datasetname), val('logcounts'), file("QC.*.rds") into LOGCOUNTS_ENTROPY
	
        """
        R_QC_dataset.R\
		--input ${datain}\
		--bt_thres ${params.QC_rds.batch_thres}\
		--ct_thres ${params.QC_rds.celltype_thres}\
		--min_genes ${params.QC_rds.min_genes}\
		--min_cells ${params.QC_rds.min_cells}\
		--output QC.${datasetname}.rds
        """
}

// convert Sce objects to H5ad for the python tools
process h5ad2sce {
    	publishDir "${params.output_dir}/${datasetname}" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 2.GB + 10.GB * (task.attempt - 1) }
	tag "convert $datasetname sce2h5ad"
	
	input:
        set val(datasetname), file(datain) from SCE2H5AD_INPUT
	
	output:
 	set val(datasetname), file("*.h5ad") into BBKNN_METHOD, SCANORAMA_METHOD
	
	"""
        #export NUMBA_CACHE_DIR=/lustre/scratch117/cellgen/cellgeni/ruben/batchbench/NUMBA_CACHE
	export NUMBA_CACHE_DIR=~/NUMBA_CACHE
	convert_sce2h5ad.R\
		 --input ${datain}\
		 --assay_name ${params.assay_name}\
		 --output QC.${datasetname}.h5ad
	""" 	
}

// run BBKNN method
process BBKNN {
    	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "BBKNN $datasetname"
        
	input:
        set val(datasetname), file(datain) from BBKNN_METHOD, CLUST_SEURAT
        
	output:
        set val(datasetname), val('bbknn'), file('bbknn.*.h5ad') into BBKNN_ENTROPY
        
        """
        bbknn_method.py\
		--input ${datain}\
		--batch_key ${params.batch_key}\
		--n_pcs ${params.BBKNN.n_pcs}\
		--n_neighbours ${params.BBKNN.n_neighbours}\
		--output bbknn.${datasetname}.h5ad
        """
}

// run Scanorama method
process Scanorama {
    	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
     	tag "Scanorama $datasetname"
     	
	input:
     	set val(datasetname), file(datain) from SCANORAMA_METHOD
     	
	output:
     	set val(datasetname), val('scanorama'), file('scanorama.*.h5ad') into SCANORAMA_ENTROPY, CLUST_SC3, CLUST_SEURAT, MARKERS
     	
    	"""
     	scanorama_method.py\
		--input ${datain}\
		--batch_key ${params.batch_key}\
		--output scanorama.${datasetname}.h5ad
     	"""
 }

// merge python tool channels 
BBKNN_ENTROPY.mix(SCANORAMA_ENTROPY).into{PY_TOOLS_ENTROPY; PY_TOOLS_UMAP}

// run Harmony method
process Harmony {
    	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
	tag "Harmony $datasetname"
	
	input:
	set val(datasetname), file(datain) from HARMONY_METHOD, CLUST_SEURAT
	
	output:
	set val(datasetname), val('harmony'), file('harmony.*.rds') into HARMONY_ENTROPY, HARMONY_2H5AD
	
	"""
	harmony_method.R\
		--input_object ${datain}\
		--batch_key ${params.batch_key}\
		--assay_name ${params.assay_name}\
		--corrected_emb ${params.corrected_emb}\
		--n_pcs ${params.harmony.n_pcs}\
		--output_object harmony.${datasetname}.rds
	"""
}

// run Limma method
process Limma {
	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "Limma $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from LIMMA_METHOD, CLUST_SC3, CLUST_SEURAT, MARKERS
    	
	output:
    	set val(datasetname), val('limma'), file('limma.*.rds') into LIMMA_ENTROPY, LIMMA_2H5AD
    	
    	"""
    	limma_method.R\
		--input_object ${datain}\
		--assay_name ${params.assay_name}\
		--corrected_assay ${params.corrected_assay}\
		--output_object limma.${datasetname}.rds
    	"""
}

// run ComBat method
process Combat{
	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "Combat $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from COMBAT_METHOD, CLUST_SC3, CLUST_SEURAT, MARKERS
    	
	output:
    	set val(datasetname), val('ComBat'), file('ComBat.*.rds') into COMBAT_ENTROPY, COMBAT_2H5AD
    	"""
    	combat_method.R\
		--input_object ${datain}\
		--assay_name ${params.assay_name}\
		--corrected_assay ${params.corrected_assay}\
		--batch_key ${params.batch_key}\
		--output_object ComBat.${datasetname}.rds
    	"""
}

// run Seurat 3 method
process Seurat_3{
 	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "Seurat 3 $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from SEURAT3_METHOD, CLUST_SC3, CLUST_SEURAT, MARKERS
    	
	output:
    	set val(datasetname), val('Seurat3'), file('Seurat3.*.rds') into SEURAT3_ENTROPY, SEURAT3_2H5AD
    	
	
    	"""
    	seurat3_method.R\
		--input_object ${datain}\
		--assay_name ${params.assay_name}\
		--corrected_assay ${params.corrected_assay}\
		--batch_key ${params.batch_key}\
		--hvg_method ${params.Seurat_3.hvg_method}\
		--n_features ${params.Seurat_3.n_features}\
		--n_anchors ${params.Seurat_3.n_anchors}\
		--output_object Seurat3.${datasetname}.rds
    	"""
}

// run mnnCorrect method
process mnnCorrect{
	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "mnnCorrect $datasetname"
    	label "long_running"
    	
	input:
    	set val(datasetname), file(datain) from MNNCORRECT_METHOD, CLUST_SC3, CLUST_SEURAT, MARKERS
    	
	output:
    	set val(datasetname), val('mnnCorrect'), file('mnnCorrect.*.rds') into MNNCORRECT_ENTROPY, MNNCORRECT_H5AD

    	"""
    	mnnCorrect_method.R\
		--input_object ${datain}\
		--batch_key ${params.batch_key}\
		--k_num ${params.mnnCorrect.k}\
		--sigma ${params.mnnCorrect.sigma}\
		--svd_dim ${params.mnnCorrect.svd_dim}\
		--assay_name ${params.assay_name}\
		--corrected_assay ${params.corrected_assay}\
		--cos_norm ${params.mnnCorrect.cos_norm}\
		--output_object mnnCorrect.${datasetname}.rds
    	"""
}

// run fastMNN method
process fastMNN{
 	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "fastMNNt $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from FASTMNN_METHOD, CLUST_SEURAT 
    	
	output:
    	set val(datasetname), val('fastMNN'), file('fastMNN.*.rds') into FASTMNN_ENTROPY, FASTMNN_2H5AD
    	
    	"""
    	fastMNN_method.R\
		--input_object ${datain}\
		--batch_key ${params.batch_key}\
		--k_num ${params.fastMNN.k}\
		--n_pcs ${params.fastMNN.n_pcs}\
		--assay_name ${params.assay_name}\
		--corrected_emb ${params.corrected_emb}\
		--cos_norm ${params.fastMNN.cos_norm}\
		--output_object fastMNN.${datasetname}.rds
    	"""
}

//Send rds objets to the rds --> h5ad converter, to compute UMAP in python (where UMAP is implemented.
HARMONY_2H5AD.mix(LIMMA_2H5AD, COMBAT_2H5AD, SEURAT3_2H5AD, MNNCORRECT_H5AD, FASTMNN_2H5AD).set{SCE2H5AD_UMAP}

// convert SCE to H5ad objects to compute UMAP only in Python
process rds_to_h5ad_converter {
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "rds_to_h5ad_converter $datain $method $datasetname"
    	label "fast_running"
    	
	input:
    	set val(datasetname), val(method), file(datain) from SCE2H5AD_UMAP
    	
	output:
    	set val(datasetname), val(method), file('*.h5ad') into R_TOOLS_UMAP 
    	
    	"""
   	rds_h5ad_converter.R\
		 --input ${datain}\
		 --assay_name ${params.assay_name}\
		 --corrected_assay ${params.corrected_assay}\
		 --corrected_emb ${params.corrected_emb}\
		 --output ${method}.${datasetname}.h5ad
    	"""
}

//Merge rds objects channels to compute entropy in R 
LOGCOUNTS_ENTROPY.mix(LIMMA_ENTROPY, HARMONY_ENTROPY, COMBAT_ENTROPY, SEURAT3_ENTROPY, MNNCORRECT_ENTROPY, FASTMNN_ENTROPY).set{R_TOOLS_ENTROPY}

// compute Shannon entropy in R
process R_entropy {
    	publishDir "${params.output_dir}/${datasetname}/entropy" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 20.GB + 20.GB * (task.attempt - 1) }
    	tag "Entropy $datain $method $datasetname"
    	label "long_running"

    	input:
    	set val(datasetname), val(method), file(datain) from R_TOOLS_ENTROPY

    	output:
    	file('*.csv')
    	
    	"""
    	R_entropy.R\
		--input_object ${datain}\
		--assay_name ${params.assay_name}\
		--corrected_assay ${params.corrected_assay}\
		--corrected_emb ${params.corrected_emb}\
		--batch_key ${params.batch_key}\
		--celltype_key ${params.celltype_key}\
		--k_num ${params.R_entropy.k_num}\
		--dim_num ${params.R_entropy.dim_num}\
		--output_entropy entropy_${method}.${datasetname}.csv
    	"""
}

// compute Shannon entropy in py
process py_entropy {
    	publishDir "${params.output_dir}/${datasetname}/entropy" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "entropy (python) $datain $method $datasetname"

    	input:
    	set val(datasetname), val(method), file(datain) from PY_TOOLS_ENTROPY
    	output:
    	file('*.csv')

    	"""
    	entropy_py.py\
		--input ${datain}\
		--batch_key ${params.batch_key}\
		--celltype_key ${params.celltype_key}\
		--n_neighbours ${params.Py_entropy.n_neighbours}\
		--n_pcs ${params.Py_entropy.n_pcs}\
		--output_entropy  entropy_${method}.${datasetname}.csv 
    	"""
}

process clust_SC3{
    	publishDir "${params.output_dir}/${datasetname}/clustering" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "entropy (python) $datain $method $datasetname"

    	input:
    	set val(datasetname), val(method), file(datain) from CLUST_SC3
    	output:
    	file('*.csv')

    	"""
    	clust_SC3.R\
		--input_object ${datain}\
		--assay_name ${params.assay_name}\
		--corrected_assay ${params.corrected_assay}\
		--celltype_key ${params.celltype_key}\
		--output_clusters sc3_clusters_${method}.${datasetname}.csv 
	
    	"""
}


process clust_Seurat{
    	publishDir "${params.output_dir}/${datasetname}/clustering" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "entropy (python) $datain $method $datasetname"

    	input:
    	set val(datasetname), val(method), file(datain) from CLUST_SEURAT
    	output:
    	file('*.csv')

    	"""
    	clust_Seurat.R\
		--input_object ${datain}\
		--assay_name ${params.assay_name}\
		--corrected_assay ${params.corrected_assay}\
		--celltype_key ${params.celltype_key}\
		--output_clusters seurat_clusters_${method}.${datasetname}.csv 
	
    	"""
}

process find_markers{
    	publishDir "${params.output_dir}/${datasetname}/markers" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "entropy (python) $datain $method $datasetname"

    	input:
    	set val(datasetname), val(method), file(datain) from MARKERS
    	output:
    	file('*.csv')

    	"""
    	clust_Seurat.R\
		--input_object ${datain}\
		--assay_name ${params.assay_name}\
		--corrected_assay ${params.corrected_assay}\
		--celltype_key ${params.celltype_key}\
		--output_markers markers_${method}.${datasetname}.csv 
	
    	"""
}


PY_TOOLS_UMAP.mix(R_TOOLS_UMAP).set{ ALL_UMAP }

process py_UMAP {
    	publishDir "${params.output_dir}/${datasetname}/UMAP" , pattern: '*.csv' 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "UMAP $datain $method $datasetname"
    	
	input:
    	set val(datasetname), val(method), file(datain) from ALL_UMAP
    	
	output:
    	file('*.csv')
    	
    	"""
    	Py_UMAP.py\
		--input_object ${datain}\
		--n_neighbours ${params.Py_UMAP.n_neighbours}\
		--n_pcs ${params.Py_UMAP.n_pcs}\
		--output_umap umap.${method}.${datasetname}.csv  
    	"""
}
