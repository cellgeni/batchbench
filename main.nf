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
if(params.QC_rds.run == "True"){
process QC_rds{
    	publishDir "${params.output_dir}/${datasetname}", mode: 'copy' 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 10.GB * (task.attempt - 1) }
	tag "QC $datasetname"

        input:
        set val(datasetname), file(datain) from QC_RDS_CH
        
	output:
	set val(datasetname), file('QC.*.rds') into SCE2H5AD_INPUT, HARMONY_METHOD, LIMMA_METHOD, COMBAT_METHOD, SEURAT3_METHOD, MNNCORRECT_METHOD, FASTMNN_METHOD
        set val(datasetname), val('logcounts'), val('exp_matrix'), file("QC.*.rds") into LOGCOUNTS_ENTROPY, LOGCOUNTS_UMAP, LOGCOUNTS_CLUST_SC3, LOGCOUNTS_2SEURAT, LOGCOUNTS_MARKERS
	
        """
        QC_data.R\
		--input_object ${datain}\
		--bt_thres ${params.QC_rds.batch_thres}\
		--ct_thres ${params.QC_rds.celltype_thres}\
		--min_genes ${params.QC_rds.min_genes}\
		--min_cells ${params.QC_rds.min_cells}\
		--output_object QC.${datasetname}.rds
        """
	}
}
// Conv_1. Convert Sce objects to H5ad for the python tools
process h5ad2sce {
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 2.GB + 10.GB * (task.attempt - 1) }
	tag "convert $datasetname sce2h5ad"
	
	input:
        set val(datasetname), file(datain) from SCE2H5AD_INPUT
	
	output:
 	set val(datasetname), file("*.h5ad") into BBKNN_METHOD, SCANORAMA_METHOD
	
	"""
	export NUMBA_CACHE_DIR=~/NUMBA_CACHE
	convert_sce2h5ad.R\
		 --input ${datain}\
		 --assay_name ${params.assay_name}\
		 --output QC.${datasetname}.h5ad
	""" 	
}

// run BBKNN method
if(params.BBKNN.run == "True"){
process BBKNN {
    	publishDir "${params.output_dir}/${datasetname}/Corrected_objects", mode: 'copy' 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "BBKNN $datasetname"
        
	input:
        set val(datasetname), file(datain) from BBKNN_METHOD 
        
	output:
        set val(datasetname), val('bbknn'), val('graph'), file('bbknn.*.h5ad') into BBKNN_2SEURAT, BBKNN_2SCE, BBKNN_UMAP
        
        """
        bbknn_method.py\
		--input ${datain}\
		--batch_key ${params.batch_key}\
		--n_pcs ${params.BBKNN.n_pcs}\
		--n_neighbours ${params.BBKNN.n_neighbours}\
		--output bbknn.${datasetname}.h5ad
        """
	}
} else {
	BBKNN_2SEURAT = Channel.empty()
	BBKNN_2SCE = Channel.empty()
	BBKNN_UMAP = Channel.empty()
}

// run Scanorama method
if(params.scanorama.run == "True"){
process Scanorama {
    	publishDir "${params.output_dir}/${datasetname}/Corrected_objects", mode: 'copy' 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
     	tag "Scanorama $datasetname"
     	
	input:
     	set val(datasetname), file(datain) from SCANORAMA_METHOD
     	
	output:
     	set val(datasetname), val('scanorama'), val('exp_matrix'), file('scanorama.*.h5ad') into SCANORAMA_2SCE, SCANORAMA_2SEURAT, SCANORAMA_UMAP 
     	
    	"""
     	scanorama_method.py\
		--input ${datain}\
		--batch_key ${params.batch_key}\
		--output scanorama.${datasetname}.h5ad
     	"""
	}
} else {
	SCANORAMA_2SCE = Channel.empty()
	SCANORAMA_2SEURAT = Channel.empty()
	SCANORAMA_UMAP = Channel.empty()
}

PY_TOOLS_2SEURAT = BBKNN_2SEURAT.mix(SCANORAMA_2SEURAT)

// run Harmony method
if(params.harmony.run == "True"){
process Harmony {
    	publishDir "${params.output_dir}/${datasetname}/Corrected_objects", mode: 'copy' 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
	tag "Harmony $datasetname"
	
	input:
	set val(datasetname), file(datain) from HARMONY_METHOD
	
	output:
	set val(datasetname), val('harmony'), val('embedding'), file('harmony.*.rds') into HARMONY_ENTROPY, HARMONY_UMAP, HARMONY_2SEURAT, HARMONY_2H5AD
	
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
} else {
	HARMONY_ENTROPY = Channel.empty()
	HARMONY_UMAP = Channel.empty()
	HARMONY_2SEURAT = Channel.empty()
	HARMONY_2H5AD = Channel.empty()
}

// run Limma method
if(params.Limma.run == "True"){
process Limma {
	publishDir "${params.output_dir}/${datasetname}/Corrected_objects", mode: 'copy' 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "Limma $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from LIMMA_METHOD
    	
	output:
    	set val(datasetname), val('limma'), val('exp_matrix'), file('limma.*.rds') into LIMMA_ENTROPY, LIMMA_UMAP,  LIMMA_CLUST_SC3, LIMMA_2H5AD, LIMMA_2SEURAT 
    	
    	"""
    	limma_method.R\
		--input_object ${datain}\
		--assay_name ${params.assay_name}\
		--corrected_assay ${params.corrected_assay}\
		--output_object limma.${datasetname}.rds
    	"""
	}
} else {
	LIMMA_ENTROPY = Channel.empty()
	LIMMA_UMAP = Channel.empty()
	LIMMA_CLUST_SC3 = Channel.empty()
	LIMMA_2H5AD = Channel.empty()
	LIMMA_2SEURAT = Channel.empty()
}

// run ComBat method
if(params.Combat.run == "True"){
process Combat{
	publishDir "${params.output_dir}/${datasetname}/Corrected_objects", mode: 'copy' 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "Combat $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from COMBAT_METHOD
    	
	output:
    	set val(datasetname), val('ComBat'), val('exp_matrix'), file('ComBat.*.rds') into COMBAT_ENTROPY, COMBAT_UMAP, COMBAT_CLUST_SC3, COMBAT_2H5AD, COMBAT_2SEURAT 
    	"""
    	combat_method.R\
		--input_object ${datain}\
		--assay_name ${params.assay_name}\
		--corrected_assay ${params.corrected_assay}\
		--batch_key ${params.batch_key}\
		--output_object ComBat.${datasetname}.rds
    	"""
	}
} else {
	COMBAT_ENTROPY = Channel.empty()
	COMBAT_UMAP = Channel.empty()
	COMBAT_CLUST_SC3 = Channel.empty()
	COMBAT_2H5AD = Channel.empty()
	COMBAT_2SEURAT = Channel.empty()
}

// run Seurat 3 method
if(params.Seurat_3.run == "True"){
process Seurat_3{
 	publishDir "${params.output_dir}/${datasetname}/Corrected_objects", mode: 'copy' 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "Seurat 3 $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from SEURAT3_METHOD
    	
	output:
    	set val(datasetname), val('Seurat3'), val('exp_matrix'), file('Seurat3.*.rds') into SEURAT3_2H5AD, SEURAT3_2SCE, SEURAT3_CLUST_SEURAT, SEURAT3_MARKERS 
    	
	
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
} else {
	SEURAT_ENTROPY = Channel.empty()
	SEURAT3_2H5AD = Channel.empty()
	SEURAT_2SCE = Channel.empty()
	SEURAT3_CLUST_SEURAT = Channel.empty()
	SEURAT3_MARKERS = Channel.empty()
}

// run mnnCorrect method
if(params.mnnCorrect.run == "True"){
process mnnCorrect{
	publishDir "${params.output_dir}/${datasetname}/Corrected_objects", mode: 'copy' 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "mnnCorrect $datasetname"
    	label "long_running"
    	
	input:
    	set val(datasetname), file(datain) from MNNCORRECT_METHOD
    	
	output:
    	set val(datasetname), val('mnnCorrect'), val('exp_matrix'), file('mnnCorrect.*.rds') into MNNCORRECT_ENTROPY, MNNCORRECT_UMAP, MNNCORRECT_2H5AD, MNNCORRECT_2SEURAT, MNNCORRECT_CLUST_SC3 

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
} else {
	MNNCORRECT_ENTROPY = Channel.empty()
	MNNCORRECT_UMAP = Channel.empty()
	MNNCORRECT_2H5AD = Channel.empty()
	MNNCORRECT_2SEURAT = Channel.empty()
	MNNCORRECT_CLUST_SC3 = Channel.empty()
}

// run fastMNN method
if(params.fastMNN.run == "True"){
process fastMNN{
 	publishDir "${params.output_dir}/${datasetname}/Corrected_objects", mode: 'copy' 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "fastMNNt $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from FASTMNN_METHOD
    	
	output:
    	set val(datasetname), val('fastMNN'), val('embedding'), file('fastMNN.*.rds') into FASTMNN_ENTROPY, FASTMNN_UMAP, FASTMNN_2H5AD, FASTMNN_2SEURAT 
    	
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
} else {
	FASTMNN_ENTROPY = Channel.empty()
	FASTMNN_UMAP = Channel.empty()
	FASTMNN_2H5AD = Channel.empty()
	FASTMNN_2SEURAT = Channel.empty()
}

// Conv_1. Convert H5AD objects to SCE 
process conv_h5ad2sce {
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 2.GB + 10.GB * (task.attempt - 1) }
	tag "convert $datasetname sce2h5ad"
	
	input:
        set val(datasetname), val(method), val(space_corrected), file(datain) from BBKNN_2SCE  
        set val(datasetname), val(method), val(space_corrected), file(datain) from SCANORAMA_2SCE 
	
	output:
 	set val(datasetname), val(method), val(space_corrected), file('*.rds') into PY_METHODS_ENTROPY, PY_METHODS_CLUST_SC3 
	
	"""
	convert_h5ad2sce.R\
		 --input ${datain}\
		 --corrected_assay ${params.corrected_assay}\
		 --method ${method}\
		 --output_object ${method}.${datasetname}.rds
	""" 	
	}

// Only Scanorama method can be clustered by SC3
SCANORAMA_CLUST_SC3  = PY_METHODS_CLUST_SC3.filter{ it[2] == "scanorama" }

// Conv_2. Convert SEURAT object to SCE
process conv_seurat2sce {
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 2.GB + 10.GB * (task.attempt - 1) }
	tag "convert $datasetname seurat2sce"
	
	input:
        set val(datasetname), val(method), val(space_corrected), file(datain) from SEURAT3_2SCE 
	
	output:
 	set val(datasetname), val(method), val(space_corrected), file('*.rds') into SEURAT_ENTROPY, SEURAT_CLUST_SC3, SEURAT_UMAP
	
	"""
	convert_seurat2sce.R\
		 --input ${datain}\
		 --corrected_assay ${params.corrected_assay}\
		 --output_object sce_obj.${method}.${datasetname}.rds
	""" 	
	}


// Conv_3. Convert H5AD object to SEURAT
process conv_h5ad2seurat{
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
     	tag "h5ad2seurat $datasetname"
     	
	input:
     	set val(datasetname), val(method), val(space_corrected), file(datain) from PY_TOOLS_2SEURAT 
     	
	output:
     	set val(datasetname), val(method), val(space_corrected), file('*.rds') into PY_METHODS_CLUST_SEURAT, PY_METHODS_MARKERS 
	"""
	convert_h5ad2seurat.R\
		--input ${datain}\
		--corrected_assay ${params.corrected_assay}\
		--method ${method}\
		--output_object seurat_obj.${method}.${datasetname}.rds
	""" 
	}


// Only Scanorama method can be used to compute marker genes 
SCANORAMA_MARKERS  = PY_METHODS_MARKERS.filter{ it[2] == "scanorama" }

//Merge SCE object channels to convert to Seurat
LOGCOUNTS_2SEURAT.mix(HARMONY_2SEURAT, LIMMA_2SEURAT, COMBAT_2SEURAT, MNNCORRECT_2SEURAT, FASTMNN_2SEURAT).set{CONV_SCE2SEURAT}

//Conv_4. Convert SCE to Seurat 
process conv_sce2seurat{
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
     	tag "h5ad2sce $datasetname"
     	
	input:
	set val(datasetname), val(method), val(space_corrected), file(datain) from CONV_SCE2SEURAT 
     	
	output:
     	set val(datasetname), val(method), val(space_corrected), file('*.rds') into SCE_CLUST_SEURAT, SCE_MARKERS
	"""
	convert_sce2seurat.R\
		--input ${datain}\
		--assay_name ${params.assay_name}\
		--corrected_assay ${params.corrected_assay}\
		--method ${method}\
		--output_object seurat_obj.${method}.${datasetname}.rds
	""" 
	}

// filter out Harmony and fastMNN (low-D embedding) from MARKERS Channel
SCE_MARKERS_FILT  = SCE_MARKERS.filter{ it[2] != "harmony" || it[2] != "fastMNN" }

// Merge input channels for entropy
LOGCOUNTS_ENTROPY.mix(HARMONY_ENTROPY, LIMMA_ENTROPY, COMBAT_ENTROPY, MNNCORRECT_ENTROPY, FASTMNN_ENTROPY, PY_METHODS_ENTROPY,  SEURAT_ENTROPY).set{ENTROPY}

// compute entropy 
if(params.entropy.run == "True"){

// calculate Shannon entropy 
process entropy {
 	publishDir "${params.output_dir}/${datasetname}/entropy", mode: 'copy' 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 2.GB + 10.GB * (task.attempt - 1) }
	tag "compute entropy $datain"
	
	input:
        set val(datasetname), val(method), val(space_corrected), file(datain) from ENTROPY 
	
	output:
 	file("*.csv") 
	
	"""
	R_entropy.R\
		--input_object ${datain}\
		--assay_name ${params.assay_name}\
		--corrected_assay ${params.corrected_assay}\
		--corrected_emb ${params.corrected_emb}\
		--method ${method}\
		--batch_key ${params.batch_key}\
		--celltype_key ${params.celltype_key}\
		--k_num ${params.entropy.k_num}\
		--dim_num ${params.entropy.dim_num}\
		--output_entropy entropy.${method}.${datasetname}.csv
	""" 	
	}
}

//Merge all input channels to SC3 clustering
LOGCOUNTS_CLUST_SC3.mix(LIMMA_CLUST_SC3, COMBAT_CLUST_SC3, MNNCORRECT_CLUST_SC3, SCANORAMA_CLUST_SC3, SEURAT_CLUST_SC3).set{ CLUST_SC3 }

// run SC3 clustering
if(params.clust_SC3.run == "True"){

process clust_SC3{
    	publishDir "${params.output_dir}/${datasetname}/clustering", mode: 'copy' 
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
		--method ${method}\
		--celltype_key ${params.celltype_key}\
		--biology ${biology}\
		--output_clusters SC3_clusters_${method}.${datasetname}.csv\
		--output_rowdata SC3_features_${method}.${datasetname}.csv 
	
    	"""
	}
}


//Merge all Seurat clustering input channels
SEURAT3_CLUST_SEURAT.mix(PY_METHODS_CLUST_SEURAT, SCE_CLUST_SEURAT).set{ CLUST_SEURAT }

// run Seurat clustering
if(params.clust_seurat.run == "True"){
process clust_Seurat{
    	publishDir "${params.output_dir}/${datasetname}/clustering", mode: 'copy' 
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
		--corrected_emb ${params.corrected_emb}\
		--method ${method}\
		--celltype_key ${params.celltype_key}\
		--output_clusters louvain_leiden_clusters_${method}.${datasetname}.csv 
	
    	"""
	}
}

// Merge all Seurat FindMarkers input channel 
SEURAT3_MARKERS.mix(SCANORAMA_MARKERS, SCE_MARKERS_FILT).set{ MARKERS }

// run marker genes
if(params.find_markers.run == "True"){
process find_markers{
    	publishDir "${params.output_dir}/${datasetname}/markers", mode: 'copy' 
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
}


//Merge all input channels to UMAP
LOGCOUNTS_UMAP.mix(HARMONY_UMAP, LIMMA_UMAP, COMBAT_UMAP, SEURAT_UMAP, MNNCORRECT_UMAP, FASTMNN_UMAP).set{SCE_UMAP}
// run UMAP

if(params.Py_UMAP.run == "True"){
// Conv_5 Convert RDS (SCE and Seurat) to H5ad objects to compute UMAP only in Python
process rds_to_h5ad_converter {
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'ignore' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "rds_to_h5ad_converter $datain $method $datasetname"
    	label "fast_running"
    	
	input:
    	set val(datasetname), val(method), file(datain) from SCE_UMAP
    	
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

//Merge Python and converted R outputs for UMAP computation
BBKNN_UMAP.mix(SCANORAMA_UMAP, R_TOOLS_UMAP).set{ ALL_UMAP }

process py_UMAP {
    	publishDir "${params.output_dir}/${datasetname}/UMAP" , mode: 'copy', pattern: '*.csv' 
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
}
