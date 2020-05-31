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
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
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
    	  echo "Please check existence of $pyfile and $Rfile"
    	  false
    	fi
    	#ln -s $pyfile .
    	ln -s $Rfile .
    	'''
}


// preform QC based on min_genes expressed per cell, min_cells with expression per gene, remove batch and cell types representing less than bt_thres and ct_thres proportion of totalcells 
process QC_rds{
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 10.GB + 10.GB * (task.attempt - 1) }
	tag "QC $datasetname"

        input:
        set val(datasetname), file(datain) from QC_RDS_CH
        
	output:
	set val(datasetname), file('QC.*.rds') into SCE2H5AD_INPUT, HARMONY_METHOD, LIMMA_METHOD, COMBAT_METHOD, SEURAT3_METHOD, MNNCORRECT_METHOD, FASTMNN_METHOD
        set val(datasetname), val('logcounts'), file("QC.*.rds") into LOGCOUNTS_ENTROPY
	
	shell:
        '''
        R_QC_dataset.R --input !{datain} --bt_thres 0 --ct_thres 0 --min_genes 0 --min_cells 0 --output QC.!{datasetname}.rds
        '''
}

// convert Sce objects to H5ad for the python tools
process h5ad2sce {
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 2.GB + 10.GB * (task.attempt - 1) }
	tag "$convert h5ad2sce"
	
	input:
        set val(datasetname), file(datain) from SCE2H5AD_INPUT
	
	output:
 	set val(datasetname), file("*.h5ad") into BBKNN_METHOD, SCANORAMA_METHOD
	
	shell:
	'''
	convert_sce2h5ad.R --input !{datain} --assay_name logcounts --output QC.!{datasetname}.h5ad
	''' 	
}

// run BBKNN method
process BBKNN {
    	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
        tag "BBKNN $datasetname"
        
	input:
        set val(datasetname), file(datain) from BBKNN_METHOD
        
	output:
        set val(datasetname), val('bbknn'), file('bbknn.*.h5ad') into BBKNN_ENTROPY
        
	shell:
        '''
        bbknn_method.py --input !{datain} --output bbknn.!{datasetname}.h5ad
        '''
}

// run Scanorama method
process Scanorama {
    	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
     	tag "Scanorama $datasetname"
     	
	input:
     	set val(datasetname), file(datain) from SCANORAMA_METHOD
     	
	output:
     	set val(datasetname), val('scanorama'), file('scanorama.*.h5ad') into SCANORAMA_ENTROPY
     	
	shell:
    	'''
     	scanorama_method.py --input !{datain} --output scanorama.!{datasetname}.h5ad
     	'''
 }

// merge python tool channels 
BBKNN_ENTROPY.mix(SCANORAMA_ENTROPY).into{PY_TOOLS_ENTROPY; PY_TOOLS_UMAP}

// run Harmony method
process Harmony {
    	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
	tag "Harmony $datasetname"
	
	input:
	set val(datasetname), file(datain) from HARMONY_METHOD
	
	output:
	set val(datasetname), val('harmony'), file('harmony.*.rds') into HARMONY_ENTROPY, HARMONY_2H5AD
	
	shell:
	'''
	harmony_method.R --input !{datain} --output harmony.!{datasetname}.rds
	'''
}

// run Limma method
process Limma {
	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "Limma $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from LIMMA_METHOD
    	
	output:
    	set val(datasetname), val('limma'), file('limma.*.rds') into LIMMA_ENTROPY, LIMMA_2H5AD
    	shell:
    	'''
    	limma_method.R --input !{datain} --output limma.!{datasetname}.rds
    	'''
}

// run ComBat method
process Combat{
	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "Combat $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from COMBAT_METHOD
    	
	output:
    	set val(datasetname), val('ComBat'), file('ComBat.*.rds') into COMBAT_ENTROPY, COMBAT_2H5AD
    	
	shell:
    	'''
    	combat_method.R --input !{datain} --output ComBat.!{datasetname}.rds
    	'''
}

// run Seurat 3 method
process Seurat_3{
 	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "Seurat 3 $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from SEURAT3_METHOD
    	
	output:
    	set val(datasetname), val('Seurat3'), file('Seurat3.*.rds') into SEURAT3_ENTROPY, SEURAT3_2H5AD
    	
	shell:
    	'''
    	seurat3_method.R --input !{datain} --output Seurat3.!{datasetname}.rds
    	'''
}

// run mnnCorrect method
process MnnCorrect{
	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "mnnCorrect $datasetname"
    	label "long_running"
    	
	input:
    	set val(datasetname), file(datain) from MNNCORRECT_METHOD
    	
	output:
    	set val(datasetname), val('mnnCorrect'), file('mnnCorrect.*.rds') into MNNCORRECT_ENTROPY, MNNCORRECT_H5AD

    	shell:
    	'''
    	mnnCorrect_method.R --input !{datain} --output mnnCorrect.!{datasetname}.rds
    	'''
}

// run fastMNN method
process fastMNN{
 	publishDir "${params.output_dir}/${datasetname}/Corrected_objects" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "fastMNNt $datasetname"
    	
	input:
    	set val(datasetname), file(datain) from FASTMNN_METHOD
    	
	output:
    	set val(datasetname), val('fastMNN'), file('fastMNN.*.rds') into FASTMNN_ENTROPY, FASTMNN_2H5AD
    	
	shell:
    	'''
    	fastMNN_method.R --input !{datain} --output fastMNN.!{datasetname}.rds
    	'''
}

//Send rds objets to the rds --> h5ad converter, to compute UMAP in python (where UMAP is implemented.
HARMONY_2H5AD.mix(LIMMA_2H5AD, COMBAT_2H5AD, SEURAT3_2H5AD, MNNCORRECT_H5AD, FASTMNN_2H5AD).set{SCE2H5AD_UMAP}

// In Future will go away!
// convert SCE to H5ad objects to compute UMAP only in Python
process rds_to_h5ad_converter {
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "rds_to_h5ad_converter $datain $method $datasetname"
    	label "fast_running"
    	
	input:
    	set val(datasetname), val(method), file(datain) from SCE2H5AD_UMAP
    	
	output:
    	set val(datasetname), val(method), file('*.h5ad') into R_TOOLS_UMAP 
    	
    	shell:
    	'''
   	rds_h5ad_converter.R --input !{datain} --output !{method}.!{datasetname}.h5ad
    	'''
}

//Merge rds objects channels to compute entropy in R 
LOGCOUNTS_ENTROPY.mix(LIMMA_ENTROPY, HARMONY_ENTROPY, COMBAT_ENTROPY, SEURAT3_ENTROPY, MNNCORRECT_ENTROPY, FASTMNN_ENTROPY).set{R_TOOLS_ENTROPY}

// compute Shannon entropy in R
process R_entropy {
    	publishDir "${params.output_dir}/${datasetname}/entropy" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 20.GB + 20.GB * (task.attempt - 1) }
    	tag "Entropy $datain $method $datasetname"
    	label "long_running"

    	input:
    	set val(datasetname), val(method), file(datain) from R_TOOLS_ENTROPY

    	output:
    	file('*.csv')
    	shell:
    	'''
    	R_entropy_after.R --input !{datain} --output entropy_!{method}.!{datasetname}.csv
    	'''
}

// compute Shannon entropy in py
process py_entropy {
    	publishDir "${params.output_dir}/${datasetname}/entropy" 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "entropy (python) $datain $method $datasetname"

    	input:
    	set val(datasetname), val(method), file(datain) from PY_TOOLS_ENTROPY
    	output:
    	file('*.csv')

    	shell:
    	'''
    	entropy_py.py --input !{datain} --output_entropy  entropy_!{method}.!{datasetname}.csv 
    	'''
}

PY_TOOLS_UMAP.mix(R_TOOLS_UMAP).set{ ALL_UMAP }

process py_UMAP {
    	publishDir "${params.output_dir}/${datasetname}/UMAP" , pattern: '*.csv' 
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt <= process.maxRetries ? 'retry' : 'finish' }
	memory = { 10.GB + 20.GB * (task.attempt - 1) }
    	tag "UMAP $datain $method $datasetname"
    	
	input:
    	set val(datasetname), val(method), file(datain) from ALL_UMAP
    	
	output:
    	file('*.csv')
    	
	shell:
    	'''
    	Py_UMAP.py --input !{datain} --output umap.!{method}.!{datasetname}.csv  
    	'''
}
