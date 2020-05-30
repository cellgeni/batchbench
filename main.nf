#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.metadata = null
params.datadir  = null  
params.outdir   = null

Channel.fromPath(params.metadata)
  .flatMap{ it.readLines() }
  .view()
  .set { ch_datalist }


process get_datasets {

    tag "$datasetname"

    input:
    val datasetname from ch_datalist
    output:
    //set val(datasetname), file('*.h5ad') into ch_QC_h5ad
    set val(datasetname), file('*rds') into ch_QC_rds
    
    shell:
    '''
    #pyfile="!{params.datadir}/!{datasetname}.h5ad"
    Rfile="!{params.datadir}/!{datasetname}.rds"
    if [[ ! -e $Rfile ]]; then
      echo "Please check existence of $pyfile and $Rfile"
      false
    fi
    #ln -s $pyfile .
    ln -s $Rfile .
    '''
}


// preform QC based on min_genes expressed per cell, min_cells with expression per gene, bt_thres and ct thres 
process QC_rds{
   	//memory { 5.GB + 50.MB * Math.max(4*(1<<20), ((1<<20) + (int) ((task.attempt+1) ** 0.5 * datain.size()))).intdiv(1<<20) }
	tag "QC $datasetname"
        input:
        set val(datasetname), file(datain) from ch_QC_rds
        output:
        set val(datasetname), file('QC.*.rds') into ch_sce2h5ad, ch_R_harmony, ch_R_limma, ch_R_combat, ch_R_Seurat_v3, ch_R_mnnCorrect, ch_R_fastMNN
        set val(datasetname), val('logcounts'), file("QC.*.rds") into ch_logcounts
	shell:
        '''
        R_QC_dataset.R --input !{datain} --bt_thres 0 --ct_thres 0 --min_genes 0 --min_cells 0 --output QC.!{datasetname}.rds
        '''
}

// convert Sce objects to H5ad for the python tools
process h5ad2sce {
	tag "$convert h5ad2sce"
	input:
        set val(datasetname), file(datain) from ch_sce2h5ad
	output:
 	set val(datasetname), file("*.h5ad") into ch_bbknn, ch_scanorama
	shell:
	'''
	convert_sce2h5ad.R --input !{datain} --assay_name logcounts --output QC.!{datasetname}.h5ad
	''' 	
}

// run BBKNN method
process BBKNN {
    	publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
        tag "BBKNN $datasetname"
        input:
        set val(datasetname), file(datain) from ch_bbknn
        output:
        set val(datasetname), val('bbknn'), file('bbknn.*.h5ad') into ch_bbknn_entropy
        shell:
        '''
        bbknn_method.py --input !{datain} --output bbknn.!{datasetname}.h5ad
        '''
}

// run Scanorama method
process Scanorama {
    	publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
     	tag "Scanorama $datasetname"
     	input:
     	set val(datasetname), file(datain) from ch_scanorama
     	output:
     	set val(datasetname), val('scanorama'), file('scanorama.*.h5ad') into ch_scanorama_entropy
     	shell:
    	'''
     	scanorama_method.py --input !{datain} --output scanorama.!{datasetname}.h5ad
     	'''
 }

// merge python tool channels 
ch_bbknn_entropy.mix(ch_scanorama_entropy).into{ch_py_tools_entropy; ch_py_tools_UMAP}

// run Harmony method
process Harmony {
    	publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
	tag "Harmony $datasetname"
	input:
	set val(datasetname), file(datain) from ch_R_harmony
	output:
	set val(datasetname), val('harmony'), file('harmony.*.rds') into ch_harmony_entropy, ch_harmony_converter
	shell:
	'''
	harmony_method.R --input !{datain} --output harmony.!{datasetname}.rds
	'''
}

// run Limma method
process Limma {
    publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
    tag "Limma $datasetname"
    input:
    set val(datasetname), file(datain) from ch_R_limma
    output:
    set val(datasetname), val('limma'), file('limma.*.rds') into ch_R_limma_entropy, ch_limma_converter
    shell:
    '''
    limma_method.R --input !{datain} --output limma.!{datasetname}.rds
    '''
}

// run ComBat method
process Combat{
    publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
    tag "Combat $datasetname"
    input:
    set val(datasetname), file(datain) from ch_R_combat
    output:
    set val(datasetname), val('ComBat'), file('ComBat.*.rds') into ch_R_ComBat_entropy, ch_ComBat_converter
    shell:
    '''
    combat_method.R --input !{datain} --output ComBat.!{datasetname}.rds
    '''
}

// run Seurat 3 method
process Seurat_3_anchors{
    publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
    tag "Seurat 3 $datasetname"
    input:
    set val(datasetname), file(datain) from ch_R_Seurat_v3
    output:
    set val(datasetname), val('Seurat3'), file('Seurat3.*.rds') into ch_R_Seurat_v3_entropy, ch_Seurat_v3_converter
    shell:
    '''
    seurat3_method.R --input !{datain} --output Seurat3.!{datasetname}.rds
    '''
}

// run mnnCorrect method
process MnnCorrect{
    publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
    tag "mnnCorrect $datasetname"
    label "long_running"
    input:
    set val(datasetname), file(datain) from ch_R_mnnCorrect
    output:
    set val(datasetname), val('mnnCorrect'), file('mnnCorrect.*.rds') into ch_R_mnnCorrect_entropy, ch_mnnCorrect_converter

    shell:
    '''
    mnnCorrect_method.R --input !{datain} --output mnnCorrect.!{datasetname}.rds
    '''
}

// run fastMNN method
process fastMNN{
    publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
    tag "fastMNNt $datasetname"
    input:
    set val(datasetname), file(datain) from ch_R_fastMNN
    output:
    set val(datasetname), val('fastMNN'), file('fastMNN.*.rds') into ch_R_fastMNN_entropy, ch_fastMNN_converter
    shell:
    '''
    fastMNN_method.R --input !{datain} --output fastMNN.!{datasetname}.rds
    '''
}

//Here I'm sending all the rds corrected objets, to the rds --> h5ad converter. So that we can integrate UMAP into a single python script.
ch_harmony_converter.mix(ch_limma_converter, ch_ComBat_converter, ch_Seurat_v3_converter, ch_mnnCorrect_converter, ch_fastMNN_converter).set{ch_rds_to_h5ad_converter}

// In Future will go away!
// convert SCE to H5ad objects to compute UMAP only in Python
process rds_to_h5ad_converter {
    	tag "rds_to_h5ad_converter $datain $method $datasetname"
    	label "fast_running"
    	input:
    	set val(datasetname), val(method), file(datain) from ch_rds_to_h5ad_converter
    	output:
    	set val(datasetname), val(method), file('*.h5ad') into ch_R_tools_UMAP 
    	
    	shell:
    	'''
   	rds_h5ad_converter.R --input !{datain} --output !{method}.!{datasetname}.h5ad
    	'''
}

// merge channels that go to R entropy
ch_logcounts.mix(ch_R_limma_entropy, ch_harmony_entropy, ch_R_ComBat_entropy, ch_R_Seurat_v3_entropy, ch_R_mnnCorrect_entropy, ch_R_fastMNN_entropy).set{ch_R_entropy}

// compute Shannon entropy in R
process R_entropy {
    	publishDir "${params.outdir}/${datasetname}/entropy" 
    	tag "Entropy $datain $method $datasetname"
    	label "long_running"

    	input:
    	set val(datasetname), val(method), file(datain) from ch_R_entropy

    	output:
    	file('*.csv')
    	shell:
    	'''
    	R_entropy_after.R --input !{datain} --output entropy_!{method}.!{datasetname}.csv
    	'''
}

// compute Shannon entropy in py
process py_entropy {
    	publishDir "${params.outdir}/${datasetname}/entropy" 
    	tag "entropy (python) $datain $method $datasetname"

    	input:
    	set val(datasetname), val(method), file(datain) from ch_py_tools_entropy
    	output:
    	file('*.csv')

    	shell:
    	'''
    	entropy_py.py --input !{datain} --output_entropy  entropy_!{method}.!{datasetname}.csv 
    	'''
}

ch_py_tools_UMAP.mix(ch_R_tools_UMAP).set{ ch_UMAP }

process py_UMAP {
	//memory = { 30.GB + 30.GB * (task.attempt - 1) }
    	tag "UMAP $datain $method $datasetname"
    	publishDir "${params.outdir}/${datasetname}/UMAP" , pattern: '*.csv' 
    	input:
    	set val(datasetname), val(method), file(datain) from ch_UMAP
    	output:
    	file('*.csv')
    	//file('*.h5ad')
    	shell:
    	'''
    	#export NUMBA_CACHE_DIR=/lustre/scratch117/cellgen/cellgeni/ruben/batchbench/NUMBA_CACHE
    	Py_UMAP.py --input !{datain} --output umap.!{method}.!{datasetname}.csv  
    	'''
    //Py_UMAP.py --input !{datain} --output umap.!{method}.!{datasetname}.csv --output_h5ad corrected.!{method}.!{datasetname}.h5ad 
}
