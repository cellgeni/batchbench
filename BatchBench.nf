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
    set val(datasetname), file('*.h5ad') into ch_QC_h5ad
    set val(datasetname), file('*rds') into ch_QC_rds
    
    shell:
    '''
    pyfile="!{params.datadir}/!{datasetname}.h5ad"
    Rfile="!{params.datadir}/!{datasetname}.rds"
    if [[ ! -e $pyfile || ! -e $Rfile ]]; then
      echo "Please check existence of $pyfile and $Rfile"
      false
    fi
    ln -s $pyfile .
    ln -s $Rfile .
    '''
   //# '''
   //# sourcedir="!{params.datadir}"
   //# pyfile="!{datasetname}.h5ad"
   //# Rfile="!{datasetname}.rds"
   // 
   // #if [[ ! -e "$sourcedir/$pyfile" ]]; then
   //   #R --vanilla --input "$sourcedir/$Rfile" --input_layer "logcounts" --output $pyfile  < !{baseDir}/bin/rds_h5ad_converter.R
   //   #ln -s $pyfile "$sourcedir/$pyfile"
   // #elif [[ ! -e "sourcedir/$Rfile" ]]; then
   //   #R --vanilla --input "$sourcedir/$pyfile" --output $Rfile < !{baseDir}/bin/h5ad_rds_converter.R
   //   #ln -s $Rfile "$sourcedir/$Rfile" 
   // #elif [[ ! -e $Rfile && ! -e $pyfile ]]; then
   //  # echo "Please check existence of $pyfile and $Rfile"
   //   #false
   // #fi
   //# '''   
}

process QC_h5ad {
        //tag "$datasetname $method"
        tag "$datasetname"
   	//memory { 5.GB + 100.MB * Math.max(4*(1<<20), ((1<<20) + (int) ((task.attempt+1) ** 0.5 * datain.size()))).intdiv(1<<20) }
    	publishDir "${params.outdir}/${datasetname}/", pattern: "*.txt" 
    	publishDir "${params.outdir}/${datasetname}/Corrected_objects", pattern: "*.h5ad" 
        input:
        set val(datasetname), file(datain) from ch_QC_h5ad
        output: 
        //set val(datasetname), file('QC.*.h5ad') into ch_bbknn, ch_scanorama
        set val(datasetname),  file('QC.*.h5ad') into ch_bbknn, ch_scanorama
        set val(datasetname), val("logcounts"), file('QC.*.h5ad') into ch_logcounts
       	file('QC_info.*.txt')

        shell:
        '''
	#export NUMBA_CACHE_DIR=/lustre/scratch117/cellgen/cellgeni/ruben/batchbench/NUMBA_CACHE
        Py_QC_dataset.py --input !{datain} --bt_thres 0 --ct_thres 0 --min_genes 0 --min_cells 0 --output QC.!{datasetname}.h5ad > QC_info.!{datasetname}.txt
        '''
}


//Save_logcounts
//ch_logcounts_save
//  .subscribe {
//      dataname  = it[0]
//      method    = it[1]
//      filename  = it[2]
//      filename.copyTo("${params.outdir}/${dataname}/Corrected_objects/${method}.${dataname}.h5ad")
//      }

process QC_rds{
   	//memory { 5.GB + 50.MB * Math.max(4*(1<<20), ((1<<20) + (int) ((task.attempt+1) ** 0.5 * datain.size()))).intdiv(1<<20) }
	tag "$datasetname"
        input:
        set val(datasetname), file(datain) from ch_QC_rds
        output:
        set val(datasetname), file('QC.*.rds') into ch_entropy_before, ch_distance_before, ch_umap_before, ch_R_harmony, ch_R_limma, ch_R_combat, ch_R_Seurat_v2, ch_R_Seurat_v3, ch_R_mnnCorrect, ch_R_fastMNN
        shell:
        '''
        R_QC_dataset.R --input !{datain} --bt_thres 0 --ct_thres 0 --min_genes 0 --min_cells 0 --output QC.!{datasetname}.rds
        '''
}

process BBKNN {
    	//memory { 180.MB * (((1<<20) + datain.size()).intdiv(1<<20)) }
   	//memory { 150.MB * Math.max(4*(1<<20), ((1<<20) + (int) ((task.attempt+1) ** 0.5 * datain.size()))).intdiv(1<<20) }
        tag "BBKNN (py) $datasetname"
        input:
        set val(datasetname), file(datain) from ch_bbknn
    	publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
        output:
        set val(datasetname), val('bbknn'), file('bbknn.*.h5ad') into ch_bbknn_entropy
        shell:
        '''
        #export NUMBA_CACHE_DIR=/lustre/scratch117/cellgen/cellgeni/ruben/batchbench/NUMBA_CACHE
        bbknn_script.py --input !{datain} --output bbknn.!{datasetname}.h5ad
        '''
}

process Scanorama {
     	//memory { 180.MB * (((1<<20) + datain.size()).intdiv(1<<20)) }
     	//memory { 200.MB * Math.max(4*(1<<20), ((1<<20) + (int) ((task.attempt+1) ** 0.5 * datain.size()))).intdiv(1<<20) }
     	tag "Scanorama (py) $datasetname"
     	input:
     	set val(datasetname), file(datain) from ch_scanorama
    	publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
     	output:
     	set val(datasetname), val('scanorama'), file('scanorama.*.h5ad') into ch_scanorama_entropy
     	shell:
    	'''
     	#export NUMBA_CACHE_DIR=/lustre/scratch117/cellgen/cellgeni/ruben/batchbench/NUMBA_CACHE
     	scanorama_script.py --input !{datain} --output scanorama.!{datasetname}.h5ad
     	'''
 }

ch_logcounts.mix(ch_bbknn_entropy, ch_scanorama_entropy).into{ch_py_tools_entropy; ch_py_tools_UMAP}

process Harmony {
    	//memory { 100.MB * (((1<<20) + datain.size()).intdiv(1<<20)) }
    	//memory { 100.MB * Math.max(4*(1<<20), ((1<<20) + (int) ((task.attempt+1) ** 0.5 * datain.size()))).intdiv(1<<20) }
	tag "Harmony (R) $datasetname"
	input:
	set val(datasetname), file(datain) from ch_R_harmony
    	publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
	output:
	set val(datasetname), val('harmony'), file('harmony.*.rds') into ch_harmony_entropy, ch_harmony_converter
	shell:
	'''
	harmony_script.R --input !{datain} --output harmony.!{datasetname}.rds
	'''
}

process Limma {
    //memory { 280.MB * (((1<<20) + datain.size()).intdiv(1<<20)) }
    //memory { 100.MB * Math.max(4*(1<<20), ((1<<20) + (int) ((task.attempt+1) ** 0.5 * datain.size()))).intdiv(1<<20) }
    tag "Limma (R) $datasetname"
    input:
    set val(datasetname), file(datain) from ch_R_limma
    publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
    output:
    set val(datasetname), val('limma'), file('limma.*.rds') into ch_R_limma_entropy, ch_limma_converter
    shell:
    '''
    limma_script.R --input !{datain} --output limma.!{datasetname}.rds
    '''
}

process Combat{
    //memory { 180.MB * (((1<<20) + datain.size()).intdiv(1<<20)) }
    //memory { 100.MB * Math.max(4*(1<<20), ((1<<20) + (int) ((task.attempt+1) ** 0.5 * datain.size()))).intdiv(1<<20) }
    tag "Combat (R) $datasetname"
    input:
    set val(datasetname), file(datain) from ch_R_combat
    publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
    output:
    set val(datasetname), val('ComBat'), file('ComBat.*.rds') into ch_R_ComBat_entropy, ch_ComBat_converter
    shell:
    '''
    combat_script.R --input !{datain} --output ComBat.!{datasetname}.rds
    '''
}

process Seurat_3_anchors{
    //memory { 400.MB * (((1<<20) + datain.size()).intdiv(1<<20)) }
    //memory { 100.MB * Math.max(4*(1<<20), ((1<<20) + (int) ((task.attempt+1) ** 0.5 * datain.size()))).intdiv(1<<20) }
    tag "Seurat3 (R) $datasetname"
    input:
    set val(datasetname), file(datain) from ch_R_Seurat_v3
    publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
    output:
    set val(datasetname), val('Seurat3'), file('Seurat3.*.rds') into ch_R_Seurat_v3_entropy, ch_Seurat_v3_converter
    shell:
    '''
    Seurat_v3_anchors.R --input !{datain} --output Seurat3.!{datasetname}.rds
    '''
}
process MnnCorrect{
    //memory { 400.MB * (((1<<20) + datain.size()).intdiv(1<<20)) }
    //memory { 50.MB * Math.max(4*(1<<20), ((1<<20) + (int) ((task.attempt+1) ** 0.5 * datain.size()))).intdiv(1<<20) }
    tag "mnnCorrect (R) $datasetname"
    label "long_running"
    input:
    set val(datasetname), file(datain) from ch_R_mnnCorrect
    publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
    output:
    set val(datasetname), val('mnnCorrect'), file('mnnCorrect.*.rds') into ch_R_mnnCorrect_entropy, ch_mnnCorrect_converter

    shell:
    '''
    mnnCorrect_script.R --input !{datain} --output mnnCorrect.!{datasetname}.rds
    '''
}

process fastMNN{
    //memory { 400.MB * (((1<<20) + datain.size()).intdiv(1<<20))*task.attempt }
    //memory { 100.MB * Math.max(4*(1<<20), ((1<<20) + (int) ((task.attempt+1) ** 0.5 * datain.size()))).intdiv(1<<20) }
    tag "fastMNNt (R) $datasetname"
    input:
    set val(datasetname), file(datain) from ch_R_fastMNN
    publishDir "${params.outdir}/${datasetname}/Corrected_objects" 
    output:
    set val(datasetname), val('fastMNN'), file('fastMNN.*.rds') into ch_R_fastMNN_entropy, ch_fastMNN_converter
    shell:
    '''
    fastMNN.R --input !{datain} --output fastMNN.!{datasetname}.rds
    '''
}
//Here I'm sending all the rds corrected objets, to the rds --> h5ad converter. So that we can integrate UMAP into a single python script.
ch_harmony_converter.mix(ch_limma_converter, ch_ComBat_converter, ch_Seurat_v3_converter, ch_mnnCorrect_converter, ch_fastMNN_converter).set{ch_rds_to_h5ad_converter}

process rds_to_h5ad_converter {
	//memory = { 30.GB + 30.GB * (task.attempt - 1) }
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


ch_R_limma_entropy.mix(ch_harmony_entropy, ch_R_ComBat_entropy, ch_R_Seurat_v3_entropy,  ch_R_mnnCorrect_entropy, ch_R_fastMNN_entropy).set{ch_R_entropy}

process R_entropy {
	//memory = { 30.GB + 30.GB * (task.attempt - 1) }
    	tag "entropy_after (R) $datain $method $datasetname"
    	label "long_running"

    	publishDir "${params.outdir}/${datasetname}/entropy" 

    	input:
    	set val(datasetname), val(method), file(datain) from ch_R_entropy

    	output:
    	file('*.csv')
    	shell:
    	'''
    	R_entropy_after.R --input !{datain} --output entropy_!{method}.!{datasetname}.csv
    	'''
}

process py_entropy {
	//memory = { 10.GB + 30.GB * (task.attempt - 1) }
    	tag "entropy (python) $datain $method $datasetname"
    	publishDir "${params.outdir}/${datasetname}/entropy" 

    	input:
    	set val(datasetname), val(method), file(datain) from ch_py_tools_entropy
    	output:
    	file('*.csv')

    	shell:
    	'''
    	#export NUMBA_CACHE_DIR=/lustre/scratch117/cellgen/cellgeni/ruben/batchbench/NUMBA_CACHE
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
