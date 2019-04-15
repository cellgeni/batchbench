#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

params.metadata = null
params.datadir  = null
params.outdir   = "New_results_2"

// awesome comment


Channel.fromPath(params.metadata)
  .flatMap{ it.readLines() }
  .view()
  .set { ch_datalist }


process get_datasets {

    tag "$datasetname"

    input:
    val datasetname from ch_datalist

    output:
    //set val(datasetname), file('*.h5ad') into ch_bbknn, ch_scanorama
    //set val(datasetname), file('*.rds') into ch_R_harmony, ch_R_limma, ch_R_combat, ch_R_multiCCA, ch_R_mnnCorrect
    set val(datasetname), file('*.h5ad') into ch_QC_h5ad
    set val(datasetname), file('*rds') into ch_QC_rds
    //set val(datasetname), file('*.rds') into ch_entropy_before
    //set val(datasetname), file('*.rds') into ch_distance_before
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
}
process QC_h5ad {
        tag "$datasetname"
        input:
        set val(datasetname), file(datain) from ch_QC_h5ad
        output: 
        set val(datasetname), file('QC.*.h5ad') into ch_bbknn, ch_scanorama 
        shell:
        '''
        Py_QC_dataset_v2.py --input !{datain} --batch_threshold 10 --cell_type_threshold 1 --output QC.!{datasetname}.h5ad
        '''
}

process QC_rds {
        tag "$datasetname"
        input:
        set val(datasetname), file(datain) from ch_QC_rds
        output:
        set val(datasetname), file('QC.*.rds') into ch_entropy_before, ch_distance_before, ch_R_harmony, ch_R_limma, ch_R_combat, ch_R_multiCCA, ch_R_mnnCorrect
        shell:
        '''
        R_QC_dataset_v2.R --input !{datain} --batch_threshold 10 --cell_type_threshold 1 --output QC.!{datasetname}.rds
        '''
}

process py_bbknn {
        tag "bbknn (py) $datasetname"
        input:
        set val(datasetname), file(datain) from ch_bbknn
        output:
        set val(datasetname), val('bbknn'), file('bbknn.*.h5ad') into ch_bbknn_entropy
        shell:
        '''
        bbknn_script.py --input !{datain} --output bbknn.!{datasetname}.h5ad
        '''
}

process py_scanorama {

    tag "scanorama (py) $datasetname"

    input:
    set val(datasetname), file(datain) from ch_scanorama

    output:
    set val(datasetname), val('scanorama'), file('scanorama.*.h5ad') into ch_scanorama_entropy

    shell:
    '''
    scanorama_script.py --input !{datain} --output scanorama.!{datasetname}.h5ad
    '''
}

ch_bbknn_entropy.mix(ch_scanorama_entropy).into{ch_py_entropy; ch_py_distances}


process R_harmony {
	
	  tag "harmony (R) $datasetname"

	  input:
	  set val(datasetname), file(datain) from ch_R_harmony

	  output:
	  set val(datasetname), val('harmony'), file('harmony.*.rds') into ch_harmony_entropy

	  shell:
	  '''
	  harmony_script.R --input !{datain} --output harmony.!{datasetname}.rds
	  '''
}

process R_limma {

    tag "limma (R) $datasetname"

    input:
    set val(datasetname), file(datain) from ch_R_limma

    output:
    set val(datasetname), val('limma'), file('limma.*.rds') into ch_R_limma_entropy

    shell:
    '''
    limma_script.R --input !{datain} --output limma.!{datasetname}.rds
    '''
}

process R_combat{

    tag "combat (R) $datasetname"

    input:
    set val(datasetname), file(datain) from ch_R_combat

    output:
    set val(datasetname), val('ComBat'), file('ComBat.*.rds') into ch_R_ComBat_entropy

    shell:
    '''
    combat_script.R --input !{datain} --output ComBat.!{datasetname}.rds
    '''
}


process R_multiCCA{

    tag "multiCCA (R) $datasetname"
    label "long_running"

    input:
    set val(datasetname), file(datain) from ch_R_multiCCA

    output:
    set val(datasetname), val('multiCCA'), file('multiCCA.*.rds') into ch_R_multiCCA_entropy

    shell:
    '''
    multiCCA_script.R --input !{datain} --output multiCCA.!{datasetname}.rds
    '''
}


process R_mnnCorrect{
    tag "mnnCorrect (R) $datasetname"
    label "long_running"

    input:
    set val(datasetname), file(datain) from ch_R_mnnCorrect

    output:
    set val(datasetname), val('mnnCorrect'), file('mnnCorrect.*.rds') into ch_R_mnnCorrect_entropy

    shell:
    '''
    mnnCorrect_script.R --input !{datain} --output mnnCorrect.!{datasetname}.rds
    '''
                    
}

ch_R_limma_entropy.mix(ch_harmony_entropy, ch_R_ComBat_entropy, ch_R_multiCCA_entropy, ch_R_mnnCorrect_entropy).into{ch_R_entropy; ch_R_distances}

process R_entropy {

    tag "entropy_after (R) $datain $method $datasetname"
    label "long_running"

    publishDir "${params.outdir}/${datasetname}/entropy", mode: 'copy'

    input:
    set val(datasetname), val(method), file(datain) from ch_R_entropy

    output:
    file('*.epy')

    shell:
    '''
    R_entropy_after_v2.R --input !{datain} --output entropy_!{method}.!{datasetname}.epy
    '''
}

process py_entropy {

    tag "entropy (python) $datain $method $datasetname"

    publishDir "${params.outdir}/${datasetname}/entropy", mode: 'copy'

    input:
    set val(datasetname), val(method), file(datain) from ch_py_entropy

    output:
    file('*.epy')

    shell:
    '''
    entropy_py_v2.py --input !{datain} --output entropy_!{method}.!{datasetname}.epy
    '''
}

process py_distances {
  
    tag "distances (python) $datain $method $datasetname"

    publishDir "${params.outdir}/${datasetname}/distances", mode: 'copy'

    input:
    set val(datasetname), val(method), file(datain) from ch_py_distances

    output:
    file('*.epy')

    shell:
    '''
    Py_distance_after_v2.py --input !{datain} --output distances_!{method}.!{datasetname}.epy
    '''
}

process R_entropy_before {
    
    tag "entropy_before (R) $datain $datasetname"
    
    publishDir "${params.outdir}/${datasetname}/entropy", mode: 'copy'

    input:
    set val(datasetname), file(datain) from ch_entropy_before

    output:
    file('*.epy')

    shell:
    '''
    R_entropy_before_v2.R --input !{datain} --output entropy_before_.!{datasetname}.epy
    '''
}

process R_distance_before {

    tag "distance_before (R) $datain $datasetname"

    publishDir "${params.outdir}/${datasetname}/distances", mode: 'copy'

    input:
    set val(datasetname), file(datain) from ch_distance_before

    output:
    file('*.epy')

    shell:
    '''
    distances_before_v2.R --input !{datain} --output distances_before_.!{datasetname}.epy
    '''
}

process R_distance_after {

    tag "distance_after (R) $datain $method $datasetname"

    publishDir "${params.outdir}/${datasetname}/distances", mode: 'copy'

    input:
    set val(datasetname), val(method), file(datain) from ch_R_distances

    output:
    file('*.epy')

    shell:
    '''
    distances_after_v2.R --input !{datain} --output distances_!{method}.!{datasetname}.epy
    '''
}


