#/bin/bash

set -euo pipefail

unset http_proxy
unset https_proxy

source_dir="/home/ubuntu/BatchBench/runbb/git_batchbench"
data_dir="${source_dir}/test_data"
output_dir="${source_dir}/results"
report_dir="${output_dir}/reports"
$HOME/bin/nextflow run "${source_dir}/main.nf"\
        -profile docker,local\
        --datadir ${data_dir}\
        --metadata "${data_dir}/dataset_list.txt"\
        --outdir ${output_dir}\
        -with-report "${report_dir}/report.html"\
        -with-trace "${report_dir}/trace.txt"\
        -with-timeline "${report_dir}/timeline.html"\
        -with-dag "${report_dir}/flowchart.png"\
        #-resume
