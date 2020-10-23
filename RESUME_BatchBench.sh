#/bin/bash

set -euo pipefail

unset http_proxy
unset https_proxy

source_dir="/home/ubuntu/batchbench"
output_dir="${source_dir}/results"
report_dir="${output_dir}/reports"
/usr/local/bin/nextflow run "${source_dir}/main.nf"\
        -profile docker,local\
        -with-report "${report_dir}/report.html"\
        -with-trace "${report_dir}/trace.txt"\
        -with-timeline "${report_dir}/timeline.html"\
        -resume
