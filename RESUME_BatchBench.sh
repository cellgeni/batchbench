#/bin/bash

set -euo pipefail

unset http_proxy
unset https_proxy

source=/home/ubuntu/BatchBench.nf

$HOME/bin/nextflow run $source \
        -profile docker,local\
        --datadir /home/ubuntu/BatchBench/data\
        --metadata /home/ubuntu/BatchBench/data/dataset_list.txt \
        --outdir /home/ubuntu/BatchBench/results\
        -with-report /home/ubuntu/BatchBench/reports/report.html\
        -with-trace /home/ubuntu/BatchBench/reports/trace.txt\
        -with-timeline /home/ubuntu/BatchBench/reports/timeline.html\
        -with-dag  /home/ubuntu/BatchBench/reports/flowchart.png\
        -resume