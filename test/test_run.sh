#!/bin/bash

SCRIPT_PATH=$(dirname "$(realpath "\$0")")
cd $SCRIPT_PATH

python ../GFvoter_quick_rzt.py -i testdata.fastq.gz -o test_out -t pacbio -s real
