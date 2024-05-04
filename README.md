# GFvoter
 A software tool for accurate detection of gene fusion in long-read transcriptome sequencing.


## Description
GFvoter is a method that accurately identifies gene fusion from long-read transcriptome sequencing data with multi-voting strategies. The five voters consist of two alignment programs, Minimap2, Winnowmap, two existing gene fusion detection methods, LongGF, JAFFAL, and a custom scoring mechanism.And both in real and simulated data, GFvoter performs very well overall in fusion detection.
The input should be a fastq.gz file (PacBio Iso-Seq, Nanopore, or mixed platform). The output is a list of confident gene fusions. By default, FusionSeeker uses hg38.fa and gencode.v22.chr_patch_hapl_scaff.annotation.gtf as reference.


## Installation
    GFvoter can be used on the linux system.
### Prerequisites
    conda install -c bioconda minimap2=2.24 samtools=1.9 winnowmap=2.03 -y
    conda install -c bioconda longgf
### Download the GFvoter package
    git clone https://github.com/xiaolan-z/GFvoter.git 
    cd ./GFvoter
## General usage
GFvoter can be applied with built-in Human reference genome (hg3.fa) and annotation (gencode.v22.chr_patch_hapl_scaff.annotation.gtf):
    python GFvoter.py  -i <path_to_directory with fastq files>/*.gz -s <real or simulate> -t <pacbio or nanopore> -o GFvoter_out/

A test dataset is available to verify successful installation:
    python GFvoter.py  -i ../../testdata.fastq.gz -s simulate -t pacbio -o GFvoter_out_test/

optional arguments:
    -i   Input **.fastq.gz file
    -t   Input read type(pacbio,nanopore)
    -s   The source of input read (real,simulate)
    -o   The directory path to save the result files of GFvoter
    -n   Minimal reads supporting an gene fusion event. Default=3  
    -l   The minimum overlap length between two aligments of one read. Default=15
    -el  The minimum exon overlap length of an alignment record against the reference genome. Default=400
    -sn  The number of secondary alignment of a read. Default=0
    -rp  The propotion of one alignment to the . Default=0.75
    -gp  The propotion of one alignment to the genome. Default=0.3

## Output files
The output directory includes:
    process_output_file/   Intermediate files during running GFvoter.py.
    reported_fusions.txt   A list of gene fusions reported by GFvoter from input read.




