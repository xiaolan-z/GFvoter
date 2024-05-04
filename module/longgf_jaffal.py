import argparse
import subprocess
import os

def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_fastq_gz', type=str, required=False, default="None_file",
                        help='input_path')
    parser.add_argument('-t', action='store', dest='datatype',
                         type=str, required=True, help='datatype')
    # parser.add_argument('-po', action='store', dest='process_output_file', type=str, required=False, default="process_output_file",
    #                     help='process_output_file_path.')   
    parser.add_argument('-o', action='store', dest='output_file', type=str, required=False, default="./GFvoter_out/",
                        help='output directory.')  
    parser.add_argument('-s', action='store', dest='choice', type=str, required=False, default="real",
                        help='the source of input read.')                                                                                    
    args = parser.parse_args()
    return args

def run_commands(commands):
    for cmd in commands:
        result = subprocess.run(cmd, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            output = result.stdout
        else:
            error = result.stderr
            print("Error:", error)

args = initialization_parameters()
input_fastq_gz=args.input_fastq_gz
datatype=args.datatype
# process_output_file=args.process_output_file
output_file=args.output_file
if output_file[-1]!='/':
    output_file+='/'
choice=args.choice


# os.makedirs(process_output_file, exist_ok=True)
# os.chdir('./process_output_file')
os.makedirs(output_file+'process_output_file/', exist_ok=True)
os.chdir(output_file+'process_output_file/')

commands=[
        f"gunzip -c {input_fastq_gz} > data.fastq",
         "minimap2 -ax splice ../../hg38.fa  data.fastq > data.sam",
         "samtools sort -o data.bam data.sam",
         "samtools index data.bam"  ]
if datatype=="pacbio":
   commands1 = [
       "samtools sort -n data.bam -o sorted.bam",
       "LongGF sorted.bam ../../gencode.v22.chr_patch_hapl_scaff.annotation.gtf 100 50 100 0 0 1 > LongGF.run.on.XXX.date.log",
       "grep 'SumGF' LongGF.run.on.XXX.date.log > longgf_results.txt"
     ]
if datatype=="nanopore":
   commands1 = [
       "samtools sort -n data.bam -o sorted.bam",
       "LongGF sorted.bam ../../gencode.v22.chr_patch_hapl_scaff.annotation.gtf 100 50 100 0 0 2 > LongGF.run.on.XXX.date.log",
       "grep 'SumGF' LongGF.run.on.XXX.date.log > longgf_results.txt"
     ]

commands2=[
         f"../../JAFFA/tools/bin/bpipe run  ../../JAFFA/JAFFAL.groovy {input_fastq_gz}"
         ]

run_commands(commands)
run_commands(commands1)
run_commands(commands2)

l="longgf_results.txt"
LongGF_fusions=[]
with open(l,"r") as f:
     lines=f.readlines()
     for line in lines:
         columns=line.strip().split()
         if len(columns)>=2:
            LongGF_fusions.append(columns[1])

csvfile="jaffa_results.csv"
    # 存储第二列数据的空列表
JAFFAL_fusions = []
with open(csvfile,'r') as f:
    lines=f.readlines()
    for line in lines[1:]:
        elements=line.strip().split(',')
        if len(elements)>=2:
            JAFFAL_fusions.append(elements[1])

common_fusion=[]
for fusion in JAFFAL_fusions:
    fusion1, fusion2 = fusion.split(':')
    fusion_reversed = fusion2 + ":" + fusion1  # 融合基因前后顺序反转
    if fusion in LongGF_fusions or fusion_reversed in LongGF_fusions:
       common_fusion.append(fusion)
