import argparse
import subprocess
import os
from collections import Counter

def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_fastq_gz', type=str, required=False, default="None_file",
                        help='input_path')
    parser.add_argument('-t', action='store', dest='datatype',
                         type=str, required=True, help='datatype')
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
output_file=args.output_file
if output_file[-1]!='/':
    output_file+='/'
choice=args.choice


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
LongGF_gene=[]
with open(l,"r") as f:
     lines=f.readlines()
     selected_columns=[1,2,3,4]
     for line in lines:
         elements=line.strip().split()
         if len(elements)>=2:
            selected_elements = [elements[i] for i in selected_columns]
            element1=selected_elements[0]
            element2=selected_elements[1]
            element3=selected_elements[2].split(":")[0]
            element4=selected_elements[2].split(":")[1]
            element5=selected_elements[3].split(":")[0]
            element6=selected_elements[3].split(":")[1]
            fusion=[element1,element2,element3,element4,element5,element6]
            LongGF_fusions.append(fusion)     
            LongGF_gene.append(element1)
csvfile="jaffa_results.csv"
JAFFAL_fusions = []
JAFFAL_gene=[]
with open(csvfile,'r') as f:
    lines=f.readlines()
    selected_columns1=[1,10,2,3,5,6]
    for line in lines[1:]:
        elements=line.strip().split(',')
        if len(elements)>=2:
            selected_elements = [elements[i] for i in selected_columns1]
            JAFFAL_fusions.append(selected_elements)
            JAFFAL_gene.append(selected_elements[0])

def find_duplicates(lst):
    counts = Counter(lst)
    duplicates = [item for item, count in counts.items() if count > 1]
    return duplicates
def nonrepeart(fusion_geneset,geneinfo):
   duplicates=find_duplicates(fusion_geneset)
   fusions_nonrepreat=[]
   if len(duplicates)>0:
      for fusion in duplicates:
          infom=[]
          read_num=0
          breakpoint_sum1=0
          breakpoint_sum2=0
          for item in geneinfo:
              if fusion==item[0]:
                 infom.append(item)
              else:
                 fusions_nonrepreat.append(item)
          for item in infom:
              read_num+=int(item[1])
              breakpoint_sum1+=int(item[3])
              breakpoint_sum2+=int(item[5])
              chrom1=item[2]
              chrom2=item[4]
          breakpoint1=breakpoint_sum1/len(infom)
          breakpoint2=breakpoint_sum2/len(infom)
          fusion_info=[fusion,read_num,chrom1,breakpoint1,chrom2,breakpoint2]
          fusions_nonrepreat.append(fusion_info)
   else:
        fusions_nonrepreat=geneinfo
   return fusions_nonrepreat
LongGF_results=nonrepeart(LongGF_gene,LongGF_fusions)
JAFFAL_results=nonrepeart(JAFFAL_gene,JAFFAL_fusions)


