import argparse
import subprocess
import os
import datetime
from collections import Counter
import subprocess
from concurrent.futures import ThreadPoolExecutor

def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_fastq_gz', type=str, required=True,
                        help='Indicates a path to the input **.fastq.gz file')
    parser.add_argument('-t', action='store', dest='datatype',type=str, required=True, 
                        help='Indicates the type of the input data, pacbio or nanopore')
    parser.add_argument('-s', action='store', dest='choice', type=str, required=True, 
                        help='The input data is, real or simulate') 
    parser.add_argument('-o', action='store', dest='output_file', type=str, required=False, default="./GFvoter_out/",
                        help='Indicates a path to the output directory.')
    parser.add_argument('-n', action='store', dest='min_sup_read', type=int, required=False, default=1,
                        help='Indicates the minimum number of supporting reads of each reported fusion.')
    parser.add_argument('-l', action='store', dest='min_ovlp_len', type=int, required=False, default=15,
                        help='Indicates the minimum overlap length between two aligments of one read.')
    parser.add_argument('-el', action='store', dest='min_exonovlp_len', type=int, required=False, default=400,
                        help='Indicates the minimum exon overlap length of an alignment record against the reference genome.')
    parser.add_argument('-sn', action='store', dest='secondary_alignment', type=int, required=False, default=0,
                        help='Indicates the number of secondary alignment of a read.') 
    parser.add_argument('-rp', action='store', dest='read_proportion', type=float, required=False, default=0.75,
                        help='Indicates the propotion of one alignment to the read.')
    parser.add_argument('-gp', action='store', dest='gene_proportion', type=float, required=False, default=0.3,
                        help='Indicates the propotion of one alignment to the genome.')
    parser.add_argument('-sp', action='store', dest='sup_read', type=int, required=False, default=2,
                        help='Indicates the minimum number of supporting reads of each candidate in Scoring process.')                                  
    parser.add_argument('-score', action='store', dest='min_score', type=int, required=False, default=400,
                        help='Indicates the score threshold set during the scoring process.')          
    parser.add_argument('-poll', action='store', dest='min_poll', type=int, required=False, default=6,
                        help='Indicates the minimum number of votes for each reported fusion.')
    parser.add_argument('-ground_truth', action='store', dest='ground_truth', type=str, required=False, default="a",
                        help='Indicates a path to the custom known_fusions.')                                                                                                  
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

def run_command(command):
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"command '{command}' failed to execute. Errors:{result.stderr}")
    return result


def run_commands_parallel(commands1, commands2):
    with ThreadPoolExecutor(max_workers=2) as executor:
        future1 = executor.submit(run_command_list, commands1)
        future2 = executor.submit(run_command_list, commands2)

        future1.result()
        future2.result()


def run_command_list(commands):
    for command in commands:
        run_command(command)

GFvoter_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
GFvoter_path=os.path.normpath(GFvoter_path)

args = initialization_parameters()
input_fastq_gz=args.input_fastq_gz
datatype=args.datatype
output_file=args.output_file
if output_file[-1]!='/':
    output_file+='/'
choice=args.choice

start_time = datetime.datetime.now()

print(f"---------------------- Starting pineline at {start_time} ----------------------------")

os.makedirs(output_file+'process_output_file/', exist_ok=True)
os.chdir(output_file+'process_output_file/')
process_path=os.getcwd()
repetitive_k15=GFvoter_path+"/doc/repetitive_k15.txt"

commands1=[
        f"gunzip -c {input_fastq_gz} > data.fastq",
        f"minimap2 -ax splice {GFvoter_path}/hg38.fa  data.fastq > data.sam",
        "samtools sort -o data.bam data.sam",
        "samtools index data.bam"  
        ]

commands2=[
         f"winnowmap -W {repetitive_k15} -ax splice {GFvoter_path}/hg38.fa {input_fastq_gz} > winnowmap.sam",
          "samtools view -bS winnowmap.sam > winnowmap.bam",
          "samtools sort winnowmap.bam -o winnowmap.sorted.bam",
          "samtools index winnowmap.sorted.bam"
          ]

if datatype=="pacbio":
    commands3 = [
        "samtools sort -n data.bam -o sorted.bam",
       f"LongGF sorted.bam {GFvoter_path}/gencode.v22.chr_patch_hapl_scaff.annotation.gtf 100 50 100 0 0 1 > LongGF.run.on.XXX.date.log",
        "grep 'SumGF' LongGF.run.on.XXX.date.log > longgf_results.txt"
        ]
if datatype=="nanopore":
    commands3 = [
         "samtools sort -n data.bam -o sorted.bam",
        f"LongGF sorted.bam {GFvoter_path}/gencode.v22.chr_patch_hapl_scaff.annotation.gtf 100 50 100 0 0 2 > LongGF.run.on.XXX.date.log",
        "grep 'SumGF' LongGF.run.on.XXX.date.log > longgf_results.txt"
        ]

commands4=[
        f"{GFvoter_path}/JAFFA-version-2.3/tools/bin/bpipe run {GFvoter_path}/JAFFA-version-2.3/JAFFAL.groovy {input_fastq_gz}"
        ]


print("------------------------- Run minimap and  winnowmap -----------------------------")
run_commands_parallel(commands1, commands2)
print("---------------------- Searching for candidate fusions ----------------------------")
run_commands_parallel(commands3, commands4)


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

LongGF_results=nonrepeart(LongGF_gene,LongGF_fusions)
JAFFAL_results=nonrepeart(JAFFAL_gene,JAFFAL_fusions)
