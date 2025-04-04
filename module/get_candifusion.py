import argparse
import subprocess
import os
import itertools
import re
from  collections import defaultdict
from collections import Counter
import multiprocessing
import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed


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
    parser.add_argument('-f', action='store', dest='fusion_scoring', type=str, required=False, default="F",
                        help='output the scoring of each fusion method during the voting process (F or T)')           
    parser.add_argument('-poll', action='store', dest='min_poll', type=int, required=False, default=6,
                        help='Indicates the minimum number of votes for each reported fusion.')
    parser.add_argument('-ground_truth', action='store', dest='ground_truth', type=str, required=False, default="a",
                        help='Indicates a path to the custom known_fusions.')                    
    args = parser.parse_args()
    return args


def search_fusion_gene(fusion_gene,records):
    fusion_gene_records = []
    for line in records:
        fields = line.strip().split()
        gene1 = fields[1]
        gene2 = fields[9]
        if gene1 == fusion_gene.split(':')[0] and gene2 == fusion_gene.split(':')[1]:
           fusion_gene_records.append(line.strip())
    return fusion_gene_records


def merged_file(file1,file2,merged_file):
    with open(merged_file,'w') as f:
         with open(file1,'r') as f1:
              content1=f1.readlines()
              for line in content1:
                  f.write(line)
         with open(file2,'r') as f2:
              content2=f2.readlines()
              for line in content2:
                  f.write(line)

def run_commands_parallel(commands1, commands2):
    with ThreadPoolExecutor(max_workers=2) as executor:
        future1 = executor.submit(run_command_list, commands1)
        future2 = executor.submit(run_command_list, commands2)
        future1.result()
        future2.result()

def run_command_list(commands):
    for command in commands:
        run_command(command)


def run_command(command):
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"command '{command}' failed to execute. Errors:{result.stderr}")
    return result

# def read_to_gene(info_list, time):
#     elements = info_list.split()
#     read_chromosome=elements[1]
#     read_start_pos=int(elements[3])
#     read_end_pos=int(elements[4])
#     left_soft=int(elements[5])
#     right_soft=int(elements[6])
#     rl=int(elements[10])
#     a=elements[8]
#     b=elements[9]
#     if a =='secondary':
#         return None
#     elif time > 1:
#          with open(gene_list, 'r') as file:
#               for gene_info in file:
#                   elements = line.strip().split(', ')
#                   genex, geney = '', ''
#                   if gene_info[1]==read_chromosome:
#                      if read_start_pos >= int(gene_info[2]) and read_end_pos <= int(gene[3]):
#                         genex=gene[0]
#                         geney=gene[4]
#                         return [genex,read_chromosome, read_start_pos, read_end_pos, geney, left_soft,right_soft,rl,a,b]
#                         break
#     else:
#         return None


def fusion_record_read_file(align_file):
    with open(align_file, 'r') as file:
        lines = file.readlines()
    first_elements = [line.split()[0] for line in lines[1:]]
    count = Counter(first_elements)
    count_list = [count[element] for element in first_elements]
    return first_elements,lines[1:],count_list

def fusion_record(fusionrecord,datadict):
    specified_strings1 = ['primary', 'not_supplementary']  
    specified_strings2 = ['primary', 'supplementary']  
    new_dict = {}
    for key, value in datadict.items():
        specified_strings1_values = [item for item in value if item[-2:] == specified_strings1]
        specified_strings2_values = [item for item in value if item[-2:] == specified_strings2]
        if len(value) == 2:
           if len(specified_strings1_values)==len(specified_strings2_values)==1:
              new_dict[key] = value
        elif len(value) > 2:
           paired_values = list(itertools.product(specified_strings1_values, specified_strings2_values))
           for i, pair in enumerate(paired_values):
               new_key = f"{key}_{i+1}"
               new_dict[new_key] = list(pair)
    with open(fusionrecord, 'w') as output_file:
         for key, values in new_dict.items():
             new_values = [item[:-2] for item in values]
             line = f"{key} {' '.join(map(str, new_values))}\n"
             elements = re.findall(r'\[[^\]]*\]|\S+', line)
             if len(elements) == 3:
                list1 = eval(elements[1])
                list2 = eval(elements[2])
                if list1[0] != list2[0]:
                   if int(list1[5]) < int(list2[5]) and int(list1[6]) > int(list2[6]):
                      list1 = list(map(str, list1))
                      list2 = list(map(str, list2))
                      new_line = elements[0] + ' ' + ' '.join(list1) + ' ' + ' '.join(list2) + '\n'
                      output_file.write(new_line)
                   elif  int(list1[5]) > int(list2[5]) and int(list1[6]) < int(list2[6]):
                         list1 = list(map(str, list1))
                         list2 = list(map(str, list2)) 
                         new_line = elements[0] + ' ' + ' '.join(list2) + ' ' + ' '.join(list1) + '\n'
                         output_file.write(new_line)


def get_result():
    total_candifusions=[]
    total_records=[]
    seen_candifusions = set()
    with open(fusion_records, 'r') as f_in:
        for line in f_in:
            elements = line.strip().split()
            if len(elements) >= 10 and elements[1] != elements[9]:
                new_element = elements[1] + ":" + elements[9]
                total_records.append(line.strip())
                if new_element not in seen_candifusions:
                   seen_candifusions.add(new_element)
                   total_candifusions.append(new_element)

    minimap2_results=set()
    winnowmap_results=set()
    with open(fusion_records1, 'r') as f:
        for line in f:
            elements = line.strip().split()
            if len(elements) >= 10 and elements[1] != elements[9]:
                new_element = elements[1] + ":" + elements[9]
                if new_element not in minimap2_results:
                    minimap2_results.add(new_element)
    with open(fusion_records2, 'r') as f:
        for line in f:
            elements = line.strip().split()
            if len(elements) >= 10 and elements[1] != elements[9]:
                new_element = elements[1] + ":" + elements[9]
                if new_element not in winnowmap_results:
                    winnowmap_results.add(new_element)
    return total_candifusions,total_records,minimap2_results,winnowmap_results


GFvoter_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
GFvoter_path = os.path.normpath(GFvoter_path)
args = initialization_parameters()
input_fastq_gz=args.input_fastq_gz
datatype=args.datatype
output_file_path=args.output_file
align_file="alignment.info"
fusion_records="candi_fusiongene_records.txt"
fusion_records1="fusion_records1"
fusion_records2="fusion_records2"


if __name__ == "__main__":

    run_command(f"chmod +x {GFvoter_path}/get_alignment_info")
    run_command(f"chmod +x {GFvoter_path}/read2gene")
    run_command(f"chmod +x {GFvoter_path}/mergeAndsort")

    align_file1="alignment.info1"
    align_file2="alignment.info2"
    minimap_result="data.bam"
    winnowmap_result="winnowmap.sorted.bam"
    
    run_commands_parallel([f"{GFvoter_path}/get_alignment_info {minimap_result} {align_file1}"], [f"{GFvoter_path}/get_alignment_info {winnowmap_result} {align_file2}"])
    run_commands_parallel([f"{GFvoter_path}/read2gene {GFvoter_path}/doc/sum_gene_info_sort_nonrepeat.txt1 {align_file1} {fusion_records1}"], [f"{GFvoter_path}/read2gene {GFvoter_path}/doc/sum_gene_info_sort_nonrepeat.txt1 {align_file2} {fusion_records2}"])
    merged_file(align_file1,align_file2,align_file)
    run_command(f"{GFvoter_path}/mergeAndsort fusion_records1 fusion_records2 candi_fusiongene_records.txt")
