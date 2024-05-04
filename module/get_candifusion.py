import argparse
import subprocess
import os
import itertools
import re

def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_fastq_gz', type=str, required=False, default="None_file",
                        help='input_path')
    parser.add_argument('-t', action='store', dest='datatype',
                         type=str, required=True, help='datatype')
    parser.add_argument('-s', action='store', dest='choice', type=str, required=False, default="real",
                        help='the source of input read.') 
    parser.add_argument('-o', action='store', dest='output_file', type=str, required=False, default="./GFvoter_out/",
                        help='output directory.')
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


def get_gene(read_chromosome, read_start_pos, read_end_pos, gene_info):
    gene_list = gene_info.get(read_chromosome, [])
    genex=''
    geney=''
    for gene in gene_list:
        if read_start_pos >= int(gene[0]) and read_end_pos <= int(gene[1]):
           genex=gene[2]
           geney=gene[3]
    return genex,geney


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


def fusion_record(align_file,fusion_records):
    from  collections import defaultdict
    with open(align_file, 'r') as file,open(gene_file, 'r') as gene_data:
         lines = file.readlines()
         gene_info = eval(gene_data.read())
         data_dict = defaultdict(list)
         first_elements = [line.split()[0] for line in lines[1:]]  # 获取每行的第一个元素
         for line, key in zip(lines[1:], first_elements):
             elements = line.split()
             read_chromosome=elements[1]
             read_start_pos=int(elements[3])
             read_end_pos=int(elements[4])
             left_soft=int(elements[5])
             right_soft=int(elements[6])
             rl=int(elements[10])
             a=elements[8]
             b=elements[9]
             if a=='secondary':
                continue
             elif first_elements.count(key) > 1:  # 如果第一个元素出现的次数大于1
                genex,geney=get_gene(read_chromosome, read_start_pos, read_end_pos, gene_info)
                if genex !='':
                   data_dict[key].append([genex,read_chromosome, read_start_pos, read_end_pos,geney,left_soft,right_soft,rl,a,b])  # 则将该行其他元素组成的列表添加到对应的值中
    specified_strings1 = ['primary', 'not_supplementary']  # 第一组指定的两个字符串
    specified_strings2 = ['primary', 'supplementary']  # 第二组指定的两个字符串
    new_dict = {}
    for key, value in data_dict.items():
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
# 打印字典
    with open(fusion_records, 'w') as output_file:
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


args = initialization_parameters()
input_fastq_gz=args.input_fastq_gz
repetitive_k15="../../repetitive_k15.txt"

commands1=[
         f"winnowmap -W {repetitive_k15} -ax splice ../../hg38.fa {input_fastq_gz} > winnowmap.sam",
          "samtools view -bS winnowmap.sam > winnowmap.bam",
          "samtools sort winnowmap.bam -o winnowmap.sorted.bam",
          "samtools index winnowmap.sorted.bam"
          ]

commands2=["../../get_alignment_info data.bam alignment.info1"]
commands3=["../../get_alignment_info winnowmap.sorted.bam alignment.info2"]
run_commands(commands1)
run_commands(commands2)
run_commands(commands3)

gene_file = "../../sum_gene_info_sort_nonrepeat.txt"
align_file1="alignment.info1"
align_file2="alignment.info2"
fusion_records="candi_fusiongene_records.txt"
fusion_records1="fusion_records1"
fusion_records2="fusion_records2"
align_file="alignment.info"
fusion_record(align_file1,fusion_records1)
fusion_record(align_file2,fusion_records2)
merged_file(align_file1,align_file2,align_file)
merged_file(fusion_records1,fusion_records2,fusion_records)

total_candifusions= set()
total_records=[]
with open(fusion_records, 'r') as f_in:
     for line in f_in:
         total_records.append(line.strip())
         elements = line.strip().split()
         if len(elements) >= 10 and elements[1] != elements[9]:
            new_element = elements[1] + ":" + elements[9]
            if new_element not in total_candifusions:
               total_candifusions.add(new_element)
