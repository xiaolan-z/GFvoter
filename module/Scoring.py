import numpy as np
import argparse
import subprocess
import os
import get_candifusion
import math

def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_fastq_gz', type=str, required=False, default="None_file",
                        help='input_path')
    parser.add_argument('-t', action='store', dest='datatype',
                         type=str, required=True, help='datatype')
    parser.add_argument('-n', action='store', dest='min_sup_read', type=int, required=False, default=2,
                        help='the number of supporting reads of gene fusion.')
    parser.add_argument('-l', action='store', dest='min_ovlp_len', type=int, required=False, default=15,
                        help='the minimum overlap length between two aligments of one read.')
    parser.add_argument('-el', action='store', dest='min_exonovlp_len', type=int, required=False, default=400,
                        help='the minimum exon overlap length of an alignment record against the reference genome.')
    parser.add_argument('-sn', action='store', dest='secondary_alignment', type=int, required=False, default=0,
                        help='the number of secondary alignment of a read.') 
    parser.add_argument('-rp', action='store', dest='read_proportion', type=float, required=False, default=0.75,
                        help='the propotion of one alignment to the read.')
    parser.add_argument('-gp', action='store', dest='gene_proportion', type=float, required=False, default=0.3,
                        help='the propotion of one alignment to the genome.')  
    parser.add_argument('-o', action='store', dest='output_file', type=str, required=False, default="./GFvoter_out/",
                        help='output directory.')
    parser.add_argument('-s', action='store', dest='choice', type=str, required=False, default="real",
                        help='the source of input read.')                                                                                    
    args = parser.parse_args()
    return args

args = initialization_parameters()
min_sup_read=args.min_sup_read
min_ovlp_len=args.min_ovlp_len
min_exonovlp_len=args.min_exonovlp_len
secondary_alignment=args.secondary_alignment
read_proportion=args.read_proportion
gene_proportion=args.gene_proportion

def get_reads_num(fusion_gene,records):
    fusion_gene_records = []
    for line in records:
        fields = line.strip().split()
        gene1 = fields[1]
        gene2 = fields[9]
        if gene1 == fusion_gene.split(':')[0] and gene2 == fusion_gene.split(':')[1]:
           fusion_gene_records.append(line.strip())
    reads_num=len(fusion_gene_records)
    return fusion_gene_records,int(reads_num)

def get_distance(record):
    distance=0
    leftread_soft_left=int(record.split()[6])
    leftread_soft_right=int(record.split()[7])
    rightread_soft_left=int(record.split()[14])
    rightread_soft_right=int(record.split()[15])
    rl1=int(record.split()[8])
    rl2=int(record.split()[16])
    distance=abs(rightread_soft_left-leftread_soft_left-rl1)
    return int(distance)

def get_secondary(record,alignment_data):
    sec=0
    read_name = record.split()[0]
    i=0
    for line in alignment_data:
        if line.split()[0] == read_name and line.split()[8]=='secondary':
           i+=1
    if i==0:
       sec=1
    else:
       sec=0
    return int(sec)


def get_exon_ovlp(read_chromosome,genename,read_start,read_end,gene_info):
    gene_list = gene_info.get(read_chromosome, [])
    for gene in gene_list:
        if gene[2].strip('"')==genename:
           exons=gene[4]
           exon_ovlp = 0
           for exon in exons:
               exon_start = int(exon[1])
               exon_end = int(exon[2])
               if read_start <= exon_end and read_end >= exon_start:
                  overlap_start = max(read_start, exon_start)
                  overlap_end = min(read_end, exon_end)
                  exon_ovlp += overlap_end - overlap_start + 1
           break
    return int(exon_ovlp)
 
def calculate_variance(data):
    variance = round(np.var(data),3)
    return variance

def get_geneinfo(chromname,gene_info,genename):
    gene_list = gene_info.get(chromname, [])
    for gene in gene_list: 
        if gene[2].strip('"')==genename:
           gene_start=gene[0]
           gene_end=gene[1]
           break
    return gene_start,gene_end 

def calculate_ratio(num1, num2):
    if num2 != 0:
        ratio = round(num1 / num2, 3)
        return ratio
    else:
        return 0
def get_twoexonovlp_and_fourratio(gene1,gene2,gene_info,record):
    r_start1=int(record.split()[3])
    r_end1=int(record.split()[4])
    soft_left1=record.split()[6]
    r_start2=int(record.split()[11])
    r_end2=int(record.split()[12])
    soft_right2=record.split()[15]
    chromname1=record.split()[2]
    chromname2=record.split()[10]
    rl1=record.split()[8]
    rl2=record.split()[16]
    gene_start,x=get_geneinfo(chromname1,gene_info,gene1)
    y,gene_end=get_geneinfo(chromname2,gene_info,gene2)
    a1=int(rl1)
    a11=int(get_exon_ovlp(chromname1,gene1,r_start1,r_end1,gene_info))
    b1=int(rl1)+int(soft_left1)
    b11=int(get_exon_ovlp(chromname1,gene1,int(gene_start),r_end1,gene_info))
    a2=int(rl2)
    b2=int(rl2)+int(soft_right2)
    a22=int(get_exon_ovlp(chromname2,gene2,r_start2,r_end2,gene_info))
    b22=int(get_exon_ovlp(chromname2,gene2,r_start2,int(gene_end),gene_info))
    exon_ovlp1=a11
    exon_ovlp2=a22
    prl = calculate_ratio(a1,b1)
    pgl= calculate_ratio(a11,b11)
    prr=calculate_ratio(a2,b2)
    pgr= calculate_ratio(a22,b22)
    return int(exon_ovlp1),int(exon_ovlp2),round(prl,2),round(pgl,2),round(prr,2),round(pgr,2)
 

def get_score(reads_num,distance,sec,exon_ovlp1,exon_ovlp2,prl,pgl,prr,pgr):
    m=round(min(prl,prr)/read_proportion-1,1)
    n=round(min(pgl,pgr)/gene_proportion-1,1)
    if exon_ovlp1>min_exonovlp_len and exon_ovlp2>min_exonovlp_len and distance<=min_ovlp_len:
      s=(reads_num-min_sup_read)*(round(math.log10(exon_ovlp1),1)+round(math.log10(exon_ovlp2),1))*math.pow(2,10*min(m,n))*sec*(16-distance)
    else:
       s=0
    return round(s)


