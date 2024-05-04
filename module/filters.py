import numpy as np
import argparse
import subprocess
import os
import get_candifusion

def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_fastq_gz', type=str, required=False, default="None_file",
                        help='input_path')
    parser.add_argument('-t', action='store', dest='datatype',
                         type=str, required=True, help='datatype')
    parser.add_argument('-n', action='store', dest='min_sup_read', type=int, required=False, default=3,
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

def get_info(record):
    gapsize=0
    leftread_soft_left=int(record.split()[6])
    leftread_soft_right=int(record.split()[7])
    rightread_soft_left=int(record.split()[14])
    rightread_soft_right=int(record.split()[15])
    rl1=int(record.split()[8])
    rl2=int(record.split()[16])
    gapsize=abs(rightread_soft_left-leftread_soft_left-rl1)
    return gapsize

def merged_ovlp_read(detected_fusions,total_records):
    ovlp_read_candifus=set()
    good_records=[]
    for element in detected_fusions:
        fusion_gene = element.strip("'")
        records = get_candifusion.search_fusion_gene(fusion_gene,total_records)
        if int(len(records))>=min_sup_read:
           for record in records:
               gapsize=get_info(record)
               if int(gapsize) <= min_ovlp_len:
                  ovlp_read_candifus.add(fusion_gene)
                  good_records.append(record)
    return ovlp_read_candifus,good_records

def count_secondary(records,alignment_data):
    record_list=[]
    for record in records:
        read_name = record.split()[0]
        i=0
        for line in alignment_data:
            if line.split()[0] == read_name and line.split()[8]=='secondary':
               i+=1
        if i==secondary_alignment:
           record_list.append(record)
    return record_list 
           
   
def calculate_exon_overlap(read_chromosome,genename,read_start,read_end,gene_info):
    gene_list = gene_info.get(read_chromosome, [])
    for gene in gene_list:
        if gene[2].strip('"')==genename:
           exons=gene[4]
           overlap_length = 0
           for exon in exons:
               exon_start = int(exon[1])
               exon_end = int(exon[2])
               if read_start <= exon_end and read_end >= exon_start:
                  overlap_start = max(read_start, exon_start)
                  overlap_end = min(read_end, exon_end)
                  overlap_length += overlap_end - overlap_start + 1
           break
    return overlap_length

def choose_exon_overlap(ovlp_read_candifus,good_records1,alignment_data,gene_info):
     sec_exonovlp_candifus=set()
     good_records2=[]
     for element in ovlp_read_candifus:
         fusion_gene = element.strip()
         records = get_candifusion.search_fusion_gene(fusion_gene,good_records1)
         record_list = count_secondary(records,alignment_data)
         if len(record_list) > 0:
            for record in record_list:
                read_chromosome1 = record.split()[2]
                read_chromosome2 = record.split()[10]
                genename1 = record.split()[1]
                genename2 = record.split()[9]
                read_start1 = int(record.split()[3])
                read_start2 = int(record.split()[11])
                read_end1 = int(record.split()[4])
                read_end2 = int(record.split()[12])
                exon_ovlplength1 = calculate_exon_overlap(read_chromosome1, genename1, read_start1, read_end1,gene_info)
                exon_ovlplength2 = calculate_exon_overlap(read_chromosome2, genename2, read_start2, read_end2,gene_info)
                if int(exon_ovlplength1) > min_exonovlp_len and int(exon_ovlplength2) > min_exonovlp_len :
                   sec_exonovlp_candifus.add(fusion_gene)
                   good_records2.append(record)

     return sec_exonovlp_candifus,good_records2

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
        return "除数不能为零"

def choose_fourratio(fusiongene,fusiongenes,good_records2,gene_info):
    records=get_candifusion.search_fusion_gene(fusiongene,good_records2)
    gene1=fusiongene.split(':')[0]
    gene2=fusiongene.split(':')[1]
    breakpoint1=[]
    breakpoint2=[]
    for record in records:
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
        a11=int(calculate_exon_overlap(chromname1,gene1,r_start1,r_end1,gene_info))
        b1=int(rl1)+int(soft_left1)
        b11=int(calculate_exon_overlap(chromname1,gene1,int(gene_start),r_end1,gene_info))
        a2=int(rl2)
        b2=int(rl2)+int(soft_right2)
        a22=int(calculate_exon_overlap(chromname2,gene2,r_start2,r_end2,gene_info))
        b22=int(calculate_exon_overlap(chromname2,gene2,r_start2,int(gene_end),gene_info))
        align_leftread = calculate_ratio(a1,b1)
        align_leftgene = calculate_ratio(a11,b11)
        align_rightread =calculate_ratio(a2,b2)
        align_rightgene = calculate_ratio(a22,b22)
        align_readexon1=calculate_ratio(a11,a1)
        align_readexon2=calculate_ratio(a22,a2)
        if align_leftread>=read_proportion and align_rightread>=read_proportion and align_leftgene>=gene_proportion and align_rightgene>=gene_proportion:
           fusiongenes.add(fusiongene)
           break
    return fusiongenes  