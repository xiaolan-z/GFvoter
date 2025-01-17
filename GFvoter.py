import datetime
import argparse
import os
import sys
from pathlib import Path
current_dir = Path(__file__).resolve().parent
module_dir = current_dir / 'module'
sys.path.append(str(module_dir))

import longgf_jaffal
import get_candifusion
import Scoring
import time
import json



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


def union(set1,set2):
    union_fusions=set(set1)
    for line in set2:
        element = line.strip()
        element1, element2 = line.strip().split(':')
        new = element2 + ':' + element1
        if element not in set1 and new not in set1:
            union_fusions.add(element)
    return union_fusions

    
def voter(union_set):
    voters_dict={}
    voter_fusions=set()
    for fusion in union_set:
        voters=0
        fusion1, fusion2 = fusion.split(':')
        fusion_reversed = fusion2 + ":" + fusion1
        if fusion in minimap2_results or fusion_reversed in minimap2_results:
           voters+=1
        if fusion in winnowmap_results or fusion_reversed in winnowmap_results:
           voters+=1
        if fusion in LongGF_gene or fusion_reversed in LongGF_gene:
           voters+=3
        if fusion in JAFFAL_gene or fusion_reversed in JAFFAL_gene:
           voters+=3
        if fusion in Scoring_list or fusion_reversed in Scoring_list:
           voters+=6
        voters_dict.update({fusion:voters})
        if voters>=min_poll:
           voter_fusions.add(fusion)
    sorted_dict = dict(sorted(voters_dict.items(), key=lambda x: x[1], reverse=True))
    with open("./voter_list",'w') as f:
         for key, value in sorted_dict.items():
             f.write(f"{key}\t{value}\n")
    return voter_fusions


def scorforfusion():
   with open(align_file, 'r') as f,open(gene_file, 'r') as gene_data:
      alignment_dict={}
      for line in f:
         parts = line.split()
         if parts[0] not in alignment_dict:
            alignment_dict[parts[0]] =[]
         alignment_dict[parts[0]].append(parts[8])
      gene_info = json.load(gene_data)
   Scoring_dict={}
   Scoring_list=[]
   count=0
   current_line_index=0
   for fusion in total_candifusions:
      gene1=fusion.split(":")[0]
      gene2=fusion.split(":")[1]
      records,reads_num,current_line_index=Scoring.get_reads_num(gene1,gene2,total_records,current_line_index)
      suma=0
      i=0
      for record in records:
         sec=Scoring.get_secondary(record,alignment_dict)
         distance=Scoring.get_distance(record)
         exon_ovlp1,exon_ovlp2,prl,pgl,prr,pgr=Scoring.get_twoexonovlp_and_fourratio(gene1,gene2,gene_info,record)
         s=Scoring.get_score(reads_num,distance,sec,exon_ovlp1,exon_ovlp2,prl,pgl,prr,pgr)
         if s>0:
            i+=1
            suma+=s
      if i>0:
         score=round(suma/i) 
      else:
         score=0
      Scoring_dict.update({fusion:score})
      if score>min_score:
         Scoring_list.append(fusion)
      count += 1
    
      if count % 5000 == 0:
         current_time = datetime.datetime.now()
   sorted_dict = dict(sorted(Scoring_dict.items(), key=lambda x: x[1], reverse=True))
   with open("./Scoring_list",'w') as f:
      for key, value in sorted_dict.items():
         f.write(f"{key}\t{value}\n")
   return Scoring_list



def main():
    voter_fusions=voter(union_set)
    GFvoter_results=[]
    for fusion in voter_fusions:
        fusion1, fusion2 = fusion.split(':')
        fusion_reversed = fusion2 + ":" + fusion1
        read_num=0
        breakpoint_sum1=[]
        breakpoint_sum2=[]
        for item in LongGF_results:
            if fusion==item[0]:
               if int(read_num)<int(item[1]):
                  read_num=item[1]
               breakpoint_sum1.append(int(item[3]))
               breakpoint_sum2.append(int(item[5]))
               chrom1=item[2]
               chrom2=item[4]
            elif fusion_reversed==item[0]:
                 if int(read_num)<int(item[1]):
                    read_num=item[1]
                 breakpoint_sum1.append(int(item[5]))
                 breakpoint_sum2.append(int(item[3]))
                 chrom1=item[2]
                 chrom2=item[4]
            else:
                continue
        for item in JAFFAL_results:
            if fusion==item[0]:
               if int(read_num)<int(item[1]):
                  read_num=item[1]
               breakpoint_sum1.append(int(item[3]))
               breakpoint_sum2.append(int(item[5]))
               chrom1=item[2]
               chrom2=item[4]
            elif fusion_reversed==item[0]:
                 if int(read_num)<int(item[1]):
                    read_num=item[1]
                 breakpoint_sum1.append(int(item[5]))
                 breakpoint_sum2.append(int(item[3]))
                 chrom1=item[2]
                 chrom2=item[4]
            else:
                continue
        for item in Scoring_results:
            if fusion==item[0]:
               if int(read_num)<int(item[1]):
                  read_num=item[1]
               breakpoint_sum1.append(int(item[3]))
               breakpoint_sum2.append(int(item[5]))
               chrom1=item[2]
               chrom2=item[4]
            elif fusion_reversed==item[0]:
                 if int(read_num)<int(item[1]):
                    read_num=item[1]
                 breakpoint_sum1.append(int(item[5]))
                 breakpoint_sum2.append(int(item[3]))
                 chrom1=item[2]
                 chrom2=item[4]
            else:
                continue
        breakpoint1=sum(breakpoint_sum1)/len(breakpoint_sum1)
        breakpoint2=sum(breakpoint_sum2)/len(breakpoint_sum2)
        fusion_info=[fusion,read_num,chrom1,breakpoint1,chrom2,breakpoint2]
        GFvoter_results.append(fusion_info)  

    known_fusions=[]
    if ground_truth=="a":
      if choice == 'real':
          known_fusion = f'{GFvoter_path}/doc/known_fusions.txt'
      elif choice == 'simulate':
          known_fusion = f'{GFvoter_path}/doc/smi_known_fusions.txt'
    else:
         known_fusion=ground_truth
    with open(known_fusion, 'r') as f:
        for line in f:
            genes = line.split()  
            fusion = ":".join(genes)  
            known_fusions.append(fusion)
    GFvoter_known_fusion=[]
    for fusion in voter_fusions:
        fusion1, fusion2 = fusion.split(':') 
        fusion_reversed = fusion2 + ":" + fusion1  
        if fusion in known_fusions or fusion_reversed in known_fusions:
           GFvoter_known_fusion.append(fusion)
    with open("./GFvoter_gene",'w') as f:
         for item in voter_fusions:
             f.write(str(item) + '\n')

      
    sorted_fusions = sorted(GFvoter_results, key=lambda x: int(x[1]), reverse=True)
    count = sum(1 for item in sorted_fusions if float(item[1]) < 2)
    print(len(sorted_fusions))
    if count >= 3*(len(sorted_fusions)-count) and len(sorted_fusions)!=count:
        sorted_fusions = [item for item in sorted_fusions if float(item[1]) >= 2]
    reported_fusion = "../reported_fusions.txt"
    with open(reported_fusion, "w") as file:
         file.write("ID\tGene Fusion\tsupporting_reads\tchrom1\tbreakpoint1\tchrom2\tbreakpoint2\tKnown\n")  
         for idx, item in enumerate(sorted_fusions, start=1):
             gene_fusion = item[0]
             supporting_reads = item[1]
             chrom1 = item[2]
             breakpoint1 = round(item[3])
             chrom2 = item[4]
             breakpoint2 = round(item[5])
             known_status = "yes" if gene_fusion in GFvoter_known_fusion else "" 
             if int(supporting_reads)>=min_sup_read:
                file.write(f"{idx}\t{gene_fusion}\t{supporting_reads}\t{chrom1}\t{breakpoint1}\t{chrom2}\t{breakpoint2}\t{known_status}\n")   


if __name__ == "__main__":
   args = initialization_parameters()
   input_fastq_gz=args.input_fastq_gz
   min_sup_read=args.min_sup_read
   choice=args.choice
   output_file=args.output_file
   datatype=args.datatype
   ground_truth=args.ground_truth
   GFvoter_path = os.path.abspath(os.path.dirname(__file__))
   if output_file[-1]!='/':
      output_file+='/'
   min_score=args.min_score
   min_poll=args.min_poll
   gene_file = f"{GFvoter_path}/doc/sum_gene_info_sort_nonrepeat.txt"
   os.system(f"python {module_dir}/get_candifusion.py -i {input_fastq_gz} -t {datatype} -o {output_file} -s {choice}")
   align_file = get_candifusion.align_file
   LongGF_gene=longgf_jaffal.LongGF_gene
   JAFFAL_gene=longgf_jaffal.JAFFAL_gene
   JAFFAL_results=longgf_jaffal.JAFFAL_results
   LongGF_results=longgf_jaffal.LongGF_results
 

   total_candifusions,total_records,minimap2_results,winnowmap_results=get_candifusion.get_result()
   union_set=union(union(total_candifusions,LongGF_gene),JAFFAL_gene)
   print("--------------------------- Scoring for candidate fusions ---------------------------------")
   Scoring_list=scorforfusion()
   Scoring_results=[]
   for fusiongene in Scoring_list:
      records,read_num=Scoring.get_reads_num1(fusiongene,total_records)
      gene1=fusiongene.split(":")[0]
      gene2=fusiongene.split(":")[1]
      breakpoint1_values=[]
      breakpoint2_values=[]
      for record in records:
         elements=record.split()
         chromosome1=elements[2]
         chromosome2=elements[10]
         breakpoint1_values.append(int(elements[4]))
         breakpoint2_values.append(int(elements[11]))
      breakpoint1=sum(breakpoint1_values)/int(read_num)
      breakpoint2=sum(breakpoint2_values)/int(read_num)
      fusion=[gene1+":"+gene2,read_num,chromosome1,round(breakpoint1),chromosome2,round(breakpoint2)]
      Scoring_results.append(fusion)
   main()
   end_time = datetime.datetime.now()
   print(f"---------------------- Pineline finished at {end_time} ----------------------------")

   


