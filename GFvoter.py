import datetime
start_time = datetime.datetime.now()
print(start_time)
import argparse
import sys
sys.path.append("./module/")
import longgf_jaffal
sys.path.append("../../module/")
import get_candifusion
import Scoring                                                                                                              



def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='store', dest='input_fastq_gz', type=str, required=False, default="None_file",
                        help='input_path')
    parser.add_argument('-t', action='store', dest='datatype',
                         type=str, required=False, default="pacbio", help='datatype')
    parser.add_argument('-s', action='store', dest='choice', type=str, required=False, default="real",
                        help='the source of input read.') 
    parser.add_argument('-o', action='store', dest='output_file', type=str, required=False, default="./GFvoter_out/",
                        help='output directory.')               
    parser.add_argument('-score', action='store', dest='min_score', type=int, required=False, default=400,
                        help='the minimum score for fusion genes.')          
    parser.add_argument('-poll', action='store', dest='min_poll', type=int, required=False, default=6,
                        help='the minimum number of votes for fusion genes.')                                                         
    args = parser.parse_args()
    return args

args = initialization_parameters()
choice=args.choice
output_file=args.output_file
if output_file[-1]!='/':
    output_file+='/'
min_score=args.min_score
min_poll=args.min_poll
gene_file = "../../sum_gene_info_sort_nonrepeat.txt"
align_file = get_candifusion.align_file
total_candifusions=get_candifusion.total_candifusions
total_records=get_candifusion.total_records
minimap2_results=get_candifusion.minimap2_results
winnowmap_results=get_candifusion.winnowmap_results
LongGF_gene=longgf_jaffal.LongGF_gene
JAFFAL_gene=longgf_jaffal.JAFFAL_gene
JAFFAL_results=longgf_jaffal.JAFFAL_results
LongGF_results=longgf_jaffal.LongGF_results


def union(set1,set2):
    union_fusions=set(set1)
    for line in set2:
        element = line.strip()
        element1, element2 = line.strip().split(':')
        new = element2 + ':' + element1
        if element not in set1 and new not in set1:
            union_fusions.add(element)
    return union_fusions
union_set=union(union(total_candifusions,LongGF_gene),JAFFAL_gene)

    
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
    
    
with open(align_file, 'r') as f,open(gene_file, 'r') as gene_data:
     alignment_data = f.readlines()
     gene_info = {}
     gene_info = eval(gene_data.read())
Scoring_dict={}
Scoring_list=[]
for fusion in total_candifusions:
    gene1=fusion.split(":")[0]
    gene2=fusion.split(":")[1]
    records,reads_num=Scoring.get_reads_num(fusion,total_records)
    suma=0
    i=0
    for record in records:
        sec=Scoring.get_secondary(record,alignment_data)
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
sorted_dict = dict(sorted(Scoring_dict.items(), key=lambda x: x[1], reverse=True))
with open("./Scoring_list",'w') as f:
     for key, value in sorted_dict.items():
        f.write(f"{key}\t{value}\n")



Scoring_results=[]
for fusiongene in Scoring_list:
    records,read_num=Scoring.get_reads_num(fusiongene,total_records)
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
    if choice == 'real':
        known_fusion = '../../known_fusions.txt'
    elif choice == 'simulate':
        known_fusion = '../../smi_known_fusions.txt'
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
    reported_fusion = "../GFvoter_fusions.txt"
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
             file.write(f"{idx}\t{gene_fusion}\t{supporting_reads}\t{chrom1}\t{breakpoint1}\t{chrom2}\t{breakpoint2}\t{known_status}\n")   


if __name__ == "__main__":
   main()
end_time = datetime.datetime.now()
print(end_time)


