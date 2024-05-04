import datetime
start_time = datetime.datetime.now()
print(start_time)
import argparse
import sys
sys.path.append("./module/")
import longgf_jaffal
sys.path.append("../../module/")
import get_candifusion
import filters                                                                                                                              



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
    args = parser.parse_args()
    return args

args = initialization_parameters()
choice=args.choice
output_file=args.output_file
if output_file[-1]!='/':
    output_file+='/'
gene_file = "../../sum_gene_info_sort_nonrepeat.txt"
align_file = get_candifusion.align_file
total_candifusions=get_candifusion.total_candifusions
total_records=get_candifusion.total_records
common_fusion=longgf_jaffal.common_fusion



def main():
    # 读取文件并处理数据
    with open(align_file, 'r') as f, open(gene_file, 'r') as gene_data:
        alignment_data = f.readlines()
        gene_info = eval(gene_data.read())
        ovlp_read_candifus, good_records1 = filters.merged_ovlp_read(total_candifusions, total_records)
        sec_exonovlp_candifus, good_records2 = filters.choose_exon_overlap(ovlp_read_candifus, good_records1, alignment_data, gene_info)
        fusiongenes = set()
        for fusiongene in sec_exonovlp_candifus: 
            fusiongenes = filters.choose_fourratio(fusiongene, fusiongenes, good_records2,gene_info)
    
        
    gene_fusions=set(common_fusion)
    for line in fusiongenes:
        element = line.strip()
        element1, element2 = line.strip().split(':')
        new = element2 + ':' + element1
        if element not in gene_fusions and new not in gene_fusions:
            gene_fusions.add(element)

    # 读取已知融合基因数据
    if choice == 'real':
        known_fusion = '../../known_fusions.txt'
    elif choice == 'simulate':
        known_fusion = '../../smi_known_fusions.txt'
    known_fusions = []
    with open(known_fusion, 'r') as f:
        for line in f:
            genes = line.split()  # 使用空格拆分基因字符串
            fusion = ":".join(genes)  # 以冒号连接两个基因
            known_fusions.append(fusion)

    # 生成报告
    

    report_fusion = "../reported_fusions.txt"
    with open(report_fusion, "w") as file:
        file.write("ID\tGene Fusion\tKnown\n")  # 写入表头

        # 遍历融合基因列表，写入文件
        for idx, gene in enumerate(gene_fusions, start=1):
            known_status = "yes" if gene in known_fusions else ""  # 判断是否为已知融合基因
            file.write(f"{idx}\t{gene}\t{known_status}\n")  # 写入每行数据

if __name__ == "__main__":
    main()
end_time = datetime.datetime.now()
print(end_time)
