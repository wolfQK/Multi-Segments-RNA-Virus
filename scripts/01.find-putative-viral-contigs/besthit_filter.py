### 获取最优比对
import re
import sys
import argparse
import pandas as pd

# 传参
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--blast_out', type=str, default=None) # 输入blast结果; 比对NR或比对NT均可
parser.add_argument('--best_hit', type=str, default=None) # 输出besthit
parser.add_argument('--filtered_id', type=str, default=None) # 输出最优比对不是病毒的id
args = parser.parse_args()

# 定义基于最优比对过滤掉潜在的非病毒序列
def find_best_hit(blast_out=None, best_hit=None, filtered_id=None):
    if blast_out and best_hit and filtered_id:
        contig_list = []
        besthit_out = []
        filtered_id_list = []
        with open(blast_out, 'r') as blast_info:
            for line in blast_info:
                infos = line.split('\t')
                contig = infos[0]
                #第一个出现的就是最优比对
                if contig not in contig_list:
                    contig_list.append(contig)
                    besthit_info = "\t".join(infos)
                    besthit_out.append(besthit_info)
                    #判断最优比对是否比对上病毒;若没有,则需要被过滤掉
                    if (not re.search('Viruses', besthit_info)):
                        filtered_id_list.append(contig)
                    
        with open(best_hit, 'w') as outfile:
            outfile.writelines(besthit_out)
        
        filterId = pd.DataFrame(filtered_id_list)
        filterId.to_csv(filtered_id, sep='\t', index=None, header=None)
    else:
        print("please input file name of blast_out & best_hit & filtered_id !")

find_best_hit(blast_out=args.blast_out, best_hit=args.best_hit, filtered_id=args.filtered_id)