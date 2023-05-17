import re
import os
import sqlite3
import pandas as pd
import sys

"""
提供核酸|蛋白序列的Accession,查询其对应的分类信息('Accession', 'Taxid', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Virus_name', 'Lineage')
"""

# 参数
RdRp_Accession_file = sys.argv[1]
out_txt = sys.argv[2]

# 定义查询结果的函数
def select_data(table, key, name_list):
	seq = ','.join(['?'] * len(name_list))
	sql = f'select * from {table} where {key} in ({seq})'
	cur.execute(sql, name_list)
	return cur.fetchall()

# 查询结果并处理
def accession2taxonomy(accession):
	found_info = select_data('accession2taxonomy', 'accession', [accession])
	results = found_info[0] if found_info else None
	return results

#virus_dbfile = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/Multi-segmentsVirusFinder/data/RdRp/virus.prot.accession2taxonomy.db"
virus_dbfile = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/Multi-segmentsVirusFinder/data/RdRp/virus.nucl.accession2taxonomy.db"

# 基于RdRp的accession获取分类信息
## 创建数据库和连接
conn = sqlite3.connect(virus_dbfile)
cur = conn.cursor()
## 查询
out_list = []
with open(RdRp_Accession_file, 'r+') as file:
	for line in file:
		accession = line.strip()
		found = accession2taxonomy(accession)
		if found:
			out_list.append(found)

header = ['Accession', 'Taxid', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Virus_name', 'Lineage']
out_df = pd.DataFrame(data = out_list, columns = header)
out_df.to_csv(out_txt, sep='\t', index=False)
## 提交更改
conn.commit()
## 关闭游标和连接
cur.close()
conn.close()
