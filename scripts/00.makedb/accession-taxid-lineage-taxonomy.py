import re
import os
import sqlite3
from itertools import islice
import pandas as pd

"""
构建病毒的核酸和蛋白  Accession-lineage-Taxonomy 的sqlite数据库
"""

#accession2taxid_file = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/database/NCBI-BLAST/2021-05/Taxonomy/accession2taxid/nucl_gb.accession2taxid"
accession2taxid_file = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/Multi-segmentsVirusFinder/data/segmented_virus/metainfo/NCBI_Virus.segments.Accession2Taxid.txt"
accession2taxid_dbfile = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/Multi-segmentsVirusFinder/data/RdRp/nucl.accession2taxid.db"
#accession2taxid_file = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/database/NCBI-BLAST/2021-05/Taxonomy/accession2taxid/prot.accession2taxid.FULL"
#accession2taxid_dbfile = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/Multi-segmentsVirusFinder/data/RdRp/prot.accession2taxid.db"

taxid2lineage_file = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/database/NCBI-BLAST/2021-05/Taxonomy/new_taxdump/fullnamelineage.dmp"
#taxid2lineage_dbfile = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/segmented_virus_lib/RdRp_lib/raw/taxid2lineage.db"
taxid2lineage_dbfile = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/Multi-segmentsVirusFinder/data/RdRp/virus.taxid2lineage.db"
#taxid2taxonomy_dbfile = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/segmented_virus_lib/RdRp_lib/raw/taxid2taxonomy.db"
taxid2taxonomy_dbfile = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/Multi-segmentsVirusFinder/data/RdRp/virus.taxid2taxonomy.db"

virus_taxid2lineage_file = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/database/NCBI-BLAST/2021-05/Taxonomy/new_taxdump/virus.lineage.dmp"
#get_virus_lineage = f'grep -i virus {taxid2lineage_file} > {virus_taxid2lineage_file}'
#os.system(get_virus_lineage)

virus_taxid2lineage_dbfile = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/Multi-segmentsVirusFinder/data/RdRp/virus.taxid2lineage.db"
virus_taxid2taxonomy_dbfile = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/Multi-segmentsVirusFinder/data/RdRp/virus.taxid2taxonomy.db"

virus_dbfile = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/Multi-segmentsVirusFinder/data/RdRp/virus.nucl.accession2taxonomy.db"
virus_dbfile = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/software/Anaconda3-2021.05/envs/Multi-segmentsVirusFinder/data/RdRp/virus.prot.accession2taxonomy.db"

# 创建第一个表accession-taxid并保存到数据库
## 创建连接和游标
conn = sqlite3.connect(accession2taxid_dbfile)
cur = conn.cursor()

tab_accession2taxid = '''
	create table if not exists
	accession2taxid(
	accession text,
	taxid number);
	'''
cur.execute(tab_accession2taxid)

with open(accession2taxid_file, 'r+') as file:
	for line in islice(file, 1, None):
		accession, taxid = line.split()	
		cur.execute('insert into accession2taxid values(?,?)', [accession, taxid])
# 提交更改
conn.commit()
# 关闭游标和连接
cur.close()
conn.close()

"""
# 创建第二个表taxid2lineage并保存到数据库
## 创建连接和游标
#conn = sqlite3.connect(taxid2lineage_dbfile)
conn = sqlite3.connect(virus_taxid2lineage_dbfile)
cur = conn.cursor()

tab_taxid2lineage = '''
	create table if not exists
	taxid2lineage(
	taxid number,
	species text,
	lineage text,
	others text
	);
	'''
cur.execute(tab_taxid2lineage)

lineage_exp = r'(?P<taxid>\d+)\s+\|\s+(?P<virus_name>\w+.*)\s+\|\s+(?P<lineage>\s+|\w+.*)\s+\|(?P<others>.*)'
compile_lineage_exp = re.compile(lineage_exp)
#with open(taxid2lineage_file, 'r+') as file:
with open(virus_taxid2lineage_file, 'r+') as file:
	for line in file:
		match_info = compile_lineage_exp.search(line.strip())
		if match_info:
			taxid = match_info.group('taxid')
			virus_name = match_info.group('virus_name') 
			lineage = match_info.group('lineage') 
			others = match_info.group('others')
			#try:
			cur.execute('insert into taxid2lineage values(?,?,?,?)', [taxid, virus_name, lineage, others])
			#except Exception:
			#	pass
# 提交更改
conn.commit()
# 关闭游标和连接
cur.close()
conn.close()

# 创建第三个表taxid2taxonomy
## 创建连接和游标
conn = sqlite3.connect(virus_taxid2taxonomy_dbfile)
cur = conn.cursor()

tab_taxid2taxonomy = '''
	create table if not exists
	taxid2taxonomy(
	taxid number,
	Kingdom text,
	subKingdom text,
	Phylum text,
	subPhylum text,
	Class text,
	subClass text,
	`Order` text,
	subOrder text,
	Family text,
	subFamily text,
	Genus text,
	Species text,
	virus_name text,
	lineage text
	);
	'''
cur.execute(tab_taxid2taxonomy)

## 设置正则表达式提取分类信息
Kingdom_exp = r'(?P<Kingdom>\w+virae);'
compile_Kingdom_exp = re.compile(Kingdom_exp)
subKingdom_exp = r'(?P<subKingdom>\w+virites);'
compile_subKingdom_exp = re.compile(subKingdom_exp)
Phylum_exp = r'(?P<Phylum>\w+viricota);'
compile_Phylum_exp = re.compile(Phylum_exp)
subPhylum_exp = r'(?P<subPhylum>\w+tina);'
compile_subPhylum_exp = re.compile(subPhylum_exp)
Class_exp = r'(?P<Class>\w+viricetes);'
compile_Class_exp = re.compile(Class_exp)
subClass_exp = r'(?P<subClass>\w+idae);'
compile_subClass_exp = re.compile(subClass_exp)
Order_exp = r'(?P<Order>\w+virales);'
compile_Order_exp =  re.compile(Order_exp)
subOrder_exp = r'(?P<subOrder>\w+virineae);'
compile_subOrder_exp = re.compile(subOrder_exp)
Family_exp = r'(?P<Family>\w+viridae);'
compile_Family_exp = re.compile(Family_exp)
subFamily_exp = r'(?P<subFamily>\w+inae);'
compile_subFamily_exp = re.compile(subFamily_exp)
Genus_exp = r'(?P<Genus>\w+virus);'
compile_Genus_exp = re.compile(Genus_exp)
Species_exp = r'virus; (?P<Species>\w+.*);\t\|'
compile_Species_exp = re.compile(Species_exp)

def get_tax(compile_exp, obj, tax_chr):
	searched = compile_exp.search(obj)
	if searched:
		tax_obj = searched.group(tax_chr)
	else:
		tax_obj = None
	return tax_obj

#with open(taxid2lineage_file, 'r+') as file:
with open(virus_taxid2lineage_file, 'r+') as file:
	for line in file:
		match_info = compile_lineage_exp.search(line.strip())
		if match_info:
			taxid = match_info.group('taxid')
			virus_name = match_info.group('virus_name')
			lineage = match_info.group('lineage')
			Kingdom = get_tax(compile_Kingdom_exp, lineage, 'Kingdom')
			subKingdom = get_tax(compile_subKingdom_exp, lineage, 'subKingdom')
			Phylum = get_tax(compile_Phylum_exp, lineage, 'Phylum')
			subPhylum = get_tax(compile_subPhylum_exp, lineage, 'subPhylum')
			Class = get_tax(compile_Class_exp, lineage, 'Class')
			subClass = get_tax(compile_subClass_exp, lineage, 'subClass')
			Order = get_tax(compile_Order_exp, lineage, 'Order')
			subOrder = get_tax(compile_subOrder_exp, lineage, 'subOrder')
			Family = get_tax(compile_Family_exp, lineage, 'Family')
			subFamily = get_tax(compile_subFamily_exp, lineage, 'subFamily')
			Genus = get_tax(compile_Genus_exp, lineage, 'Genus')
			Species = get_tax(compile_Species_exp, lineage, 'Species')
			#try:
			cur.execute('insert into taxid2taxonomy values(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)', [taxid, Kingdom, subKingdom, Phylum, subPhylum, Class, subClass, Order, subOrder, Family, subFamily, Genus, Species, virus_name, lineage])
			#except Exception:
				#print(lineage)
# 提交更改
conn.commit()
# 关闭游标和连接
cur.close()
conn.close()

"""
# 合并病毒的accession2taxid和taxid2lineage以及taxid2taxonomy
'''virus相关的taxid2taxonomy和accession2taxid'''
## 创建连接和游标
conn = sqlite3.connect(virus_dbfile)
cur = conn.cursor()

attach_database = f'''
	ATTACH DATABASE \'{accession2taxid_dbfile}\' AS 'accession2taxid';
'''
cur.execute(attach_database)
attach_database = f'''
	ATTACH DATABASE \'{virus_taxid2taxonomy_dbfile}\' AS 'taxid2taxonomy';
'''
cur.execute(attach_database)

join='''
	CREATE TABLE accession2taxonomy As 
	SELECT accession, taxid, Kingdom, Phylum, Class, `Order`, Family, Genus, Species, virus_name, lineage 
	FROM taxid2taxonomy.taxid2taxonomy LEFT OUTER JOIN accession2taxid.accession2taxid USING(taxid);
'''		
cur.execute(join)
# 提交更改
conn.commit()
# 关闭游标和连接
cur.close()
conn.close()

"""
# 仅当使用RdRp的时候用;构建的核酸的Accession2Taxonomy是有其他用处

# 定义查询结果的函数
def select_data(table, key, name_list):
	seq = ','.join(['?'] * len(name_list))
	sql = f'select * from {table} where {key} in ({seq})'
	cur.execute(sql, name_list)
	return cur.fetchall()

# 查询结果并处理
def accession2taxonomy(accession):
	found_info = select_data('accession2taxonomy', 'accession', [accession])[0]
	return found_info

# 基于RdRp的accession获取分类信息
## 创建数据库和连接
conn = sqlite3.connect(virus_dbfile)
cur = conn.cursor()
## 查询
out_list = []
RdRp_Accession_file = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/database/RdRp/qk_RdRp/RdRp.SM-blastp2NR-CDD-clusteredCoronaviridae.id"
with open(RdRp_Accession_file, 'r+') as file:
	for line in file:
		accession = line.strip()
		found = accession2taxonomy(accession)
		out_list.append(found)
header = ['Accession', 'Taxid', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Virus_name', 'Lineage']
out_df = pd.DataFrame(data = out_list, columns = header)
out_csv = "/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/database/RdRp/qk_RdRp/RdRp.SM-blastp2NR-CDD-clusteredCoronaviridae.accession2taxonomy.txt"
out_df.to_csv(out_csv, sep='\t', index=False)
## 提交更改
conn.commit()
## 关闭游标和连接
cur.close()
conn.close()

"""
