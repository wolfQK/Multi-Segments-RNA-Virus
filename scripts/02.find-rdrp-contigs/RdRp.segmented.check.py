
import re
import sqlite3
import argparse
import pandas as pd

"""
或许这个也应该被写成一个模块？
基于比对上的RdRp的分类信息,确定潜在的病毒序列是节段病毒还是非节段病毒
潜在病毒序列比对RdRp的结果:qseqid qlen sseqid slen qstart qend sstart send evalue bitscore length pident mismatch gaps
"""

### 传参
# Accession2Taxonomy库
virus_Accession2Taxonomy_db = '/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Multi-Segments-Virus/data/RdRp/db_sqlite3/virus.prot.accession2taxonomy.db'
# 节段病毒分类信息库
segmented_FGS_xlsx = '/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Multi-Segments-Virus/data/multi_segments_virus/share/ICTV_VH.checked.segmented.FGS.xlsx'

parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--blastx2RdRp_out', type=str, default=None)
parser.add_argument('--sample_out', type=str, default=None)
parser.add_argument('--related_p', type=float, default=100)
args = parser.parse_args()

# 比对RdRp的结果
blastx2RdRp_out = args.blastx2RdRp_out
# 输出判断结果
sample_out = args.sample_out
# 相差范围 
related_p = args.related_p 

# 比对RdRp的结果
blastx2RdRp_info = pd.read_csv(blastx2RdRp_out, sep='\t', names=['qseqid', 'qlen', 'sseqid', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'length', 'pident', 'mismatch', 'gaps'])
contig_id_list = list(blastx2RdRp_info['qseqid'].drop_duplicates()) 

# 连接Accession2Taxonomy数据库
conn = sqlite3.connect(virus_Accession2Taxonomy_db)
cur = conn.cursor()

# 使用来自NCBI_Virus的节段病毒分类信息
segmented_Order_list = ['Bunyavirales']
segmented_Family = pd.read_excel(segmented_FGS_xlsx, sheet_name='Family')
segmented_Family_list = list(segmented_Family['Family'])
segmented_Genus = pd.read_excel(segmented_FGS_xlsx, sheet_name='Genus')
segmented_Genus_list = list(segmented_Genus['Genus'])
segmented_Species = pd.read_excel(segmented_FGS_xlsx, sheet_name='Species')
segmented_Species_list = list(segmented_Species['Species'])

only_segmented_Family = list(segmented_Family.loc[segmented_Family['share'] == 'n']['Family'])
share_segmented_Family = list(segmented_Family.loc[segmented_Family['share'] == 'y']['Family'])

Order_exp = r'(?P<Order>\w+virales);'
compile_Order_exp =  re.compile(Order_exp)
Family_exp = r'(?P<Family>\w+viridae);'
compile_Family_exp = re.compile(Family_exp)
Genus_exp = r'(?P<Genus>\w+virus);'
compile_Genus_exp = re.compile(Genus_exp)
Species_exp = r'(?P<Species>\w+.* virus.*)'
compile_Species_exp = re.compile(Species_exp)

# 在get_taxonomy_df中用来查询数据库的函数
def select_data(table, key, name_list):
    seq = ','.join(['?'] * len(name_list))
    sql = f'select * from {table} where {key} in ({seq})'
    cur.execute(sql, name_list)
    return cur.fetchall()

def accession2taxonomy(accession):
    found_info = select_data('accession2taxonomy', 'accession', [accession])
    results = found_info[0] if found_info else None
    return results

# 对于一条contigs比对上的RdRp的分类信息取LCA
def get_LCA(Taxonomy_df=None):
    LCA_list = []
    rtaxes = ['Virus_name', 'Species', 'Genus', 'Family', 'Order',  'Class', 'Phylum', 'Kingdom']
    taxes = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Virus_name']
    for tax in taxes:
        uniq_tax = list(Taxonomy_df.loc[:,tax].drop_duplicates())
        #print(uniq_tax)
        if uniq_tax == [None]:
            #print(tax)
            continue
        elif len(uniq_tax) == 1:
            LCA_list.append(uniq_tax[0])
        else :
            break
    return LCA_list

# 获取与max bitscore相差不超过(1-related_p)的蛋白的分类信息
# 并对分类信息取LCA
# 输入contig_id_list返回这些contigs取LCA之后的科属种，以及是否为节段病毒
def get_taxonomy_LCA_segmented(related_p=related_p, contig_id_list=contig_id_list):
    out_list = []
    for contig_id in contig_id_list:
        tmp_df = blastx2RdRp_info[blastx2RdRp_info['qseqid'] == contig_id]
        threshold = (related_p/100) * tmp_df['bitscore'].max()
        tmp_use = tmp_df[tmp_df['bitscore'] >= threshold]
        RdRp_id_list = list(tmp_use['sseqid'].drop_duplicates())
        RdRp_taxonomy_list = []
        for RdRp_id in RdRp_id_list:
            found = accession2taxonomy(RdRp_id)
            if found:
                RdRp_taxonomy_list.append(found)
        header = ['Accession', 'Taxid', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Virus_name', 'Lineage']
        RdRp_taxonomy_df = pd.DataFrame(data=RdRp_taxonomy_list, columns=header)
        LCA_str = ';'.join(get_LCA(RdRp_taxonomy_df))
    
        try:
            LCA_Order = compile_Order_exp.search(LCA_str).group('Order')
        except:
            LCA_Order = None
        try:
            LCA_Family = compile_Family_exp.search(LCA_str).group('Family')
        except:
            LCA_Family = None
        segmentedFamily = LCA_Family in segmented_Family_list
        try:                                       
            LCA_Genus = compile_Genus_exp.search(LCA_str).group('Genus')
        except:
            LCA_Genus = None
        segmentedGenus = LCA_Genus in segmented_Genus_list
        try:
            tmp = LCA_str.split(';')[-1]
            LCA_Species =  tmp if compile_Species_exp.search(tmp) else None
        except:
            LCA_Species = None
        segmentedSpecies = LCA_Species in segmented_Species_list
        segmentedOrder = LCA_Order in segmented_Order_list or segmentedFamily or segmentedGenus or segmentedSpecies

        # 基于RdRp的lineage来确定其是否为节段病毒的RdRp
        if segmentedGenus: #识别到属且是节段病毒
            segmented_RdRp = 1
        elif LCA_Genus: #识别到属,但不是节段病毒
            segmented_RdRp = 0
        elif not LCA_Genus:
            segmented_RdRp = 3 #没有识别到属
            
        if segmented_RdRp == 3 and LCA_Family in only_segmented_Family: #没有识别到属但识别到科且科下全为节段病毒
            segmented_RdRp = 1
        elif segmented_RdRp == 3 and LCA_Family in share_segmented_Family: #没有识别到属但识别到科且科下包含节段病毒
            segmented_RdRp = 2
        elif segmented_RdRp == 3 and LCA_Family: #没有识别到属但识别到科且科下不包含节段病毒
            segmented_RdRp = 0
        elif segmented_RdRp == 3 and not LCA_Family :
            segmented_RdRp = 4 #没有识别到科
        
        if segmented_RdRp == 4 and LCA_Order == 'Bunyavirales': #没有识别到科和属,但确定是布尼亚病毒
            segmented_RdRp = 1
            
        out_list.append([contig_id, segmented_RdRp, LCA_Order, segmentedOrder, LCA_Family, segmentedFamily, LCA_Genus, segmentedGenus, LCA_Species, segmentedSpecies])
    return pd.DataFrame(data=out_list, columns=['contig_id', 'segmented_RdRp', 'LCA_Order', 'segmentedOrder', 'LCA_Family', 'segmentedFamily', 'LCA_Genus', 'segmentedGenus', 'LCA_Species', 'segmentedSpecies'])

sample_results = get_taxonomy_LCA_segmented()
sample_results.to_csv(sample_out, sep='\t', index=False)

# 关闭数据库和连接
cur.close()
conn.close()

