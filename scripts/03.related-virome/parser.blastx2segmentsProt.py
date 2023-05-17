import pandas as pd
import sqlite3
import argparse
import json

# 传参
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--blastx_out_file', type=str, default=None)
parser.add_argument('--related_p', type=float, default=100)
#parser.add_argument('--sample', type=str, default=None)
parser.add_argument('--RdRp_check_file', type=str, default=None)
parser.add_argument('--out_file', type=str, default=None)
parser.add_argument('--virome_out_anyhit', type=str, default=None)
parser.add_argument('--virome_out_besthit', type=str, default=None)
args = parser.parse_args()
segmented_FGS_xlsx = '/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Multi-Segments-Virus/data/multi_segments_virus/share/ICTV_VH.checked.segmented.FGS.xlsx'
prot_accession2taxonomy_db = '/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Multi-Segments-Virus/data/RdRp/db_sqlite3/virus.prot.accession2taxonomy.db'
## blastx结果
blastx_out_file = args.blastx_out_file
## 相差范围 
related_p = args.related_p
## 样本名称
#sample = args.sample # 主要是因为比对RdRp的contigs没有添加样本标签
## 比对RdRp库的结果
RdRp_check_file = args.RdRp_check_file
## 输出文件
out_file = args.out_file
virome_out_anyhit = args.virome_out_anyhit
virome_out_besthit = args.virome_out_besthit

# 1 信息读取
## 1.1 blastx结果
names = names = ['qseqid', 'qlen', 'sseqid', 'slen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'length', 'pident', 'mismatch', 'gaps']
blastx_out = pd.read_csv(blastx_out_file, sep='\t', names=names)


## 1.2 RdRp比对与判定的结果
# 基于比对RdRp的结果判断是否为节段病毒RdRp的结果
RdRp_check_results = pd.read_csv(RdRp_check_file, sep='\t') 
RdRp_contigs = set(RdRp_check_results['contig_id'])


## 1.3 节段病毒库中的基因组信息
# 使用来自ICTV_VH的节段病毒分类信息
segmented_Genus = pd.read_excel(segmented_FGS_xlsx, sheet_name='Genus')
segmented_Genus_list = list(segmented_Genus['Genus'])
segmented_Species = pd.read_excel(segmented_FGS_xlsx, sheet_name='Species')
segmented_Species_list = list(segmented_Species['Species'])

# 1.4 ICTV_VH的节段病毒基因组信息
#segmented_genome_xlsx = '/Users/fengqikai/JupyterLab/by-test/data/Final/ICTV_VH.checked.segmented_genome.use.xlsx'
#segmented_genome_info = pd.read_excel(segmented_genome_xlsx)

# 1.5 ICTV_VH的节段信息
# segments_xlsx = '/Users/fengqikai/JupyterLab/by-test/data/Final/ICTV_VH.checked.segments.gbk_info.xlsx'
#segments_info = pd.read_excel(segments_xlsx)


# 2. 函数定义
## 2.1 查询数据库
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

## 2.2 筛选出有效的比对结果
# 获取与max bitscore相差不超过(100-related_p)%的比对信息，其他的认为是无效比对
def get_useful_blast_info(blast_out=None, related_p=100):
    out_list = []
    for contig_id in tblastx_out['qseqid'].drop_duplicates():
        tmp_df = tblastx_out.query('qseqid == @contig_id')
        if int(tmp_df.qlen.drop_duplicates().max()) > 500: # 增加了长度筛选
            threshold = (related_p/100) * tmp_df['bitscore'].max()
            tmp_use = tmp_df[tmp_df['bitscore'] >= threshold]
            out_list.append(tmp_use)
        else:
            pass
    return pd.concat(out_list)

## 2.3 判断挖掘到的基因组中是否包含节段病毒的RdRp
def with_or_without_RdRp(Virome=None, RdRp_contigs=None):
    with_RdRp_Virome = dict()
    without_RdRp_Virome = dict()
    for taxid in Virome.keys():
        Virome_segments = set(Virome[taxid][0])
        if Virome_segments.intersection(RdRp_contigs):
            with_RdRp_Virome[taxid] =  list(Virome_segments)
        else:
            without_RdRp_Virome[taxid] = list(Virome_segments)
    return {'with_RdRp':with_RdRp_Virome, 'without_RdRp':without_RdRp_Virome}

## 2.4 节段数完整性判断
# 输入基因组**信息，判断构建的节段数相差多少
# 分类层级未确定
def segmentsNum_check(Virome=None):
    #known_info = segmented_genome_info['Family', 'Genus', 'Species', 'Taxid']
    known_info = segmented_genome_info['Family', 'Genus', 'Species', 'segmentsFamily', 'segmentsGenus']
    Family = Virome  #需要想好找到的基因组及其分类信息的存储格式 
    Genus = Virome
    return 



# 3. 基因组挖掘
## 3.0 比对结果预处理
### 3.0.1 获取有效的比对结果
blastx_use = get_useful_blast_info(blast_out=blastx_out, related_p=80)
### 3.0.2 获取比对上的contigs-segments
contigs_segments = blastx_use[['qseqid', 'sseqid']].drop_duplicates() # 只考虑是否有效比对，不考虑比对上了哪部分
contigs_list = list(contigs_segments['qseqid'])
segments_list = list(contigs_segments['sseqid'])
contigs_segments_list = zip(contigs_list, segments_list)

### 3.0.3 获取节段病毒核酸序列对应的Taxonomy信息
# 连接Accession2Taxonomy数据库
conn = sqlite3.connect(prot_accession2taxonomy_db)
cur = conn.cursor()
# 获取比对上的segments的分类信息
# 也就是该segments潜在的分类信息
prot_taxonomy_list = []
not_found_list = []
for contig_id, segment_id in contigs_segments_list:
    found = accession2taxonomy(segment_id)
    try :
        #len(found)
        found = list(found)
        found.insert(0, contig_id)
        prot_taxonomy_list.append(found)
    except:
        not_found_list.append(segment_id)
    header = ['contig_id', 'Accession', 'Taxid', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Virus_name', 'Lineage']
    prot_taxonomy_df = pd.DataFrame(data=prot_taxonomy_list, columns=header)

# 关闭游标和连接
cur.close()
conn.close()

prot_taxonomy_df.fillna('Unclassified', inplace=True)


### 3.0.4 整合比对信息和分类信息
contigs_segments_blastx_taxonomy = pd.merge(blastx_use, prot_taxonomy_df, how='left', left_on=['qseqid', 'sseqid'], right_on=['contig_id','Accession'])

### 3.0.5 找比对上的基因组
# 以taxid来区分基因组; 还没有验证这是否可行
Taxid_list = list(contigs_segments_blastx_taxonomy['Taxid'].drop_duplicates())

## 3.1 选比对上的基因组的每条节段的最优比对
### 3.1.1 找潜在的基因组
found_Virome = dict()
for taxid in Taxid_list:
    tmp = contigs_segments_blastx_taxonomy[contigs_segments_blastx_taxonomy['Taxid'] == taxid] # 该基因组的全部比对结果
    ref_id_list = list(tmp['Accession'].drop_duplicates()) # 该基因组全部比对上的节段
    seq_ids = []
    # 遍历该基因组比对上的节段
    for ref_id in ref_id_list:
        rtmp = tmp[tmp['Accession'] == ref_id] # 该节段的比对信息
        # 该节段的最优比对(max_bitscore)
        besthit_info = rtmp.loc[rtmp['bitscore'].idxmax(),:]
        # 构建基因组
        seq_ids.append(besthit_info['qseqid'])
    seq_info = besthit_info[['Taxid', 'Accession', 'Family', 'Order', 'Genus', 'Virus_name']]
    found_Virome[taxid] = [seq_ids, seq_info]

### 3.1.2 输出找到的基因组
out_virome = with_or_without_RdRp(Virome=found_Virome, RdRp_contigs=RdRp_contigs)

with open(virome_out_besthit, 'w') as file:
    json.dump(out_virome, file)


## 3.2 凡比对上的节段均考虑

##### 统计比对上的基因组的分类信息
contig_uniq = prot_taxonomy_df['contig_id'].drop_duplicates()
Taxid_uniq = prot_taxonomy_df['Taxid'].drop_duplicates()
Order_uniq = prot_taxonomy_df['Order'].drop_duplicates().dropna()
Family_uniq = prot_taxonomy_df['Family'].drop_duplicates().dropna()
Genus_uniq = prot_taxonomy_df['Genus'].drop_duplicates().dropna()
virusName_uniq = prot_taxonomy_df['Virus_name'].drop_duplicates()
lst = [contig_uniq, Taxid_uniq, Order_uniq, Family_uniq, Genus_uniq, virusName_uniq]
lst_name = ['contig_uniq', 'Taxid_uniq', 'Order_uniq', 'Family_uniq', 'Genus_uniq', 'virusName_uniq']
for tmp_name, tmp in zip(lst_name, lst):
    print(tmp_name, ':', len(tmp))


### 3.2.1 找潜在的基因组
#### 3.2.1.1 在Taxid层级划分基因组
Taxid_Virome = dict()
Taxid_out_list = []
for taxid in Taxid_uniq:
    tmp = prot_taxonomy_df[prot_taxonomy_df['Taxid'] == taxid]
    seq_id = list(tmp['contig_id'].drop_duplicates())
    seq_info = tmp[['Taxid', 'contig_id', 'Accession', 'Family', 'Order', 'Genus', 'Virus_name']].drop_duplicates()
    Taxid_Virome[taxid] = [seq_id, seq_info]
    Taxid_out_list.append(seq_info.drop('Accession', axis=1).drop_duplicates()) #避免同一序列比对上不同segments，但是这些segments的分类信息却相同的情况
    Taxid_out_info = pd.concat(Taxid_out_list)


#### 3.2.1.2 在Genus层级划分基因组
Genus_Virome = dict()
Genus_out_list = []
for genus in Genus_uniq:
    tmp = prot_taxonomy_df[prot_taxonomy_df['Genus'] == genus]
    seq_id = list(tmp['contig_id'].drop_duplicates())
    seq_info = tmp[['Taxid', 'contig_id', 'Accession', 'Family', 'Order', 'Genus', 'Virus_name']].drop_duplicates()
    Genus_Virome[genus] = [seq_id, seq_info]
    Genus_out_list.append(seq_info.drop(['Accession', 'Taxid', 'Family', 'Order', 'Virus_name'], axis=1).drop_duplicates()) #避免同一序列不同的Accession的情况
Genus_out_info = pd.concat(Genus_out_list)


#### 3.2.1.3 在Family层级划分基因组
Family_Virome = dict()
Family_out_list = []
for family in Family_uniq:
    tmp = prot_taxonomy_df[prot_taxonomy_df['Family'] == family]
    seq_id = list(tmp['contig_id'].drop_duplicates())
    seq_info = tmp[['Taxid', 'contig_id', 'Accession', 'Family', 'Order', 'Genus', 'Virus_name']].drop_duplicates()
    Family_Virome[family] = [seq_id, seq_info]
    Family_out_list.append(seq_info.drop(['Accession', 'Taxid', 'Order', 'Genus', 'Virus_name'], axis=1).drop_duplicates()) #避免同一序列不同的Accession的情况
Family_out_info = pd.concat(Family_out_list)
without_RdRp_Virome = with_or_without_RdRp(Virome=Taxid_Virome, RdRp_contigs=RdRp_contigs)['without_RdRp']

### 3.2.2 输出找到的基因组
##### 分别在strain，Genus，Family层级上构建了基因组，并过滤掉不包含RdRp的基因组
out_virome = {'Taxid' : with_or_without_RdRp(Virome=Taxid_Virome, RdRp_contigs=RdRp_contigs), 
              'Genus': with_or_without_RdRp(Virome=Genus_Virome, RdRp_contigs=RdRp_contigs), 
              'Family': with_or_without_RdRp(Virome=Family_Virome, RdRp_contigs=RdRp_contigs)}

with open(virome_out_anyhit, 'w') as file:
    json.dump(out_virome, file)

##### 所有有效比对的contigs-segments对应的分类信息

out_info = Taxid_out_info.sort_values(['Order', 'Family', 'Genus', 'Taxid', 'Virus_name'])[['Order', 'Family', 'Genus', 'Taxid', 'Virus_name', 'contig_id']]
out_info.to_csv(out_file, sep='\t', index=False)
