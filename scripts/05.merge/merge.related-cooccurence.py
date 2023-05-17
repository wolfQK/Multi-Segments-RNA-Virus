import os
import json
import argparse
import os.path
import pandas as pd
import networkx as nx
from collections import defaultdict, namedtuple

## 传参
parser = argparse.ArgumentParser(description='manual to this script')
parser.add_argument('--outdir', type=str, default=None)
parser.add_argument('--related_found_file', type=str, default=None)
parser.add_argument('--cluster_file', type=str, default=None)
parser.add_argument('--overlap_cutoff', type=float, default=1.0)
parser.add_argument('--RdRp_list_file', type=str, default=None)
args = parser.parse_args()

# 该样本基于近缘找到的基因组
outdir = args.outdir
related_found_file = args.related_found_file
# 聚类信息
cluster_file = args.cluster_file
# overlap_cutoff -- 该值应该会随着样本数而有所变化
overlap_cutoff = args.overlap_cutoff
# 潜在病毒序列中比对上RdRp库的RdRp list(含样本名)
RdRp_list_file = args.RdRp_list_file
RdRp_list = list(pd.read_csv(RdRp_list_file, sep='\t')['sample~contig'])

## 解析cd-hit-est聚类结果
Member = namedtuple('Member', ['sample', 'contig', 'length', 'percent_id', 'percent_id_sign', 'coverage'])
def parse_cdhit_row(row):
    if '*' in row:
        index, length, name, percent_id = row.split()
        percent_id_sign, percent_id = '0', 100
    else:
        index, length, name, _, percent_id = row.split()
    length = int(length.strip(',nt'))
    name = name.strip('>').strip('.')
    sample, contig = name.split('~')
    coverage = float(contig.split('_')[-1])
    
    if percent_id != 100:
        percent_id_sign, percent_id = percent_id.strip('%').split('/')
        percent_id = float(percent_id)
    return Member(sample=sample, contig=contig, length=length, percent_id=percent_id,
                  percent_id_sign=percent_id_sign, coverage=coverage)

## 以dict的形式记录每个cluster的信息
# key为clusterId, value为一个list, list的内容为namedtuple(每条序列的聚类信息)
all_clusters = defaultdict(list)
with open(cluster_file, 'r') as file:
    for line in file:
        if line.startswith('>Cluster'):
            cluster_id = line.split()[-1]
        else:
            member = parse_cdhit_row(line)
            all_clusters[cluster_id].append(member)

## 去除只有一条序列的cluster
clusters = {cluster_id: all_clusters[cluster_id] for cluster_id in all_clusters if len(all_clusters[cluster_id]) > 1}

## 构建contig - cluster的对应关系
# 这里默认在一个cluster下来自不同样本的contig的名称没有重复
'''
这里还是用标记了样本来源的contigs_id
'''
all_contig2cluster = {member.sample + '~' + member.contig: cluster_id for cluster_id in all_clusters for member in all_clusters[cluster_id]}
contig2cluster = {member.sample + '~' + member.contig: cluster_id for cluster_id in clusters for member in clusters[cluster_id]}

## 查看涉及到的样本 
samples = list(set([member.sample for cluster_id in clusters for member in clusters[cluster_id]]))

## 构建一个无向图，图中的顶点为cluster和sample，边表示cluster中在哪些sample中被检出
def graph_from_clusters(clusters):
    G = nx.Graph() #创建一个无向图
    for cluster_id in clusters:
        for member in clusters[cluster_id]: #遍历每个cluster下的每个member
            if member.length > 500: #筛选出大于500bp的序列
                G.add_edge(cluster_id, member.sample, attr_dict=member._asdict()) #向图中添加边：Node1 Node2 attr(该contig的聚类信息)
    return G
# 构建无向图:同一个cluster下长度大于500的cluster_id和对应的sample构成一条边
# 一个cluster会对应多个样本，一个样本中也会检出多个cluster
G = graph_from_clusters(clusters) 

##  输入一系列节点，返回他们的所有邻居节点 
def walk(nodes, G=G):
    if not isinstance(nodes, list): #将节段转换成列表
        nodes = [nodes]
    return [nbr for node in nodes for nbr in G.neighbors(node) ] 

## 找某节点在一定深度上的所有邻居节点
# cluster节点的邻居即是sample节点,反之亦然
def nbhd(start, depth=1):
    if isinstance(start, str):
        start = [start]
    n = start #起始节点
    for i in range(depth): #找距离为depth的邻居
        n = n + walk(n) 
    return set(n)
'''
至此便获得了sample与cluster的关系, 可通过调用图中的节点的一级邻居来查看
'''

## 识别dark contigs
# 基于一些样本和一些病毒簇（subset），生成一个表格：sample * cluster
# 表格的值为01或sum(coverage)或sum(length)
'''
subset是一个包含了特定cluster_id(函数定义时是可以输入多个,实际每次应该只输入一个)和sample的list
该list包含:(RdRp序列对应的cluster1, RdRp_cluster中的sample_set1, 检出RdRp的sample_set1中的contigs所在的cluster2[这些cluster2才可能与RdRp所在的cluster1共现], 这些cluster2中包含的sample_set2[用于计算overlap])
'''
def df_from_subset(subset):
    df = pd.DataFrame(columns = [cluster_id for cluster_id in clusters if cluster_id in subset], 
                      index = [sample for sample in samples if sample in subset], dtype=int).fillna(0)
    for cluster_id in clusters:
        if cluster_id not in subset:
            continue
        for member in clusters[cluster_id]:
            if member.sample not in subset:
                continue
            df.loc[member.sample, cluster_id] = 1
    return df

## 基于近缘找到的节段病毒
# 基于共现进行补充
def fish_dark_matter(bait_contigs, overlap_cutoff=0.9):
    '''
    潜在节段挖掘；
    这里的bait_contigs应该是基于比对已经识别的contigs, 需要的是找出dark contigs中可能是潜在节段的contigs
    '''
    bait_clusters = [contig2cluster[contig] for contig in bait_contigs if contig in contig2cluster] #bait_contig是RdRp对应的contigs，返回RdRp所在的clusters ------------------------------有问题-------
    neighborhood = nbhd(bait_clusters, depth = 3) #与RdRp contigs相距为3的节点
    df = df_from_subset(neighborhood) #与RdRp contigs距离为3的节点构成的sample * cluster表格；通过设置depth确实可以在一定程度上减少涉及的samples；减少计算量
    #depth=1时可获得该cluster的sample sets; depth=2时可获得sample sets中的sample中检出的Clusters（与该cluster共现的cluster就来自这些Clusters）; depth=3时可获得这些Clusters对应的sample sets; 然后就可以比较sample sets的overlap了
    
    samples_with_bait = df.loc[:, bait_clusters].sum(axis = 1) == len(bait_clusters)  #检出所有bait_contigs（通常只有1个）的样本
    n_samples_with_bait = df.loc[:, bait_clusters].sum(axis = 0).mean() #该bait_contigs在多少个样本中检出 -- 忽略了sample sets的组成
    
    clusters_containing_bait = df.loc[samples_with_bait].sum(axis = 0) >= (overlap_cutoff * n_samples_with_bait) # 检出该bait_contig的样本
    clusters_not_overflowing_bait = (df.sum() * overlap_cutoff) <= n_samples_with_bait # 检出样本数*overlap_cutoff小于检出bait_contig的样本数 -- 广泛存在（既在bait_contigs对应的sample中存在也在其他sample中存在）的contig_cluster可能是污染
    overlapping_clusters = df.columns[clusters_containing_bait & clusters_not_overflowing_bait] 
        
    df = df[overlapping_clusters]
    df = df.loc[df.sum(axis = 1) > 0,:]

    return df

# 输入近缘的结果,处理后输出新的挖掘结果
def related2clustering(file=related_found_file, level='Taxid'):
    with open(related_found_file, 'r') as load_f:
        data = json.load(load_f)
    for sample in data.keys():
        sample_virome = data.get(sample) #该样本在不同分类层级下获得的各种基因组
        # 这里仅考虑包含RdRp的基因组
        virome = sample_virome.get(level)
        with_RdRp_virome = virome.get('with_RdRp')
        dirName = outdir + '/' + sample + '/' + level + '_' + str(overlap_cutoff)
        if not os.path.exists(dirName):
            os.makedirs(dirName, )
        for tax, contigs in with_RdRp_virome.items():
            samples_clusters = fish_dark_matter(contigs)
            # 输出
            if len(samples_clusters.columns) > 1:
                merged_out_list = dict()
                baseName = sample + '.' + tax + '.segments.json'
                merged_out_file = os.path.join(dirName, baseName)
                merged_sample_contigs = set()
                tmp = samples_clusters.loc[sample]
                occur_clusters = tmp == 1
                use_tmp = tmp[occur_clusters]
                use_clusters = use_tmp.index
                for cluster_id in use_clusters:
                    members = clusters[cluster_id]
                    for member in members:
                        if member.sample == sample:
                            merged_sample_contigs.add(member.sample + '~' + member.contig)

                new_contigs = list(merged_sample_contigs.intersection(set(contigs)))
                merged_out_list[sample] = list(merged_sample_contigs)
                out_list = {'all_contigs':merged_out_list, 'new_contigs':new_contigs}
                with open(merged_out_file, 'w') as merged_out:
                    json.dump(out_list, merged_out)
            else:
                pass    
    return 0

# 在科属种层级上分别对近缘找到的节段进行补充
# 输入需要是基于近缘基因组找到的基因组
related2clustering(file=related_found_file, level='Taxid')
related2clustering(file=related_found_file, level='Genus')
related2clustering(file=related_found_file, level='Family')
