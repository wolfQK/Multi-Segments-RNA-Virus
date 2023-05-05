import json
import os.path
#import sys
import pandas as pd
import networkx as nx
from collections import defaultdict, namedtuple

'''
基于作者的contigs中筛选出来的潜在病毒序列做聚类, 构建sample_contig-cluster和cluster-protein的对应关系,以发现共现模式
基于共现模式,挖掘潜在的节段病毒序列 
cd-hit-est -M 20000 -T 4 -d 0 -c 0.95 -g 1; 这种聚类方式会让一些长度相差较大的序列聚到一起,也可能较短的序列就是较长序列的一部分
聚类参数的设置对挖掘影响也比较大,因为聚类参数直接决定着cluster中包含的样本多少,而最终的共现关系则确定的是各个cluster对应的sample_sets的overlap
新增了对广泛出现(既出现在检出RdRp的样本中也出现在其他的样本中)的cluster的过滤
'''

cluster_file = '/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Pub_data/Mosquitos-California/Analysis/clustering/Author/cd-hit/all_putative.author.ge500bp.contigs.cluster.clstr'
dirName = '/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Pub_data/Mosquitos-California/Analysis/clustering/Author/cd-hit/Clusters'
RdRp_list_file = '/hwfsxx1/ST_HN/P20Z10200N0206/USER/fengqikai/biye/Pub_data/Mosquitos-California/Analysis/results/new.author.putative.RdRp.list'
RdRp_list = list(pd.read_csv(RdRp_list_file, sep='\t')['sample~contig']) #读入文化并转成list
overlap_cutoff = 0.9

# 解析cd-hit-est聚类结果
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

# 以dict的形式记录每个cluster的信息
# key为clusterId, value为一个list, list的内容为namedtuple(每条序列的聚类信息)
all_clusters = defaultdict(list)
with open(cluster_file, 'r') as file:
    for line in file:
        if line.startswith('>Cluster'):
            cluster_id = line.split()[-1]
        else:
            member = parse_cdhit_row(line)
            if 'water' in member.sample.lower(): #跳过对照组
                continue
            all_clusters[cluster_id].append(member)

# 去除只有一条序列的cluster
clusters = {cluster_id: all_clusters[cluster_id] for cluster_id in all_clusters if len(all_clusters[cluster_id]) > 1}
# 查看涉及到的样本 
samples = list(set([member.sample for cluster_id in clusters for member in clusters[cluster_id]]))
# 构建contig - cluster的对应关系
contig2cluster = {member.sample + '~' + member.contig: cluster_id for cluster_id in clusters for member in clusters[cluster_id]}

# 构建一个无向图，图中的顶点为cluster和sample，边表示cluster中在哪些sample中被检出
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

#  输入一系列节点，返回他们的所有邻居节点 
def walk(nodes, G=G):
    if not isinstance(nodes, list): #将节段转换成列表
        nodes = [nodes]
    return [nbr for node in nodes for nbr in G.neighbors(node) ] 

# 找某节点在一定深度上的所有邻居节点
# cluster节点的邻居即是sample节点,反之亦然

def nbhd(start, depth=1):
    if isinstance(start, str):
        start = [start]
    n = start #起始节点
    for i in range(depth): #找距离为depth的邻居
        n = n + walk(n) 
    return set(n)
#至此便获得了sample与cluster的关系, 可通过调用图中的节点的一级邻居来查看


# 识别dark contigs
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

# 基于共现挖掘dark contigs
'''
这部分是挖掘的重点,但原来的代码貌似有问题,对此处进行了重写
循环输入每组bait_contigs(RdRp_cluster中的每条contigs),确定与该组RdRp共现的cluster,然后按sample来划分以确定潜在的节段病毒基因组
'''
def fish_dark_matter(bait_contigs, overlap_cutoff=overlap_cutoff):
    bait_clusters = list(set([contig2cluster[contig] for contig in bait_contigs]))  #bait_contig是RdRp对应的contigs，返回RdRp所在的clusters
    neighborhood = nbhd(bait_clusters, depth = 3) #与RdRp contigs相距为3的节点
    #depth=1时可获得该cluster的sample sets; depth=2时可获得sample sets中的sample中检出的Clusters; depth=3时可获得这些Clusters对应的sample sets; 然后就可以比较sample sets的overlap了
    df = df_from_subset(neighborhood) #与RdRp contigs距离为3的节点构成的sample * cluster表格
    
    samples_with_bait = df.loc[:, bait_clusters].sum(axis = 1) == 1 #样本中是否检出该RdRp所在的cluster
    n_samples_with_bait = sum(samples_with_bait) # sample sets size
    samples = list(df.loc[samples_with_bait, bait_clusters].index) # 检出该RdRp所在的cluster对应的sample sets
    coocur_clusters = df.loc[samples].sum() >= int(n_samples_with_bait * overlap_cutoff) # 所有可能共现的cluster在该sample sets上的检出情况 -- 需要占该sample sets的(overlap_cutoff/100)*n_sample_with_bait
    clusters_not_overflowing = (df.sum() * overlap_cutoff)  <= n_samples_with_bait
    out_df = df.loc[samples, (coocur_clusters & clusters_not_overflowing)]
    out_df = out_df.loc[out_df.sum(axis = 1) > 0,:]
    return out_df

# 遍历素有的cluster,若其中存在RdRp序列,则把该cluster下的所有序列当做bait_contigs
RdRp_clusters = set()
for cluster_id in clusters.keys():
    for member in clusters[cluster_id]:
        if member.sample + '~' + member.contig in RdRp_list:
            RdRp_clusters.add(cluster_id)
            
for cluster_id in RdRp_clusters:
    bait_contigs = []
    for member in clusters[cluster_id]:
        bait_contigs.append(member.sample + '~' + member.contig)
    samples_clusters = fish_dark_matter(bait_contigs)
    
    '''
    识别出与该RdRp cluster中出现的sample sets的overlap达到overlap_cutoff的clusters
    '''
    # 输出
    if len(samples_clusters.columns) > 1:
        merged_out_list = dict()
        baseName = 'Cluster.' + cluster_id + '.segments.json'
        merged_out_file = os.path.join(dirName, baseName)
        for sample in samples_clusters.index:
            merged_sample_contigs = []
            tmp = samples_clusters.loc[sample]
            occur_clusters = tmp == 1
            use_tmp = tmp[occur_clusters]
            use_clusters = use_tmp.index
            for cluster_id in use_clusters:
                members = clusters[cluster_id]
                for member in members:
                    if member.sample == sample:
                        merged_sample_contigs.append(member.sample + '~' + member.contig)
            merged_out_list[sample] = merged_sample_contigs
        with open(merged_out_file, 'w') as merged_out:
            json.dump(merged_out_list, merged_out)
    else:
        pass