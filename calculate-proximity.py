# caculate proximity between drugs and disease in the interactome network
#
import sys
import time
import math
import networkx as nx
from networkx.algorithms.shortest_paths.generic import shortest_path_length
from networkx.classes.function import degree, subgraph, nodes
from networkx.exception import NetworkXNoPath, NetworkXError
import random
from random import sample
import pandas as pd 
import matplotlib.pyplot as plt

def read_PPI(fp):
    ppi = []
    with open(fp, 'r') as f:
        for line in f:
            if '#' not in line:
                line = line.strip().split('\t')
                node1 = line[0]
                node2 = line[1]
                ppi.append((node1, node2))
    return ppi


def read_disease_gene(fp):
    geneset = []
    with open(fp, 'r') as f:
        for line in f:
            if '#' not in line:
                gene = line.strip()
                geneset.append(gene)
    return geneset


def read_drugs_info(fp):
    drug = {}
    with open(fp, 'r',encoding='utf-8') as f:
        for line in f:
            if '#' not in line:
                line = line.strip().split('\t')
                targets = line[1].split(',')
                targets_ = [i.upper() for i in targets if i]
                drug[line[0]] = targets_
    return drug


def weight(fp):
    w = {}
    with open(fp, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            w[line[0]] = math.log(int(line[1])+1)
    return w

def calculate_distance(graph, sourceset, targetset, weight):
    res = 0
    n = len(sourceset)
    for s in sourceset:
        ds = []
        for t in targetset:
            d = shortest_path_length(graph, s, t)
            if d == 0:  # the target is one of the disease genes
                ds.append(d-weight[t])
            else:  # the target is not a disease gene
                ds.append(d)
        res += min(ds)
    distance = res/n
    return distance


if __name__ == '__main__':
    fp_ppi = r'data\ppin.txt'
    fp_drug = r'data\drug_keep.txt'
    fp_gene = r'data\gene.txt'
    fp_wt = r'data\weight.txt'
    fp_res = r'data\d_drug.txt'
    fp_ref = r'data\d_ref.txt'

    PPI = read_PPI(fp_ppi)
    DRUG = read_drugs_info(fp_drug)
    GSET = read_disease_gene(fp_gene)

    print('start calculating...')
    # construct a protrin-protein network
    G = nx.Graph()
    G.add_edges_from(PPI)

    W=weight(fp_weight)
    RES = {}
    for i in DRUG:
        print(i)
        d = clostD(G, DRUG[i], GSET, W)
        RES[i] = d
    with open(fp_res, 'w') as f:
        for k,v  in RES.items():
            line = str(k)+'\t'+str(v)+'\t'+str(len(DRUG[k]))+'\n'
            f.write(line)
    # calculate reference distace
    print('reference distance...')
    random.seed(667)
    seq={i:int(random.weibullvariate(1,0.5))+1 for i in range(10000)}
    DRUG_REF = {i: sample(nodes(G), j) for i,j in seq.items()}

    REF = {}
    for i in DRUG_REF:
        d = clostD(G, DRUG_REF[i], GSET, W)
        REF[i] = d
    with open(fp_ref, 'w') as f:
        for k,v in REF.items():
            line = str(k)+'\t'+str(v)+'\t'+str(len(DRUG_REF[k]))+'\n'
            f.write(line)
    print('finish!')
