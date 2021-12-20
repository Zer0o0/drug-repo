# caculate proximity between drugs and disease in the interactome network
#
import math
import networkx as nx
from networkx.algorithms.shortest_paths.generic import shortest_path_length
import random
from random import sample


def read_ppin(file_ppin):
    ppin = []
    with open(file_ppin, 'r') as f:
        for line in f:
            if '#' not in line:
                line = line.strip().split('\t')
                node1 = line[0]
                node2 = line[1]
                ppin.append((node1, node2))
    return ppin


def read_disease_genes(file_genes):
    geneset = []
    with open(file_genes, 'r') as f:
        for line in f:
            if '#' not in line:
                gene = line.strip()
                geneset.append(gene)
    return geneset


def read_drugs_info(file_drugs):
    drug = {}
    with open(file_drugs, 'r', encoding='utf-8') as f:
        for line in f:
            if '#' not in line:
                line = line.strip().split('\t')
                targets = line[1].split(',')
                targets_ = [i.upper() for i in targets if i]
                drug[line[0]] = targets_
    return drug


def weight(file_weight):
    w = {}
    with open(file_weight, 'r') as f:
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
            try:
                d = shortest_path_length(graph, s, t)
                if d == 0:  # the target is one of the disease genes
                    ds.append(d-weight[t])
                else:  # the target is not a disease gene
                    ds.append(d)
            except:
                # for genes that are not in the PPIN, give them a large distance value
                ds.append(10000)
        res += min(ds)
    distance = res/n
    return distance


if __name__ == '__main__':
    fppin = r'data\ppin.txt'
    fdrug = r'data\drugs.txt'
    fgene = r'data\ad_genes.txt'
    fwt = r'data\weight.txt'
    fres = r'data\distance_of_drugs.txt'
    # fref = r'data\distance_of_null.txt'

    PPI = read_ppin(fppin)
    DRUG = read_drugs_info(fdrug)
    GSET = read_disease_genes(fgene)

    print('Start calculating...')
    # construct a protrin-protein network
    G = nx.Graph()
    G.add_edges_from(PPI)

    W = weight(fwt)
    RES = {}
    for i in DRUG:
        print(i)
        d = calculate_distance(G, DRUG[i], GSET, W)
        RES[i] = d
    with open(fres, 'w', encoding='utf-8') as f:
        f.write('#Drugs'+'\t'+'Distance'+'\t'+'Num_targets'+'\n')
        for k, v in RES.items():
            line = str(k)+'\t'+str(v)+'\t'+str(len(DRUG[k]))+'\n'
            f.write(line)
    # # calculate reference distace
    # print('Reference distance...')
    # random.seed(667)
    # seq = {i: int(random.weibullvariate(1, 0.5))+1 for i in range(10000)}
    # DRUG_REF = {i: sample(nodes(G), j) for i, j in seq.items()}

    # REF = {}
    # for i in DRUG_REF:
    #     d = calculate_distance(G, DRUG_REF[i], GSET, W)
    #     REF[i] = d
    # with open(fref, 'w', encoding='utf-8') as f:
    #     for k, v in REF.items():
    #         line = str(k)+'\t'+str(v)+'\t'+str(len(DRUG_REF[k]))+'\n'
    #         f.write(line)
    # print('Finish!')
