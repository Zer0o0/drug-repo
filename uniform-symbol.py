# uniform gene symbols
#
HGNC_NOW=[]
HGNC_BF=[]

with open(r'data\hgnc-biomart.txt','r') as f:
    for line in f:
        line=line.strip().split('\t')
        try:
            HGNC_BF.append(line[4])
            HGNC_NOW.append(line[2])
        except:
            pass


def uniform_gene_symbol(gs):
    gs=str(gs)
    if gs in HGNC_NOW:
        return gs
    elif gs in HGNC_BF:
        return HGNC_NOW[HGNC_BF.index(i)]
    else:
        print("the gene symbol is not exist!")


if __name__=='__main__':
    pass

