'''
--------------------------------------------------------------------------
                        ***IGSEA analysis***
--------------------------------------------------------------------------
Perturbations expression data sets were retrived from LINCS project(2017).
Contains: 
1) moleculars: 1281 drugs form DrugBank
2) cell lines: A375,A549,HA1E,HCC515,HEPG2,HT29,MCF7,PC3,VCAP,NEU,NPC 
3) doses: 1uM,2uM,3uM,5uM,10uM,<1uM 
4) times: 3h,6h,24h,48h,72h,96h
5) genes: 978+9196=10174
Usage example:
python igsea.py -o outdir -g disease.grp -c NEU,NPC
Session:
python3, pandas, gsea-3.0.jar, jre-8
--------------------------------------------------------------------------
'''

from optparse import OptionParser
import pandas as pd
import re
import time
import sys
import os
print(__doc__)


usage = 'usage: %prog [options] -h<--help> -g<--geneset> -o<--output> -c<--cell>'
version = '%prog 1.0'
optParser = OptionParser(usage=usage, version=version)
optParser.add_option('-g', '--geneset', action='store', type='string', dest='GENESET',
                     help='risk geneset of the target disease')
optParser.add_option('-o', '--output', action='store', type='string', dest='OUTDIR',
                     help='write output to directory',
                     default='igsea')
optParser.add_option('-c', '--cell', action='store', type='string', dest='CELL',
                     help='choose cell lines',
                     default='A375,A549,HA1E,HCC515,HEPG2,HT29,MCF7,PC3,VCAP,NEU,NPC')

(options, args) = optParser.parse_args()

grp = options.GENESET              # disease genes
outdir = options.OUTDIR            # output directory
cell = options.CELL.upper()        # cell lines
try:
    cell = cell.split(',')
except:
    cell = [cell]
CELLS = ['A375', 'A549', 'HA1E', 'HCC515', 'HEPG2',
         'HT29', 'MCF7', 'PC3', 'VCAP', 'NEU', 'NPC']

if len(args) != 0:
    optParser.error('Incorrect number of arguments!')
if not grp:
    optParser.error('Can not find geneset file of your target disease!')
if not all([i in CELLS for i in cell]):
    optParser.error('Choose supported cell lines!')

print('geneset:', grp)
print('outdir:', outdir)
print('cell lines:', cell)
print()

cwd = os.getcwd()                         # work directory
outdir = os.path.join(cwd, outdir)        # output directory
root = sys.path[0]                        #

try:
    os.mkdir(outdir)
except:
    pass


def gsea(grp, rnk, dir):  # gene set enrichment analysis
    jar = os.path.join(root, 'gsea-3.0.jar')
    command = 'java -cp %s -Xmx512m xtools.gsea.GseaPreranked \
        -gmx %s \
        -norm meandiv -nperm 1000 \
        -rnk %s \
        -scoring_scheme weighted -rpt_label my_analysis -create_svgs false -make_sets true \
        -plot_top_x 20 -rnd_seed timestamp -set_max 500 -set_min 15 -zip_report false \
        -out %s \
        -gui false' % (jar, grp, rnk, dir)
    os.system(command)


def run_igsea():
    genes_info = os.path.join(root, r'data\lincs\genes_info.txt')
    genes_info_ = pd.read_csv(genes_info, sep='\t', index_col=0)
    for i in cell:
        print('cell:', i)
        j = r'lincs\profile_'+i+'.txt'
        profile = os.path.join(root, j)
        profile_ = pd.read_csv(profile, sep='\t', index_col=0)

        genes_info_index = genes_info_.reindex(profile_.index)
        genes_symbol = genes_info_index['pr_gene_symbol']
        genes_symbol.name = '#genes_symbol'
        profile_.index = genes_symbol

        perts = list(profile_.columns)
        outdir_ = os.path.join(outdir, i)
        try:
            os.mkdir(outdir_)
        except:
            pass
        for p in perts:
            f = p.replace(':', '__').replace('-', '___')
            f_ = os.path.join(outdir_, f)
            try:
                os.mkdir(f_)
            except:
                pass
            rnk = os.path.join(f_, f+'.rnk')
            profile_[[p]].to_csv(rnk, sep='\t')
            gsea(grp, rnk, f_)
            print('-'*20)


def get_value(html):  # extract result
    dataset = re.search(
        r'<tr><td>Dataset</td><td>(.*?)</td></tr>', html).group(1)
    es = re.search(
        r'<tr><td>Enrichment Score \(ES\)</td><td>(.*?)</td></tr>', html).group(1)
    nes = re.search(
        r'<tr><td>Normalized Enrichment Score \(NES\)</td><td>(.*?)</td></tr>', html).group(1)
    pval = re.search(
        r'<tr><td>Nominal p-value</td><td>(.*?)</td></tr>', html).group(1)
    res = pd.Series([es, nes, pval], index=['es', 'nes', 'pval'])
    res.name = dataset
    return res


def run_accessation():
    target = grp+'.html'
    ann = os.path.join(root, r'data\lincs\pert_info.txt')
    ann_ = pd.read_csv(ann, sep='\t', index_col=0)

    def searchFile(root, target):
        items = os.listdir(root)
        for item in items:
            path = os.path.join(root, item)
            if path.split('\\')[-1] == target:
                targets.append(path)
            elif os.path.isdir(path):
                searchFile(path, target)
            else:
                pass

    for i in os.listdir(outdir):
        path = os.path.join(outdir, i)
        targets = []
        vals = []
        searchFile(path, target)
        for html in targets:
            with open(html, 'r') as f:
                text = f.read()
            val = get_value(text)
            vals.append(val)
        res = pd.concat(vals, axis=1).T
        sig_id = [i.replace('___', '-').replace('__', ':') for i in res.index]
        anni = ann_.loc[sig_id, :]
        res.index = sig_id
        res_ = pd.concat([anni, res], axis=1)
        res_.to_csv(os.path.join(
            outdir, 'result_of_igsea_'+i+'.txt'), sep='\t')


if __name__ == '__main__':
    try:
        print('calculating...')
        start = time.time()
        run_igsea()
        run_accessation()
        end = time.time()
        print('finish，time consuming/min：', (end-start)//60)
    except:
        print('check the session environment please!')
        exit(1)
