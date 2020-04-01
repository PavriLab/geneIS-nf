import pandas as pd
import argparse as ap
import numpy as np
import logging

def unifyAnno(anno):
    rdict = {'Distal': 'Intergenic', '5\'':'Exon', '3\'': 'Exon', 'Downstream': 'Intergenic'}
    return rdict[anno.split()[0]] if rdict.get(anno.split()[0]) else anno.split()[0]

class PromoterAnnotator():
    def __init__(self, promannos):
        self.promannos = set(promannos)

    def promoteroverlap(self, featurename):
        return 'P' + featurename if featurename in self.promannos else featurename

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)
    parser = ap.ArgumentParser(description='''combines the given tables to a master table containing all information,
    if --foldchange is set, an additional column is inserted into the table denoting if the peak is up- or downregulated
    or unchanged with respect to the given column names and a fold change threshold i.e. if the second column value is larger
    than fold change threshold + column1 the peak is upragulated (i.e. 1), if the first column value is smaller than 
    fold change threshold + column2 than the peak is downregulated (i.e. -1) else unchanged (i.e. 0)''')
    parser.add_argument('-ct', '--counttable', required = True,
                        help = 'tab-separated table holding the readquantification information')
    parser.add_argument('-a', '--annotation', required = True,
                        help = 'entrez mapped ChIPseeker annotation table')
    parser.add_argument('-o', '--out', required = True,
                        help = 'name of the output file')
    parser.add_argument('-f', '--fold', default = 1.5, type = float,
                        help = 'threshold for foldchange for peak classification, is log2 transformed before application')
    parser.add_argument('-wt', '--wtcolumn', default = 'CH12_shLacZ',
                        help = 'name of the column in --counttable holding the WT counts ')
    parser.add_argument('-kd', '--kdcolumn', default = 'CH12_shMcm6',
                        help = 'name of the column in --counttable holding the KD counts ')
    parser.add_argument('--dormantcutoff', default = 2, type = float,
                        help = 'cutoff to use to define dormant origins')
    parser.add_argument('--absentcutoff', default = 2, type = float,
                        help = 'cutoff to use to define absent origins')
    args = parser.parse_args()

    logging.info('reading counttable %s' % args.counttable)
    counts = pd.read_csv(args.counttable, sep = '\t', header = None, skiprows = 1,
                         names = ['chr', 'start', 'end', args.wtcolumn, args.kdcolumn])
    counts.sort_values(by=['chr', 'start'], inplace = True)
    counts.reset_index(drop = True, inplace = True)

    logging.info('reading mapped annotation table %s' % args.annotation)
    anno = pd.read_csv(args.annotation, sep = '\t')
    anno = anno.loc[:, ['seqnames', 'start', 'end', 'V4', 'annotation', 'genesymbol']]
    anno.columns = ['chr', 'start', 'end', 'name', 'annotation', 'symbol']
    anno.sort_values(by=['chr', 'start'], inplace =True)
    anno.reset_index(drop = True, inplace = True)
    anno.loc[:, 'annotation'] = anno.annotation.apply(unifyAnno)

    logging.info('merging tables')
    anno = anno.merge(counts[['chr', 'start', 'end', args.wtcolumn, args.kdcolumn]],
                      how = 'left', on = ['chr', 'start', 'end'])

    logging.info('writing output table')

    pseudo = any(anno[[args.wtcolumn, args.kdcolumn]] == 0)

    if pseudo:
        logging.info('zero counts detected adding pseudo count of 1')

    for col in [args.wtcolumn, args.kdcolumn]:
        anno.loc[anno[col] == 0, col] = 1
        anno['log2' + col] = np.log2((anno[col]*1000000)/anno[col].sum())

    # finding the different categories
    fc = np.log2(args.fold)
    anno['class'] = 0

    # dormant
    anno.loc[anno['log2' + args.wtcolumn] < anno[anno['log2' + args.kdcolumn] < args.dormantcutoff]['log2' + args.wtcolumn].min(), 'class'] = 1

    # upregulated
    anno.loc[(anno['class'] == 0) & (anno['log2' + args.kdcolumn] > fc + anno['log2' + args.wtcolumn]), 'class'] = 2

    # downregulated
    anno.loc[(anno['class'] == 0) & (anno['log2' + args.wtcolumn] > fc + anno['log2' + args.kdcolumn]), 'class'] = 4

    # absent
    anno.loc[anno['log2' + args.kdcolumn] < anno[anno['log2' + args.wtcolumn] < args.absentcutoff]['log2' + args.kdcolumn].min(), 'class'] = 5

    # unchanged
    anno.loc[anno['class'] == 0, 'class'] = 3


    anno = anno.loc[:, ['chr', 'start', 'end', 'name', 'log2' + args.wtcolumn, 'log2' + args.kdcolumn,
                        'class', 'annotation', 'symbol']]

    anno.to_csv(args.out, sep = '\t', index = False)