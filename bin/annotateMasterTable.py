#!/usr/bin/env python

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
    parser = ap.ArgumentParser()
    parser.add_argument('-mt', '--masterTable', required = True,
                        help = '''tab-separated table holding at least initiation site
                                  coordinates and name (i.e. chr, start, end, name)''')
    parser.add_argument('-a', '--annotation', required = True,
                        help = 'entrez mapped ChIPseeker annotation table')
    parser.add_argument('-o', '--outputFile', required = True,
                        help = 'name of the output file')
    args = parser.parse_args()

    logging.info('reading masterTable %s' % args.masterTable)
    mt = pd.read_csv(args.masterTable, sep = '\t')
    mt.sort_values(by=['chr', 'start'], inplace = True)
    mt.reset_index(drop = True, inplace = True)

    logging.info('reading mapped annotation table %s' % args.annotation)
    anno = pd.read_csv(args.annotation, sep = '\t')
    anno = anno.loc[:, ['seqnames', 'start', 'end', 'V4', 'annotation', 'genesymbol']]
    anno.columns = ['chr', 'start', 'end', 'name', 'annotation', 'symbol']
    anno.loc[anno.symbol.isna(), ['annotation', 'symbol']] = ['Intergenic', '-']
    anno.sort_values(by=['chr', 'start'], inplace =True)
    anno.reset_index(drop = True, inplace = True)
    anno.loc[:, 'annotation'] = anno.annotation.apply(unifyAnno)

    logging.info('merging tables')
    mt = mt.merge(anno.loc[:, ['name', 'annotation', 'symbol']],
                  how = 'left', on = 'name')

    logging.info('writing output table')
    mt.to_csv(args.outputFile, sep = '\t', index = False)

    for feature in mt.annotation.unique():
        mt[feature] = 0
        mt.loc[mt.annotation == feature, feature] = 1

    mt.to_csv(args.outputFile + '.plotTable', sep = '\t', index = False)
