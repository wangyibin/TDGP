#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import os
import os.path as op
import pandas as  pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

def import_fpkm(fpkm):
    names = ['geneA', 'geneB', 'A1', 'A2', 'A3',                                       
                    'B1', 'B2', 'B3']
    df = pd.read_csv(fpkm, sep='\t', header=None, 
                                index_col=None,)# names=names)
    df = df.dropna(axis=1, how='all')
    df.columns = names
    
    return df

def getLog2FC(row):
    meanA = row.loc[['A1', 'A2', 'A3']].mean()
    meanB = row.loc[['B1', 'B2', 'B3']].mean()
    log2fc = np.log2(meanA + 1) - np.log2(meanB + 1)
    
    return log2fc 

def getPvalue(row):
    valueA = row.loc[['A1', 'A2', 'A3']]
    valueB = row.loc[['B1', 'B2', 'B3']]
    pvalue = stats.ttest_ind(valueA, valueB).pvalue 
                        
    return pvalue if not np.isnan(pvalue) else 0


def main(args):
    fpkm, out = args

    df = import_fpkm(fpkm)
    df['log2fc'] = df.apply(getLog2FC, axis=1)
    df['pvalue'] = df.apply(getPvalue, axis=1)
    df['FDR'] = fdrcorrection(df.pvalue)[1]
    df.to_csv(out, header=True, index=None, 
                    sep='\t')#, float_format='%f')

if __name__ == "__main__":
    if len(sys.argv) <= 1:
        print("Usage: alleleDEG.py <fpkm> <out>")
        sys.exit()
    main(sys.argv[1:])
