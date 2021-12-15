#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
identify allele pairs expression type, including consistent and inconsistent

"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import numpy as np
import pandas as pd 

from collections import OrderedDict

def import_table(table):
    df = pd.read_csv(table, sep='\t', header=0,
                     index_col=None)
    return df 

def getExpType_single(value):
    if value > log2fc[1]:
        return "A>B"
    elif value < log2fc[0]:
        return "A<B"
    else:
        return np.nan

def getExpType(df, log2fc=(-2, 2), pvalue=0.05, FDR=0.05):
    globals()['log2fc'] = log2fc
    sign_df = df[(df.pvalue < pvalue) & (df.FDR < FDR)]
    deg_df = sign_df[(sign_df.log2fc < log2fc[0]) | (sign_df.log2fc > log2fc[1])]
    deg_df = deg_df.reset_index(drop=True)
    deg_df['type'] = deg_df.log2fc.map(getExpType_single).astype('category')
    deg_df['genes'] = deg_df.apply(lambda x: "{}-{}".format(x.geneA, x.geneB), 
                                    axis=1)

    deg_df = deg_df.drop(['geneA', 'geneB'], axis=1)
    deg_df = deg_df.set_index('genes')
    return deg_df

def getOnlyType(df):
    res_df = df.drop(df.columns.difference(['type']), axis=1)
    return res_df

def getConsistentOrInconsistent_single(row):
    types = row.tolist()
    types = list(filter(lambda x: pd.notna(x), types))
    if len(types) <= 1:
        return np.nan
    types = set(types)
    if len(types) == 1 and row.notna().all():
        return 'Consistent'
    elif len(types) == 2:
        return "Inconsistent"
    else:
        return np.nan
    
def getConsistentOrInconsistent(deg_df_list, deg_names):
    onlyTypeDf_list = list(map(getOnlyType, deg_df_list))
    # rename header
    def rename_headers(df_list, names):
        for name, df in zip(names, df_list):
            df.columns = [name]

        return df_list
    
    onlyTypeDf_list = rename_headers(onlyTypeDf_list, deg_names)
    df = pd.concat(onlyTypeDf_list, axis=1)
    df['expr_type'] = df.apply(getConsistentOrInconsistent_single, axis=1)
    df = df.dropna(subset=['expr_type'])
    
    return df

def identifyAlleleExpressionType(args):
    """
    %(prog)s <sample1.deg ...> [Options]
        identify allele pairs expression type, 
            including consistent and inconsistent

    """

    p = argparse.ArgumentParser(prog=identifyAlleleExpressionType.__name__,
                        description=identifyAlleleExpressionType.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('deg', nargs="+", help='different table')
    pOpt.add_argument('--log2fc', default="(-2, 2)", 
            help='values of log2fc [default:$(default)s]')
    pOpt.add_argument('--pvalue', default=0.05, type=float, 
            help='pvalue of different expression [default: %(default)s]')
    pOpt.add_argument('--fdr', default=0.05, type=float,
            help='fdr of different expression [default: %(default)s]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    deg_paths = args.deg
    log2fc = eval(args.log2fc)
    
    pvalue = args.pvalue
    FDR = args.fdr
    deg_names = list(map(lambda x: op.basename(x).split(".", 1)[0], 
                    deg_paths))
    df_list = list(map(import_table, deg_paths))
    deg_df_list = [getExpType(df, log2fc, pvalue, FDR) for df in df_list]
    df = getConsistentOrInconsistent(deg_df_list, deg_names)
    df.index.name = 'genes'

    df.to_csv(args.output, sep='\t', index=True, header=True, na_rep='-')

if __name__ == "__main__":
    identifyAlleleExpressionType(sys.argv[1:])