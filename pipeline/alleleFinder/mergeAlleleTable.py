#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
merge allele table from synteny and gmap, and filter by blast.
"""

import argparse
from logging import warn
import os
import os.path
import sys
import logging
import warnings
import pandas as pd
import numpy as np
import multiprocessing as mp
from pandarallel import pandarallel
from joblib import Parallel, delayed

from collections import OrderedDict, Counter
from utils import AllFrame2alleleTable, import_allele_table, remove_dup_from_allele_table
from utils import import_blast, alleleTable2AllFrame

warnings.filterwarnings("ignore")

def create_empty_allele_all_df(df):
    columns = df.columns
    empty_df = pd.DataFrame(columns=columns)
    return empty_df


def filter_by_blast(row, allele_table_df, gene_headers):

    gene1, gene2 = row.qseqid, row.sseqid
    tmp_df = allele_table_df[allele_table_df.isin([gene1, gene2]).any(1)]
    dup_gene_df = tmp_df.drop(gene_headers, axis=1)
    if len(tmp_df) < 2:
        if tmp_df.empty:
            return {0: [gene1, gene2] + [np.nan] * (len(dup_gene_df.columns) - 2) \
                   + [np.nan] * len(gene_headers)}
        else:
            idx, = tmp_df.index.to_list()
            dup_gene_nums = len(dup_gene_df.dropna())
            tmp_df_gene_list = set(tmp_df.loc[idx].to_list() + [idx])
            if gene1 not in tmp_df_gene_list:
                tmp_df.iloc[:, dup_gene_nums] = gene1
            if gene2 not in tmp_df_gene_list:
                tmp_df.iloc[:, dup_gene_nums] = gene2

            return {idx: tmp_df.loc[idx].to_list()}
    else:
        dup_gene_headers = dup_gene_df.columns
        idx_series = tmp_df.apply(lambda x: len(
            x.dropna()), 1).sort_values(ascending=False)
        idx_series = idx_series.index.to_list()
        idxmax = idx_series[0]
        idxmin = idx_series[1:]
        new_dup_genes = []
        for l in tmp_df.loc[idxmin].values.tolist():
            new_dup_genes += list(filter(lambda x: not isinstance(x, float), l))
        old_dup_genes = dup_gene_df.loc[idxmax].dropna().to_list()

        new_dup_genes = set(new_dup_genes) | set(old_dup_genes)
        new_dup_genes = sorted(new_dup_genes)
        tmp_df.loc[idxmax, dup_gene_headers] = dict(
            zip(dup_gene_headers, new_dup_genes))
#         tmp_df = tmp_df.loc[idxmax]
        idx = idxmax
    res_list = tmp_df.loc[idx].to_list()
    res_db = {idx: res_list}
    return res_db


def mergeAlleleTable(args):
    """
    %(prog)s 
    """

    p = p = argparse.ArgumentParser(prog=mergeAlleleTable.__name__,
                                    description=mergeAlleleTable.__doc__,
                                    formatter_class=argparse.RawTextHelpFormatter,
                                    conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('alleletable1', help='allele table from synteny')
    pReq.add_argument('alleletable2', help='allele table from gmap')
    pReq.add_argument('blast',
                      help='blast file after filter')
    pOpt.add_argument('--gene_headers', nargs="+",
                      default=['geneA', 'geneB', 'geneC', 'geneD'],
                      help='gene headers of table [default: %(default)s]')
    pOpt.add_argument('-t', '--threads', type=int, default=4,
                      help='threads numnber of program [default: %(default)s]')

    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
                      default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
                      help='show help message and exit.')

    args = p.parse_args(args)

    synteny_allele_table = import_allele_table(args.alleletable1)
    gmap_allele_table = import_allele_table(args.alleletable2)
    blast_df = import_blast(args.blast, onlyid=True)

    merge_df = pd.concat(
        [synteny_allele_table, gmap_allele_table], axis=0, ignore_index=True,)
    merge_df.reset_index(drop=True, inplace=True)

    merge_all_df = alleleTable2AllFrame(merge_df)

    #rmdup_df = remove_dup_from_allele_table(merge_all_df)

    # pandarallel.initialize(nb_workers=args.threads, verbose=0)
    # filter_df = blast_df.parallel_apply(filter_by_blast, axis=1, 
    #                 args=(rmdup_df, args.gene_headers))
    # res_df = pd.concat(filter_df.to_frame()[0].map(
    #     lambda x: pd.DataFrame.from_dict(x)).to_list(), axis=1).T.dropna(how='all', axis=0)
    # res_df.columns = rmdup_df.columns
    # res_df = res_df.drop_duplicates()
    # res_df.reset_index(inplace=True, drop=True)
    # rmdup_res_df = remove_dup_from_allele_table(res_df)
    res_df = AllFrame2alleleTable(merge_all_df)
    res_df.to_csv(args.output, sep='\t', header=True, 
                index=None, na_rep='.')

if __name__ == "__main__":
    mergeAlleleTable(sys.argv[1:])

