#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
add unmap gene into final allele table
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd
import numpy as np
import multiprocessing as mp
from pandarallel import pandarallel
from joblib import Parallel, delayed

from collections import OrderedDict, Counter
from utils import *

def addBlastGene2FinalTable(chrom, df, blast_df, 
                    suffix_length=2):
    
    tmp_blast_df = blast_df.loc[chrom]
    tmp_blast_df = tmp_blast_df.reset_index()
    tmp_blast_df = tmp_blast_df.set_index('sseqid')
    tmp_blast_df['chr1_hap'] = tmp_blast_df.chr1.str[:-suffix_length]
    df = df.drop('chrom', axis=1)
    i = 0
    gene_headers = df.columns[
        df.columns.str.find('gene').map(lambda x: x == 0)]
    paralogs_gene_columns = df.columns[
        df.columns.str.find('paralogs').map(lambda x: x == 0)]
    for sseqid, (qseqid, chr1_hap) in tmp_blast_df[['qseqid', 
                                                'chr1_hap']].iterrows():
        tmp_df = df[df.isin([sseqid]).any(1)]
        gene_header = 'gene{}'.format(which_hap(sseqid))
        if gene_header == 'geneNone':
            continue
        if tmp_df.empty:
            continue
        if len(tmp_df) > 1:
            i += 1
        else:
            if isinstance(tmp_df[gene_header], float) and (chr1_hap == chrom):
                tmp_df[gene_header] = qseqid
            else:
                paralog_idx = len(tmp_df[paralogs_gene_columns].dropna())
                paralog_header = paralogs_gene_columns[paralog_idx]
                tmp_df[paralog_header] = qseqid
        
        df.loc[tmp_df.index] = tmp_df
        df['chrom'] = chrom
    logging.debug('duplicated genes: {}'.format(i))
    return df

def main(args):
    """
    addBlastGene2FinalTable <final.table> <blast.out> [Options]
        Add unmap genes into allele table
    """
    p = argparse.ArgumentParser(prog=main.__name__,
                        description=main.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('table', 
            help='final allel table')
    pReq.add_argument('blast', help='blast results')
    pReq.add_argument('gff', help='gff file from MCScanX')
    pOpt.add_argument('--suffix_length', type=int, default=2,
            help='suffix length of chromosome homologus [default: %(default)s]')
    pOpt.add_argument('-t', '--threads', type=int, default=8,
            help='number of program threads[default:%(default)s]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    blast_df = import_blast(args.blast, onlyid=True)
    gff_df = import_gff(args.gff)
    gff_df.set_index('gene', inplace=True)
    blast_df = add_chrom_to_blast_df(blast_df, gff_df)
    blast_df = blast_df.drop_duplicates()
    blast_df = blast_df.set_index('chr')
    final_df = import_allele_table(args.table)
    final_all_df = finalTable2AllFrame(final_df)
    dfGroup = final_all_df.groupby('chrom')
    df = applyParallel(dfGroup, 
                    addBlastGene2FinalTable, 
                    args=(blast_df, args.suffix_length),
                    axis=0, threads=args.threads)
    df_all = AllFrame2finalTable(df)
    df_all = df_all[df_all.count(axis=1) > 2]
    df_all.to_csv(args.output, sep='\t', header=True,
                    index=None, na_rep='.')
    
if __name__ == "__main__":
    main(sys.argv[1:])