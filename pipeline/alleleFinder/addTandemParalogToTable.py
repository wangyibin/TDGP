#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import argparse
import os
import os.path
import sys
import logging
import pandas as pd
import numpy as np
import multiprocessing as mp
from pandarallel import pandarallel

from collections import OrderedDict, Counter
from TDGP.pipeline.alleleFinder.utils import *


def create_empty_stat_table(gene_headers):
    columns = ['total']
    columns += ['gene{}'.format(i) 
               for i in range(1, len(gene_headers) + 1)]
    columns += ['paralogs', 'tandem']
    ds = pd.Series(index=columns)
    return ds
    
def statFinalTableForTetra_singleChrom(chrom, df):
    gene_headers = df.columns.to_list()
    gene_headers.remove('|Paralogs')
    gene_headers.remove('Tandem')
    gene_headers.remove('chrom')
    hap_length = len(gene_headers)
    allele_df = df[gene_headers]
    ds = create_empty_stat_table(gene_headers)
    total = 0
    for i in range(hap_length):
        tmp_allele_df = allele_df[allele_df.count(axis=1) == i+1]
        ds['gene{}'.format(i+1)] = len(tmp_allele_df)
        total += ds['gene{}'.format(i+1)]
   
    
    tandem_num = df['Tandem'].str.split(',').apply(lambda x: len(x) 
                            if not isinstance(x, float) and len(x) > 0 else 0, 1).sum()
    paralog_num = df['|Paralogs'].str.split(',').apply(lambda x: len(x) 
                            if not isinstance(x, float) and len(x) > 0 else 0, 1).sum()
    
    #total = total + paralog_num + tandem_num
    ds.total = total
    ds.paralogs = paralog_num
    ds.tandem = tandem_num
    ds.name = chrom
    return ds

def statFinalTableForTetra(res, threads=4):
    gene_headers = res.columns.to_list()
    gene_headers.remove('|Paralogs')
    gene_headers.remove('Tandem')
    gene_headers.remove('chrom')
    gene_headers = ['gene{}'.format(i+1) for i in range(len(gene_headers))]
   
    stat = applyParallel(res.groupby('chrom'), statFinalTableForTetra_singleChrom, threads).T
    stat.loc['Gene with annotated alleles'] = stat[gene_headers].sum()
    stat.loc['Gene with annotated alleles']['total'] = stat.loc['Gene with annotated alleles'][gene_headers].sum()
    stat.loc['Duplicated genes'] = stat[['paralogs', 'tandem']].sum()

    stat.loc['Duplicated genes']['total'] = stat.loc['Duplicated genes'].sum()
   
    columns = ['Total no. of genes']
    columns += ["No. of genes with {} alleles".format(i) 
                       for i in range(1, len(gene_headers) + 1)]
    columns += ['No. of dispersely duplicated genes', 
                   'No. of tandem duplicated genes']
    stat.columns = columns
    stat = stat.astype(pd.Int64Dtype())
    return stat


def addTandemParalogToTable(args):
    """
    %(prog)s <allele.table> <tandem> <gff> [Options]

    """
    p = p=argparse.ArgumentParser(prog=addTandemParalogToTable.__name__,
                        description=addTandemParalogToTable.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('alleleTable', 
            help='allele table with dup_genes')
    pReq.add_argument('tandem', help='tandem list from MCScanX')
    pReq.add_argument('gff', help='gff file, which is input for MCScanX')
    pOpt.add_argument('--stat', default=None, 
            help='output the stat file of final allele table')
    pOpt.add_argument('-t', '--threads', type=int, default=8,
            help='number of program threads[default:%(default)s]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    allele_table = import_allele_table(args.alleleTable)
    gene_headers = allele_table.columns.to_list()
    gene_headers.remove('dup_gene')
    gff_df = import_gff(args.gff, rmtig=True)
    tandem_df = import_tandem(args.tandem)

    gff_df.set_index('gene', inplace=True)

    allele_table_all = alleleTable2AllFrame(allele_table)
    allele_table_all.dropna(how='all', axis=0, inplace=True)
    pandarallel.initialize(nb_workers=args.threads, verbose=0)
    res = allele_table_all.parallel_apply(formatAlleleTableToFinal, axis=1, 
                args=(gene_headers, tandem_df, gff_df))


    res.sort_values(by=['chrom'], inplace=True)
    res_only_genes = res.drop(['chrom'], axis=1)
    res_genes = alleleTable2GeneList(res_only_genes)
    with open("final.allele.table.gene.list", 'w') as out:
        out.write("\n".join(res_genes))
    res.to_csv(args.output, sep='\t', header=True, index=None, na_rep='.')
    if args.stat:
        stat = statFinalTableForTetra(res, args.threads)
        stat.to_csv(args.stat, sep='\t', header=True, index=True, na_rep='-')
        stat.to_excel(args.stat + ".xls", header=True, index=True, 
                        na_rep="-")


if __name__ == "__main__":
    addTandemParalogToTable(sys.argv[1:])
