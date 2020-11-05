#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
generate allele table from gmap results
"""
from __future__ import print_function



import argparse
import logging
import os
import os.path as op
import sys
import pandas as pd
import numpy as np


from collections import Counter
from pandarallel import pandarallel

from utils import applyParallel, import_blast
from utils import alleleTable2AllFrame, alleleTable2GeneList
from utils import AllFrame2alleleTable, which_hap
from utils import remove_dup_from_allele_table, get_dup_genes
from utils import formatAlleleTable


def getAllGeneFromGmap(gene, df):
    empty_df = pd.DataFrame(columns=['all_gene'], index=[gene])
    query_genes = df[1].to_list()
    if len(query_genes) < 2:
         return empty_df
    else:
        empty_df.loc[gene].all_gene = query_genes
    
    return empty_df


def filter_by_blast(row, blast_df):
    genes = set(row.all_gene)
    try:
        blast_genes = blast_df.loc[genes].loc[genes].index.to_list()
    except KeyError:
        return []
    hit_gene_list = set()
    for pair in blast_genes:
        hit_gene_list.add(pair[0])
        hit_gene_list.add(pair[1])
    hit_gene_list = set(hit_gene_list)
    hap = set(map(which_hap, hit_gene_list))
    if len(hap) < 2:
        return []
    hit_gene_list = list(hit_gene_list)
    
    return hit_gene_list

def formatAlleleTableFromGmap(row):
    query_genes = row.all_gene
    for query_gene in sorted(query_genes):
        gene_header = "gene{}".format(which_hap(query_gene))
        if gene_header == 'geneNone':
            continue
        if isinstance(row[gene_header], float):
            row[gene_header] = query_gene
        else:
            row['dup_gene'].append(query_gene)
    return row.drop('all_gene')


def alleleTableFromGmap(args):
    """
    %(prog)s <poly2mono_gene.list> [Options]

        generate allele table from gmap results
    """
    p = p=argparse.ArgumentParser(prog=alleleTableFromGmap.__name__,
                        description=alleleTableFromGmap.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('poly2mono', 
            help='poly2mono list (two columns)')
    pReq.add_argument('blast', help='blast result after filter')
    pOpt.add_argument('--gene_headers', nargs="+", 
            default=['geneA', 'geneB', 'geneC', 'geneD'],
            help='gene headers of table [default: %(default)s]')
    pOpt.add_argument('-t', '--threads', type=int, default=4, 
            help='number of threads [defualt: %(default)s]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    ## import blast return only two id columns and set as index
    blast_df = import_blast(args.blast, onlyid=True)
    blast_df.set_index(['qseqid', 'sseqid'], inplace=True)

    ## import gmap list

    if not op.exists:
        logging.error("No such file of `{}`".format(args.poly2mono))
        sys.exit()
    else:
        logging.debug("Load file of `{}`".format(args.poly2mono))
    
    logging.debug("get all gene from gmap ...")
    gmap_gene_df = pd.read_csv(args.poly2mono, sep='\t', 
                    header=None, index_col=0)
    group_df = gmap_gene_df.groupby(0)
    all_gene_df = applyParallel(group_df, getAllGeneFromGmap)
    all_gene_df = all_gene_df.dropna()

    ## filter gmap result by blast results
    logging.debug("filter all gene by blast result ...")
    pandarallel.initialize(nb_workers=args.threads, verbose=0)
    all_gene_df['all_gene'] = all_gene_df.parallel_apply(filter_by_blast, 1, args=(blast_df,))
    all_gene_df['all_gene'] = all_gene_df.all_gene.apply(lambda d: d if d else np.nan)
    all_gene_df.dropna(inplace=True)

    ## format allele table 

    logging.debug("format allele table ...")
    index_genes = sorted(set(all_gene_df.index))
    allele_df = pd.DataFrame(columns=args.gene_headers + ['dup_gene'], index=index_genes)
    allele_df['all_gene'] = all_gene_df['all_gene']
    allele_df['dup_gene'] = [[] for _ in range(len(index_genes))]

    res_df = allele_df.parallel_apply(formatAlleleTableFromGmap, axis=1)
    res_df.dup_gene = res_df.dup_gene.apply(lambda x: ",".join(sorted(x)) if x else np.nan)
    res_df = res_df.drop_duplicates()
    res_all_df = alleleTable2AllFrame(res_df)
    
    ## remove dup genes
    logging.debug("remove duplicated gene from allele table ...")

    rmdup_all_df = remove_dup_from_allele_table(res_all_df)
    final_rmdup_all_df = rmdup_all_df.apply(formatAlleleTable, 1, args=(args.gene_headers, ))
    final_rmdup_all_df = final_rmdup_all_df.dropna(axis=1, how='all')
    final_rmdup_df = AllFrame2alleleTable(final_rmdup_all_df)
    
    final_rmdup_all_gene_list = alleleTable2GeneList(final_rmdup_all_df)

    final_rmdup_df.to_csv(args.output, sep='\t', 
                header=True, index=None, na_rep='.')
    logging.debug("Successful, output allele table in `{}`".format(args.output.name))
    with open('alleleTableFromGmap_gene.list', 'w') as out:
        print("\n".join(sorted(set(final_rmdup_all_gene_list))), file=out)
    

if __name__ == "__main__":
    alleleTableFromGmap(sys.argv[1:]) 