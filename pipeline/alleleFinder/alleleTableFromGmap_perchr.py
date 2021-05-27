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
from joblib import Parallel, delayed
from pandarallel import pandarallel

from utils import import_blast, import_bed4, import_gff
from utils import applyParallel, checkFileExists
from utils import alleleTable2AllFrame, alleleTable2GeneList
from utils import AllFrame2alleleTable
from utils import filterAndFormatForGmap, add_chrom_to_gmap_genes
from utils import add_chrom_to_blast_df
from utils import remove_dup_from_allele_table_single


def alleleTableFromGmap(args):
    """
    %(prog)s <poly2mono_gene.list> [Options]

        generate allele table from gmap results
    """
    p = p = argparse.ArgumentParser(prog=alleleTableFromGmap.__name__,
                                    description=alleleTableFromGmap.__doc__,
                                    formatter_class=argparse.RawTextHelpFormatter,
                                    conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('poly2mono',
                      help='poly2mono list (two columns)')
    pReq.add_argument('bed', help='bed4 file of mono cds')
    pReq.add_argument('blast', help='blast result after filter')
    pReq.add_argument('gff', help='gff file from MCScanX')
    pOpt.add_argument('--gene_headers', nargs="+",
                      default=['geneA', 'geneB', 'geneC', 'geneD'],
                      help='gene headers of table [default: %(default)s]')
    pOpt.add_argument('--withMono', default=False, action='store_true',
                      help="output allele table with mono gene ID [default: %(default)s]")
    pOpt.add_argument('--withChrom', default=False, action='store_true', 
                      help='output allele table with chromosome [default: %(default)s]')
    pOpt.add_argument('-t', '--threads', type=int, default=4,
                      help='number of threads [defualt: %(default)s]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
                      default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
                      help='show help message and exit.')

    args = p.parse_args(args)

 
    list(map(checkFileExists, [args.poly2mono, args.bed,
                        args.blast, args.gff]))

    ## import gene list and add chromosome into dataframe
    logging.debug("get all gene from gmap ...")
    gmap_gene_df = pd.read_csv(args.poly2mono, sep='\t',
                               header=None, index_col=0)
    mono_bed_df = import_bed4(args.bed)
    mono_bed_df.reset_index(inplace=True)
    mono_bed_df.set_index('gene', inplace=True)
    gmap_gene_df = add_chrom_to_gmap_genes(gmap_gene_df, mono_bed_df, 
                                            threads=args.threads)

    ## import blast return only two id columns and set as index
    gff_df = import_gff(args.gff)
    gff_df.set_index('gene', inplace=True)
    blast_df = import_blast(args.blast, onlyid=True)
    blast_df = add_chrom_to_blast_df(blast_df, gff_df, threads=args.threads)
    blast_df.set_index('chr', inplace=True)

    res_df = filterAndFormatForGmap(gmap_gene_df, blast_df, args.threads, 
                                    gene_headers=args.gene_headers)
    res_df.dropna(how='all', subset=args.gene_headers + ['dup_gene'], inplace=True)
    res_df_all = alleleTable2AllFrame(res_df)
    rmdup_df_all = applyParallel(res_df_all.groupby('chr'), 
                    remove_dup_from_allele_table_single, axis=0, 
                    threads=args.threads)
    rmdup_all_gene_list = alleleTable2GeneList(rmdup_df_all)
    rmdup_df = AllFrame2alleleTable(rmdup_df_all)
    withIndex = False
    if not args.withMono:
        rmdup_df.reset_index(drop=True, inplace=True)
    else:
        withIndex = True
    
    if not args.withChrom:
        rmdup_df.drop(['chr'], inplace=True, axis=1)
    else:
        rmdup_df.set_index('chr', inplace=True)
        withIndex = True
    

    rmdup_df.to_csv(args.output, sep='\t', index=withIndex, 
                        header=True, na_rep='.')
    logging.debug(
        "Successful, output allele table in `{}`".format(args.output.name))
    with open('alleleTableFromGmap_gene.list', 'w') as out:
        print("\n".join(sorted(set(rmdup_all_gene_list))), file=out)
    

if __name__ == "__main__":
    alleleTableFromGmap(sys.argv[1:])
