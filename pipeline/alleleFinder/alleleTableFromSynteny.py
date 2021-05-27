#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
generate allele table from synteny results
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import string
import pandas as pd 

from utils import (pairs, 
                    create_empty_allele_table,
                    import_anchor,
)


def generate_allele_table(args):
    """
    %(prog)s <--path pathtosynteny> <--chrom_list Chr01g1 Chr01g2 ...> [Options]

        Generate a allele table from jcvi synteny results 

        Examples:
            %(prog)s --path ./ --chrom_list Chr01g1 Chr01g2 Chr01g3 Chr01g4 > synteny.tsv
    """
    p = argparse.ArgumentParser(prog=generate_allele_table.__name__,
                        description=generate_allele_table.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-p', '--path', required=True,
            help='path of synteny results ')
    pReq.add_argument('--chrom_list', nargs='+', required=True,
            help='list of chromosome (e.g. Chr01g1 Chr01g2 Chr01g3 Chr01g4)')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    chrom_list = args.chrom_list
    path = args.path
    gene_list = [ 'gene' + s for s in string.ascii_uppercase[:len(chrom_list)] ]
    gene_pairs = pairs(chrom_list)
    for i, pair in enumerate(gene_pairs):
        anchor_df = import_anchor("{}/{}/{}.{}.anchors".format(path, pair, *pair.split("-")))
        chrom1, chrom2 = pair.split("-")
        idx1, idx2 = chrom_list.index(chrom1), chrom_list.index(chrom2)
        gene1, gene2 = gene_list[idx1], gene_list[idx2]
        annother_genes = gene_list.copy()
        annother_genes.remove(gene1)
        annother_genes.remove(gene2)
        if i == 0:
            previous_df = create_empty_allele_table(chrom_list)
            previous_df[[gene1, gene2]] = anchor_df[[0, 1]]
            continue
    
        current_df = create_empty_allele_table(chrom_list)
        current_df[[gene1, gene2]] = anchor_df[[0, 1]]
        merge_df = previous_df.merge(current_df, right_on=gene1, left_on=gene1, sort=True, how='outer')
        
        merge_df = merge_df.drop([annother_genes[0] + "_y", annother_genes[1] + "_y", gene2 + "_x"], axis=1)
        merge_df.columns = list(map(lambda x: x.replace('_x', '').replace('_y', ''), merge_df.columns))
        merge_df = merge_df[gene_list]
        previous_df = merge_df
    merge_df = merge_df.dropna(thresh=2)
    merge_df = merge_df.drop_duplicates(subset=gene_list, keep='first')
    merge_df.to_csv(args.output, sep='\t', header=True, index=None, na_rep='NA')


if __name__ == "__main__":
    generate_allele_table(sys.argv[1:])