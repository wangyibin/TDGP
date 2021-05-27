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
from utils import *


def formatRemainTandemGene_single(row, gene_headers):
    empty_row = create_empty_series_for_final_table(gene_headers)
    chrom = row.chr
    gene1, gene2, hap = row.gene1, row.gene2, row.hap
    if hap == 'geneNone':
        return empty_row
    empty_row['chrom'] = chrom
    empty_row[hap] = gene1
    empty_row['Tandem'] = gene2
    return empty_row

def formatRemainTandemGene(tandem_df, genelist, threads=4,
                           gene_headers=['geneA', 'geneB', 'geneC', 'geneD']):
    tandem_remain_df = tandem_df[~tandem_df.isin(genelist).any(1)]
    tandem_remain_df['hap'] = tandem_remain_df['gene1'].apply(lambda x: "gene{}".format(which_hap(x)))
    pandarallel.initialize(nb_workers=threads, verbose=0)
    df = tandem_remain_df.parallel_apply(formatRemainTandemGene_single, axis=1, args=(gene_headers, ))
    df.dropna(how='all', axis=0, inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df

def addRemainTandemGene(tandem_df, final_df, threads=4):
    gene_headers = get_headers_from_final_table(final_df)
    final_all_df = finalTable2AllFrame(final_df)
    geneList = alleleTable2GeneList(final_all_df)
    tandem_remain_allele_df = formatRemainTandemGene(tandem_df, 
                                geneList, threads, gene_headers)
    

    df = pd.concat([final_df, tandem_remain_allele_df], ignore_index=True)
    df.reset_index(drop=True, inplace=True)
    return df

def remove_dup_from_final_table(final_all_df, iter_max=3, max_dups=10, 
                        gene_headers=['geneA', 'geneB', 'geneC', 'geneD']):
    res_df = final_all_df.copy()
    dup_gene_columns = final_all_df.columns[
        alleleTable_df.columns.str.find('paralogs').map(lambda x: x == 0)]
    tandem_gene_columns = final_all_df.columns[
        alleleTable_df.columns.str.find('tandem').map(lambda x: x == 0)]
    dup_genes = get_dup_genes(final_all_df[gene_headers + tandem_gene_columns])

    iter_num = 1
    while (len(dup_genes) > max_dups) and (iter_num < iter_max):
        add_dup_df2 = res_df.copy()
        rm_idx = []
        for gene in dup_genes:
            tmp_df = add_dup_df2.loc[add_dup_df2.isin([gene]).any(1)]
            if tmp_df.empty:
                continue

            idx_series = tmp_df.apply(lambda x: len(x[gene_headers].dropna()), 1).sort_values(ascending=False)
            idx_series = idx_series.index.to_list()
            
            idxmax = idx_series[0]
            idxmin = idx_series[1:] if len(idx_series) > 1 else [idx_series[0]]

            allele_genes = tmp_df.loc[idxmax, gene_headers].to_list()
            old_dup_genes_series = tmp_df.loc[idxmax].drop(gene_headers).dropna().to_list()

            if isinstance(old_dup_genes_series, float):
                old_dup_genes = []
            else:
                old_dup_genes = old_dup_genes_series

            new_dup_genes = []
            for l in tmp_df.loc[idxmin].values.tolist():
                new_dup_genes += list(filter(lambda x: not isinstance(x, float), l))
                
            new_dup_genes = set(new_dup_genes)  - set(allele_genes) | set(old_dup_genes)
            new_dup_genes = new_dup_genes - set(allele_genes)
            if len(idx_series) > 1:
                rm_idx.extend(idxmin)
                add_dup_df2.drop(idxmin, 0, inplace=True)
        
        new_dup_genes = get_dup_genes(res_df)
        remove_genes_number = len(dup_genes) - len(new_dup_genes)
        dup_genes = new_dup_genes
        res_df.drop(rm_idx, axis=0, inplace=True)
        
        logging.debug("remove {} duplicated genes in iter {}".format(
                                remove_genes_number, iter_num))
        iter_num += 1
        
    return res_df



def addTandemParalogToTable(args):
    """
    %(prog)s <allele.table> <tandem> <gff> [Options]

    """
    p = argparse.ArgumentParser(prog=addTandemParalogToTable.__name__,
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
    tandem_with_chrom_df = add_chrom_to_tandem_df(tandem_df, gff_df, 
                            threads=args.threads)

    allele_table_all = alleleTable2AllFrame(allele_table)
    allele_table_all.dropna(how='all', axis=0, inplace=True)
    
    pandarallel.initialize(nb_workers=args.threads, verbose=0)
    res = allele_table_all.parallel_apply(formatAlleleTableToFinal, axis=1, 
                args=(gene_headers, tandem_df, gff_df))

    res = addRemainTandemGene(tandem_with_chrom_df, res, threads=args.threads)
    

    res.sort_values(by=['chrom'], inplace=True)
    res_only_genes = res.drop(['chrom'], axis=1)
    res_genes = alleleTable2GeneList(finalTable2AllFrame(res_only_genes))
    
    
    with open("final.allele.table.gene.list", 'w') as out:
        out.write("\n".join(res_genes))
    res.to_csv(args.output, sep='\t', header=True, index=None, na_rep='.')
    if args.stat:
        stat = statFinalTable(res, args.threads)
        stat.to_csv(args.stat, sep='\t', header=True, index=True, na_rep='-')
        stat.to_excel(args.stat + ".xls", header=True, index=True, 
                        na_rep="-")


if __name__ == "__main__":
    addTandemParalogToTable(sys.argv[1:])
