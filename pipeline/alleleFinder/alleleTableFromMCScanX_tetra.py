#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
format allele table from MCScanX collinearity
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import numpy as np
import pandas as pd

from collections import OrderedDict, Counter
from utils import collinearity
from utils import alleleTable2AllFrame, alleleTable2AllFrame
from utils import AllFrame2alleleTable


def formatAlleleTable_tetra(cl, anchor_df):
    hap_pairs = cl.pairs(cl.hap_suffix)
    for i, pair in enumerate(hap_pairs):
        hap1, hap2 = pair.split("-")
        tmp_df = anchor_df.loc[(anchor_df.chr1.str[-cl.length:] == hap1) &
                        (anchor_df.chr2.str[-cl.length:] == hap2)]
        idx1, idx2 = cl.hap_suffix.index(hap1), cl.hap_suffix.index(hap2)

        gene1, gene2 = cl.gene_headers[idx1], cl.gene_headers[idx2]
        annother_genes = cl.gene_headers.copy()
        annother_genes.remove(gene1)
        annother_genes.remove(gene2)
        if i == 0:
                previous_df = cl.createEmptyAlleleTable()
                previous_df[[gene1, gene2]] = tmp_df[['gene1', 'gene2']]
                previous_df.drop_duplicates(gene1, keep='first')
                continue
        if i < len(cl.gene_headers) - 1:
            current_df = cl.createEmptyAlleleTable()
            current_df[[gene1, gene2]] = tmp_df[['gene1', 'gene2']]
            current_df.drop_duplicates(gene1, keep='first')
            merge_df = previous_df.merge(current_df, right_on=gene1, left_on=gene1, sort=True, how='outer')

            merge_df = merge_df.drop([annother_genes[0] + "_y", annother_genes[1] + "_y", gene2 + "_x"], axis=1)
            merge_df.columns = list(map(lambda x: x.replace('_x', '').replace('_y', ''), merge_df.columns))
            merge_df = merge_df[cl.gene_headers]
            previous_df = merge_df
            if i == 2:
                
                previous_df1 = previous_df
        
        if i >= 3 and i < 5:
            current_df = cl.createEmptyAlleleTable()
            current_df[[gene1, gene2]] = tmp_df[['gene1', 'gene2']]
            current_df.drop_duplicates(gene1, keep='first')
            merge_df = previous_df1.merge(current_df, right_on=gene1, left_on=gene1, sort=True, how='outer', indicator="merge_type")
    #         right_only = merge_df1[merge_df1.merge_type == 'right_only'][[annother_genes[0] + "_x", gene1, gene2 + "_y", annother_genes[1] + "_x"]]
    #         right_only.columns = list(map(lambda x: x.replace('_x', '').replace('_y', ''), right_only.columns))
            merge_df = merge_df[[annother_genes[0] + "_x", gene1, gene2 + "_y", annother_genes[1] + "_x"]]
            merge_df.columns = list(map(lambda x: x.replace('_x', '').replace('_y', ''), merge_df.columns))
            merge_df = merge_df[cl.gene_headers]
            previous_df1 = merge_df
            if i == 4:
                previous_df1 = previous_df1[previous_df1.geneA.isna()]
                previous_df1 = previous_df1.drop_duplicates(gene1)
                previous_df2 = pd.concat([previous_df, previous_df1], axis=0)
        
        if i == 5:
            current_df = cl.createEmptyAlleleTable()
            current_df[[gene1, gene2]] = tmp_df[['gene1', 'gene2']]
            merge_df = previous_df2.merge(current_df, right_on=gene1, left_on=gene1, sort=True, how='outer', indicator="merge_type")
            right_only = merge_df[merge_df.merge_type == 'right_only'][[annother_genes[0] + "_x", gene1, gene2 + "_y", annother_genes[1] + "_x"]]
            right_only.columns = list(map(lambda x: x.replace('_x', '').replace('_y', ''), right_only.columns))
    #         merge_df = merge_df[[annother_genes[0] + "_x", gene1, gene2 + "_y", annother_genes[1] + "_x"]]
    #         merge_df.columns = list(map(lambda x: x.replace('_x', '').replace('_y', ''), merge_df.columns))
    #         merge_df = merge_df[cl.gene_headers]
            
            
            #previous_df2 = previous_df2[previous_df2.geneA.isna() & previous_df2.geneB.isna()]
            
            previous_df2 = right_only
            
            previous_df2 = previous_df2.drop_duplicates('geneC')
            previous_df = pd.concat([previous_df, previous_df1, previous_df2], axis=0)

    #previous_df = previous_df.drop_duplicates(subset=cl.gene_headers, keep='first')
    previous_df = previous_df.sort_values(by=cl.gene_headers)
    for gene in cl.gene_headers:
        previous_df = previous_df[~previous_df[gene].duplicated(keep='first') | previous_df[gene].isna()]
    previous_df.reset_index(inplace=True)
    previous_df.drop('index', axis=1, inplace=True)

    return previous_df


def get_mapped_genes(df):
    mapped_genes = []
    for column in df.columns:
        l = df[column].dropna().to_list()
        mapped_genes += l
    mapped_genes = set(mapped_genes)
    return mapped_genes

def get_total_genes(df):
    total_genes = []
    for column in ['gene1', 'gene2']:
        l = df[column].dropna().to_list()
        total_genes += l
    total_genes = set(total_genes)
    return total_genes

def get_unmap_genes(total_genes, mapped_genes):
    unmap_genes = total_genes - mapped_genes
    return unmap_genes

def which_hap(gene, sep="|"):
    return gene.split(sep, 1)[0]
    

def add_duplicates(unmap_genes, anchor_df, add_dup_df, iter_round=3):
    """
    Params:
    unmap_genes: `set`, unmap gene set
    anchor_df: `dataframe` all gene anchor dataframe with chromosome
    add_dup_df: `dataframe` allele_table for add duplicate genes
    iter_round: 1 ,2 or 3 if 1 add gene1 from anchor_df into allele table 
                        elif 2 add gene2 to allele table
    
    """
    assert iter_round in [1, 2, 3], "iter_round must in [1, 2]"
    if iter_round == 1:
        gs = [['gene1', 'gene2']]
    elif iter_round == 2:
        gs = [['gene2', 'gene1']]
    else:
        gs = [['gene1', 'gene2'], ['gene2', 'gene1']]
    
    for pair in gs:
        g1, g2 = pair
        unmap_anchor_df = anchor_df[anchor_df[g1].isin(unmap_genes) & 
                                    ~anchor_df[g2].isin(unmap_genes)]

        unmap_anchor_df = unmap_anchor_df.reset_index()
        unmap_anchor_df.drop('index', axis=1, inplace=True)
        for idx, row in unmap_anchor_df.iterrows():
            gene1, gene2 = row[[g2, g1]]
            gene_header = "gene{}".format(which_hap(gene1))
            gene_rows = add_dup_df[add_dup_df[gene_header].isin([gene1])]
            if gene_rows.empty:
                continue
            for j, gene_row in gene_rows.iterrows():
                dup_gene = gene_row.dup_gene
                tmp_gene = row[g1]
                hap_header = "gene{}".format(which_hap(tmp_gene))
#                 if len(gene_row) > 1:
#                     print("gene {} is duplicated for {} in dataframe".format(gene1, tmp_gene))
#                     continue
                if isinstance(gene_row[hap_header], float):
                    add_dup_df.loc[j, hap_header] = tmp_gene
                    unmap_genes.discard(tmp_gene)
                    continue

                if gene_row[hap_header] == tmp_gene:
                    continue
                if not isinstance(dup_gene, float):
                    dup_gene = dup_gene.split(",")

                    if tmp_gene not in dup_gene:
                        dup_gene.append(tmp_gene)
                        add_dup_df.loc[j, 'dup_gene'] = ",".join(dup_gene)

                else:
                    add_dup_df.loc[j, 'dup_gene'] = tmp_gene
                unmap_genes.discard(tmp_gene)
    return add_dup_df


def get_dup_genes_from_add_dup_df(add_dup_df):
    """
    get duplicated genes from added duplicated gene dataframe
    which is to remove different gene pairs which parent is in dup_gene columns
    """
    dup_genes = []
    for l in add_dup_df.dup_gene.str.split(",").dropna():
        dup_genes += l
    count = Counter(dup_genes)
    dup_genes = list(filter(lambda x:count[x] > 1, count))
    return dup_genes


def remove_dup_in_dup_gene_columns(cl, add_dup_df, dup_genes):
    res_df = add_dup_df.copy()
    dup_gene_df = add_dup_df.dup_gene.str.split(",").apply(pd.Series, 1)
    dup_gene_df.columns = list(map(lambda x: "dup_gene{}".format(x), dup_gene_df.columns))
    add_dup_df2 = pd.concat([ dup_gene_df, add_dup_df], sort=False, axis=1)
    add_dup_df2.drop('dup_gene', axis=1, inplace=True)
    rm_idx = []
    for gene in dup_genes:
        tmp_df = add_dup_df2.loc[add_dup_df2.isin([gene]).any(1)]
        if tmp_df.empty:
            continue
        if len(tmp_df) < 2:
            continue
        
        idx_series = tmp_df.apply(lambda x: len(x.dropna()), 1).sort_values()
        idx_series = idx_series.index.to_list()
        for j in  range(0, len(tmp_df) - 1):
            if j == 0:
                idxmax = idx_series[j]
                idxmin = idx_series[j + 1]
            else:
                idxmin = idx_series[j + 1]
            rm_idx.append(idxmin)

            allele_genes = add_dup_df.loc[idxmax, cl.gene_headers].to_list()
            old_dup_genes_series = add_dup_df.loc[idxmax, 'dup_gene']

            if isinstance(old_dup_genes_series, float):
                old_dup_genes = []
            else:
                old_dup_genes = old_dup_genes_series.split(",")
            new_dup_genes = add_dup_df2.loc[idxmin].dropna().to_list()
            new_dup_genes = set(new_dup_genes) - set(allele_genes) | set(old_dup_genes) 
            new_dup_genes = set(new_dup_genes) - set(allele_genes)
            res_df.loc[idxmax, 'dup_gene'] = ",".join(sorted(new_dup_genes))
            add_dup_df2.drop(idxmin, 0, inplace=True)
    res_df.drop(rm_idx, axis=0, inplace=True)
    res_df.reset_index(inplace=True)
    res_df.drop('index', axis=1, inplace=True)
    return res_df
    

def final_rescue_unmap_genes(rmdup_all_df, final_unmap_df):
    """
    rescue unmap genes add into allele table
    """
    dup_gene_columns = rmdup_all_df.columns[
        rmdup_all_df.columns.str.find('dup_gene').map(lambda x: x == 0)]
    for i, row in final_unmap_df.iterrows():
        gene1, gene2 = row.gene1, row.gene2
        gene1_header = "gene{}".format(which_hap(gene1))
        gene2_header = "gene{}".format(which_hap(gene2))
        allele_row  = rmdup_all_df[rmdup_all_df.isin([gene1]).any(1)]
        idx = allele_row.index
        if allele_row.empty:
            allele_gene2_row = rmdup_all_df[rmdup_all_df.isin([gene2]).any(1)]
            dup_gene2_idx = len(allele_gene2_row[dup_gene_columns].dropna())
            dup_gene2_header = dup_gene_columns[dup_gene2_idx]
            if allele_gene2_row.empty:

                new_row = pd.DataFrame(columns=rmdup_all_df.columns)
                new_row.loc[0, gene1_header] = gene1
                new_row.loc[0, gene2_header] = gene2
                rmdup_all_df.append(new_row, ignore_index=True)
            else:
                if isinstance(allele_gene2_row[gene1_header], float):
                    allele_gene2_row.loc[:, gene1_header] = gene1
                else:
                    allele_gene2_row.loc[:, dup_gene2_header] = gene1
                rmdup_all_df.loc[allele_gene2_row.index] = allele_gene2_row
        else:
            if isinstance(allele_row[gene2_header], float):
                allele_row.loc[:, gene2_header] = gene2
            else:
                dup_gene1_idx = len(allele_row[dup_gene_columns].dropna())
                dup_gene1_header = dup_gene_columns[dup_gene1_idx]
                allele_row.loc[:, gene2_header] = gene2

            rmdup_all_df.loc[allele_row.index] = allele_row
    
    return rmdup_all_df

def alleleTableFromMCScanX_tetra(args):
    """
    """
    p = p=argparse.ArgumentParser(prog=alleleTableFromMCScanX_tetra.__name__,
                        description=alleleTableFromMCScanX_tetra.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('collinearity', 
            help='collinearity file from MCScanX')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    cl = collinearity(args.collinearity)
    anchor_df = pd.concat(cl.to_anchor_per_hap(store=False, 
                        with_chr=True).values())

    allele_table_df_init = formatAlleleTable_tetra(cl, anchor_df)
    total_genes = get_total_genes(anchor_df)
    mapped_genes = get_mapped_genes(allele_table_df_init)
    unmap_genes = get_unmap_genes(total_genes, mapped_genes)
    add_dup_df = allele_table_df_init.copy() 
    add_dup_df['dup_gene'] = [np.nan] * len(add_dup_df)
    add_dup_df = add_duplicates(unmap_genes, anchor_df, add_dup_df, iter_round=3)

    ## first round remove duplicated genes
    dup_genes_round1 = get_dup_genes_from_add_dup_df(add_dup_df)
    df = add_dup_df.copy()
    rmdup_df_round1 = remove_dup_in_dup_gene_columns(cl, df, dup_genes_round1)

    all_df_round1 = alleleTable2AllFrame(rmdup_df_round1)

    all_gene_list = alleleTable2GeneList(all_df_round1)
    ## second round remove duplicated genes
    count = Counter(all_gene_list)
    dup_genes_round2 = list(filter(lambda x: count[x] > 1, count))
    rmdup_df_round2 = remove_dup_in_dup_gene_columns(cl, rmdup_df_round1, dup_genes_round2)
    rmdup_all_df_round2 = alleleTable2AllFrame(rmdup_df_round2)
    rmdup_all_gene_list = alleleTable2GeneList(rmdup_all_df_round2)
    
    ## final rescue unmap genes
    final_unmap_genes = total_genes - set(rmdup_all_gene_list)
    final_unmap_df = anchor_df[anchor_df.isin(final_unmap_genes).any(1)]
    rmdup_df_rescue_all = final_rescue_unmap_genes(rmdup_all_df_round2, final_unmap_df)
    rmdup_df_rescue = AllFrame2alleleTable(rmdup_df_rescue_all)
    rmdup_rescue_all_gene_list = alleleTable2GeneList(rmdup_df_rescue_all)



    with open("alleleGeneFromMCScanX.gene.list", 'w') as out:
        out.write("\n".join(rmdup_rescue_all_gene_list))
    
    rmdup_df_rescue.to_csv(args.output, sep='\t', header=True, index=None, na_rep='.')



if __name__ == "__main__":
    alleleTableFromMCScanX_tetra(sys.argv[1:])