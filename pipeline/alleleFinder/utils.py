#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
util function for alleleFinder
"""
from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import re
import sys
import string
import pandas as pd
import warnings
import multiprocessing as mp
from Bio import SeqIO
from collections import OrderedDict, Counter
from itertools import combinations
from joblib import Parallel, delayed
from pandarallel import pandarallel
from pandas.core.indexing import convert_from_missing_indexer_tuple

warnings.filterwarnings('ignore')

PBS_HEADER = """#!/bin/bash
#PBS -j oe {}
#PBS -q {}
#PBS -V 
#PBS -l nodes=1:ppn={} {}
if [[ ! -z $PBS_O_WORKDIR ]]; then
    cd $PBS_O_WORKDIR
fi
"""

SGE_HEADER = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd {}
#$ -q {}
#$ -pe mpi {} {}
"""


def debug(level=logging.DEBUG):
    """
    Basic config logging format
    """
    from TDGP.apps.font import magenta, green, yellow, white
    formats = white("%(asctime)s") 
    formats += magenta(" <%(module)s:%(funcName)s>")
    formats += white(" [%(levelname)s]")
    formats += yellow(" %(message)s")
    logging.basicConfig(level=level, format=formats, datefmt="[%Y-%m-%d %H:%M:%S]")

debug()

class Cluster(object):
    """
    class of cluster operation
        in order to execute successful should set 
                the `CLUSTER` variants into ENV
    Params:
    --------
    
    Returns:
    ---------
    out: `str`: CLUSTER

    Functions:
    ---------
    get_header: get the header of cluster system with parameters
    get_raw_header: get the raw of cluster system with parameters

    """
    
    def __init__(self, cluster=None, 
                    name=None, queue=None, 
                    threads=1, array=None):
        self.CLUSTER = cluster if cluster else None
        if not self.CLUSTER:
            self.get()
        self.get_header(name, queue, threads, array)
        self.get_raw_header()

    def get(self):
        """
        To obtain the environment of `CLUSTER`,
            if not found will be set default `SGE`.
        """
        try:
            self.CLUSTER = os.environ['CLUSTER']
        except KeyError:
            self.CLUSTER = 'SGE'
            logging.warning('There is not environment `CLUSTER` in PATH')

        return self.CLUSTER


    def get_header(self, name=None, queue=None, 
                        threads=1, array=None):
        """
        According to the environment of `CLUSTER` to 
            return a header of cluster system
        """
        if self.CLUSTER.upper() == "SGE":
            name = "\n#$ -N " + name  if name else ""
            queue = queue if queue else "all.q"
            array = "\n#$ -t " + array if array else ""
            self.header = SGE_HEADER.format(name, queue, threads, array)   
        elif self.CLUSTER.upper() == "PBS":
            name = "\n#PBS -N " + name if name else ""
            queue = queue if queue else "workq"
            array = "\n#PBS -J " + array if array else ""
            self.header = PBS_HEADER.format(name, queue, threads, array)
        elif self.CLUSTER.upper() == "TORQUE":
            name = "\nPBS -N " + name if name else ""
            queue = queue if queue else "share"
            array = "\n#PBS -J " + array if array else ""
            self.header = PBS_HEADER.format(name, queue, threads, array)
        else:
            logging.warning("there is not of header "
                            "of cluster:`{}`".format(self.CLUSTER))
            sys.exit()
        return self.header

    def get_raw_header(self):
        """
        According to the environment of `CLUSTER` to 
            return a header of cluster system
        """
        if self.CLUSTER.upper() == "SGE":
            self.raw_header = SGE_HEADER
        elif self.CLUSTER.upper() == "PBS":
            self.raw_header = PBS_HEADER
        elif self.CLUSTER.upper() == "TORQUE":
            self.raw_header = PBS_HEADER
        else:
            logging.warning("there is not of header "
                            "of cluster:`{}`".format(self.CLUSTER))
            sys.exit()
        return self.raw_header


    def __str__(self):
        return self.CLUSTER

    __retr__ = __str__

class BlockHeader(object):
    def __init__(self, line):
        self.header = line.strip()
        self.list = self.header.split(" ")
        self.num = int(self.list[2].strip(":"))
        self.score = float(self.list[3].split("=")[1])
        self.evalue = float(self.list[4].split("=")[1])
        self.N = int(self.list[5].split("=")[1])
        self.pairs = self.list[6]
    def __str__(self):
        return self.header
    __repr__ = __str__


class collinearity(object):
    def __init__(self, data):
        self.hap_suffix=['g1', 'g2', 'g3', 'g4']
        self.length = len(self.hap_suffix[0])
        self.gene_headers = [ 'gene' + s 
                             for s in string.ascii_uppercase[:len(self.hap_suffix) ]]
        self.data = data
        if not op.exists(self.data):
            logging.error('No such file of `{}`'.format(self.data))
            sys.exit()
        self.load()

    def load(self):
        results = []
        with open(self.data, 'r') as fp:
            block_header = ""
            for line in fp:
                if line.startswith("#"):
                    if line[:12] == "## Alignment":
                        block_header = BlockHeader(line)
                    else:
                        continue
                else:
                    if block_header:
                        block_pairs = block_header.pairs.split("&")
                        line_list = line.strip().split()
                        line_list = re.split("\s+|-", line.strip())
                        
                        line_list = block_pairs + list(filter( lambda x: x!="", line_list))[:4]
                        results.append(line_list)
        self.df = pd.DataFrame(results, columns=['chr1', 'chr2', 'blockN', 'geneN', 
                            'gene1', 'gene2'])
        return self.df
    
    @property
    def inter_df(self):
        return self.df.loc[self.df.chr1 != self.df.chr2].reset_index()
    @property
    def inter_only_homo_df(self):
        return self.inter_df.loc[self.inter_df.chr1.str[:-self.length] ==
                            self.inter_df.chr2.str[:-self.length]].reset_index()   
    @property
    def intra_df(self):
        return self.df.loc[self.df.chr1 == self.df.chr2].reset_index()
    
    def pairs(self, samples):
        samples = list(combinations(samples, 2))

        return list(map(lambda x:"{}-{}".format(x[0], x[1]), samples))

    def to_anchor_per_hap(self, only_homo=True, with_chr=False, 
                          store=True):
        length = self.length
        anchor_res_db = OrderedDict()
        if only_homo:
            df = self.inter_only_homo_df
        else:
            df = self.inter_df
       
        for pair in self.pairs(self.hap_suffix):
            pair_list = pair.split("-")
            tmp_df = df.loc[(df.chr1.str[-length:] == pair_list[0]) &
                            (df.chr2.str[-length:] == pair_list[1])]
            if with_chr:
                anchor_res_db[pair] = tmp_df[['chr1', 'chr2', 'gene1', 'gene2']]
            else:
                anchor_res_db[pair] = tmp_df[['gene1', 'gene2']]
            
            if store:
                anchor_res_db[pair].to_csv("{}.anchors".format(pair), 
                                  sep='\t', index=None, header=None)
        if store:
            df = pd.concat(list(anchor_res_db.values()))
            df.to_csv("{}.all.anchors".format(".".join(self.hap_suffix)),
                     sep='\t', index=None, header=None)
        
        return anchor_res_db


    def createEmptyAlleleTable(self):
        """
        create a empty dataframe for allele table

        Params:
        --------
        chrom_list: `list` or `array-like` list for homologs  tags
                    e.g. ['Chr01g1', 'Chr01g2', 'Chr01g3', 'Chr01g4']
        """

        i = len(self.hap_suffix)
        columns = []
        for s in string.ascii_uppercase[:i]:
            columns.append('gene' + s)

        df = pd.DataFrame(columns=columns)
        return df
    
    def formatAlleleTable(self):
        anchor_res_db = self.to_anchor_per_hap(with_chr=True, store=False)
        for pair in anchor_res_db:
            tmp_df = anchor_res_db[pair]
            hap1, hap2 = pair.split("-")
            pass 
#             tmp_df.groupby()
            
        
    
        

def pairs(samples):
    samples = list(combinations(samples, 2))

    return list(map(lambda x:"{}-{}".format(x[0], x[1]), samples))

def run_synteny(pairs):
    sample1, sample2 = pairs.split("-")
    os.mkdir(pairs)
    os.chdir(pairs)
    header = Cluster().header
    os.system('ln -s ../data/{}* .'.format(sample1))
    os.system('ln -s ../data/{}* .'.format(sample2))
    script = 'run_{}.sh'.format(pairs)
    with open(script, 'w') as out:
        cmd = header + "\n"
        cmd += 'python -m jcvi.compara.catalog ortholog {} {} --cscore=.99 --no_strip_names --nochpf \n'.format(sample1, sample2)
        out.write(cmd)
    os.system('qsub {}'.format(script))
    os.chdir('../')


def split_bed_by_chr(bedfile, outdir='./', prefix=None):
    """
    split bed by chromosome name
    """
    # if prefix is None:
    #     prefix = op.splitext(op.basename(bedfile))[0]
           
    df = pd.read_csv(bedfile, sep='\t', header=None, index_col=None)
    for chrom, data in df.groupby(0):
        out = '{}/{}.bed'.format(outdir, chrom)
        data.to_csv(out, sep='\t',
                    header=None, index=None)
        logging.debug('output new bed in `{}`'.format(out))

def which_hap(gene, sep="|"):
    if "|" in gene:
        return gene.split(sep, 1)[0]
    else:
        return



def extract_fasta(infasta, gene_set, output_handle, exclude=False):
    """
    extract sequences by list
    
    Params:
    --------
    infasta: `str` 
            fasta file
    gene_set: `set` or `list-like` 
            gene lists
    output_handle: `handle` of output
    exclude: `bool` if True sequence will be excluded [default: False]
    Returns:
    --------
    write extracted fasta to a fasta file

    Examples:
    --------
    >>> extract_fasta("sample.fasta", gene_set, output_handle)
    """
    if infasta[-2:] == "gz":
        fp = gzip.open(infasta)
    else:
        fp = open(infasta)
    
    fa = SeqIO.parse(fp, 'fasta')
    for record in fa:
        if exclude:
            if record.id not in gene_set:
                SeqIO.write(record, output_handle, 'fasta')
        else:
            if record.id in gene_set:
                SeqIO.write(record, output_handle, 'fasta')
        



def rename_gff_by_strings_per_hap(ingff, output_handle, hap_suffix_list=['g1', 'g2', 'g3', 'g4'],
                         gene_idx=1,  prefix="", suffix=""):
    """
    rename MCScanX gff ID by a suffix/prefix string

    Params:
    --------
    ingff: `str` input gff file 
    output_handle:  `handle` of output
    gene_idx: `int` gene column number
    prefix/suffix: `str` strings of gene ID prefix/suffix [default: ""]

    Returns:
    ---------
    write renamed gff to a new file

    Examples:
    --------
    >>> rename_gff_by_strings("sample.gff", output_handle, prefix="tetra")
    """
    hap_suffix_len = len(hap_suffix_list[0])
    gene_prefix_list = []
    for i in range(len(hap_suffix_list)):
        gene_prefix_list.append(string.ascii_uppercase[i])
    gene_prefix_db = dict(zip(hap_suffix_list, gene_prefix_list))

    gff_df = pd.read_csv(ingff, sep='\t', header=None, 
                index_col=None, comment='#')
    groupby_iter = gff_df.groupby(gff_df[0].str[-hap_suffix_len:])
    results = []
    for group, df in groupby_iter:
        if group in hap_suffix_list:
            prefix = gene_prefix_db[group]
            df.loc[:, gene_idx] = "{}|".format(prefix) + df[gene_idx]
        results.append(df)
    gff_df = pd.concat(results, axis=0)
    gff_df = gff_df.sort_values(by=[0, 2])
    gff_df.to_csv(output_handle, sep='\t', header=None, index=None)


def rename_gff_by_strings(ingff, output_handle, gene_idx=1,  prefix="", suffix=""):
    """
    rename MCScanX gff ID by a suffix/prefix string

    Params:
    --------
    ingff: `str` input gff file 
    output_handle:  `handle` of output
    gene_idx: `int` gene column number
    prefix/suffix: `str` strings of gene ID prefix/suffix [default: ""]

    Returns:
    ---------
    write renamed gff to a new file

    Examples:
    --------
    >>> rename_gff_by_strings("sample.gff", output_handle, prefix="tetra")
    """
    gff_df = pd.read_csv(ingff, sep='\t', header=None, 
                index_col=None, comment='#')
    gff_df[gene_idx] = prefix + gff_df[gene_idx] + suffix
    gff_df.to_csv(output_handle, sep='\t', header=None, index=None)


def rename_fasta_by_strings(infasta, output_handle, prefix="", suffix=""):
    """
    rename fasta ID by a suffix/prefix string

    Params:
    --------
    ingff: `str` input fasta file 
    output_handle:  `handle` of output
    prefix/suffix: `str` strings of gene ID prefix/suffix [default: ""]

    Returns:
    ---------
    write renamed fasta to a new file

    Examples:
    --------
    >>> rename_fasta_by_strings("sample.fasta", output_handle, prefix="tetra")
    """
    if infasta[-2:] == "gz":
        fp = gzip.open(infasta)
    else:
        fp = open(infasta)
    
    fa = SeqIO.parse(fp, 'fasta')
    for record in fa:
        renamed_id = "{}{}{}".format(prefix, record.id, suffix)
        record.id = renamed_id
        record.description = ""
        SeqIO.write(record, output_handle, 'fasta')
        

def rename_fasta(infasta, rename_db, output_handle):
    """
    rename fasta id by a rename database

    Params:
    --------
    infasta: `str` input fasta file
    rename_db: `dict` rename database
    output_handle:  `handle` of output

    Returns:
    ---------
    write renamed fasta to a fasta file

    Examples:
    --------
    >>> rename_fasta("sample.fasta", rename_db, output_handle)
    """
    if infasta[-2:] == "gz":
        fp = gzip.open(infasta)
    else:
        fp = open(infasta)
    
    fa = SeqIO.parse(fp, 'fasta')
    for record in fa:
        if record.id in rename_db:
            renamed_id = rename_db[record.id]
            record.id = renamed_id
            record.description = ""
            SeqIO.write(record, output_handle, 'fasta')
        else:
            continue


def create_empty_allele_table(chrom_list):
    """
    create a empty dataframe for allele table

    Params:
    --------
    homoTags: `list` or `array-like` list for homologs  tags
                e.g. ['Chr01g1', 'Chr01g2', 'Chr01g3', 'Chr01g4']
    """

    i = len(chrom_list)
    columns = []
    for s in string.ascii_uppercase[:i]:
        columns.append('gene' + s)
    
    df = pd.DataFrame(columns=columns)

    return df

def import_anchor(anchor):
    """
    import anchor file
    """

    df = pd.read_csv(anchor, sep='\t', header=None, 
            index_col=None, comment="#", usecols=[0, 1])
    
    return df


def import_blast(blast, noself=True, onlyid=False):
    """
    import blast file as dataframe
    """
    header = ['qseqid', 'sseqid', 'identity', 'length', 
                'mismatch', 'gap', 'qstart', 'qend', 
                'sstart', 'send', 'evalue', 'bitscore']
    
    if not op.exists(blast):
        logging.error('No such file of `{}`.'.format(blast))
        sys.exit()
    else:
        logging.debug('Load file of `{}`.'.format(blast))
    
    if onlyid:
        header = header[:2]
        logging.debug('Only load two column of qseqid and sseqid.')
    usecols = list(range(len(header)))
    df = pd.read_csv(blast, sep='\t', header=None,
                index_col=None, names=header, usecols=usecols)
    
    if noself:
        df = df[df.qseqid != df.sseqid]
        logging.debug('Only load with self pairs.')
    
    return df


def alleleTable2GeneList(alleleTable_df):
    """
    convert allele table dataframe to gene list
    """
    genes = []
    for column in alleleTable_df.columns:
        l = alleleTable_df[column].dropna().to_list()
        genes += l
    return genes
    
def alleleTable2AllFrame(alleleTable_df):
    """
    convert allele table dataframe to a dataframe with all dup gene
    
    """
    dup_gene_df = alleleTable_df.dup_gene.str.split(",").apply(pd.Series, 1)
    dup_gene_df.columns = list(map(lambda x: "dup_gene{}".format(x), dup_gene_df.columns))
    add_dup_df2 = pd.concat([ dup_gene_df, alleleTable_df], axis=1)
    add_dup_df2.drop('dup_gene', axis=1, inplace=True)
    
    return add_dup_df2


def AllFrame2alleleTable(alleleTable_df):
    """
    convert all allele table dataframe to a dataframe with one dup gene
    
    """
    dup_gene_columns = alleleTable_df.columns[
        alleleTable_df.columns.str.find('dup_gene').map(lambda x: x == 0)]
    def func(x):
        res = x[dup_gene_columns].dropna()
        if not res.empty:
            return ",".join(res)
        else:
            return np.nan
    add_dup_df2 = alleleTable_df.copy()
    add_dup_df2.drop(dup_gene_columns, axis=1, inplace=True)
    add_dup_df2['dup_gene'] = alleleTable_df.apply(func, 1).dropna()
    
    return add_dup_df2

def remove_dup_in_dup_gene_columns(add_dup_df, dup_genes, 
                            gene_headers=['geneA', 'geneB', 'geneC', 'geneD']):
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

            allele_genes = add_dup_df.loc[idxmax, gene_headers].to_list()
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


def applyParallel(dfGrouped, func, threads=4):
    """
    parallel apply a func for pandas groupby 
    ![https://stackoverflow.com/questions/26187759/parallelize-apply-after-pandas-groupby]
    """
    results = Parallel(n_jobs=threads)(delayed(func)(name, group) 
                                       for name, group in dfGrouped)
    return pd.concat(results)



