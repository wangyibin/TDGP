#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
prune contact between allelic contigs
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd

from collections import OrderedDict
from itertools import combinations
from rich.logging import Console, RichHandler

logging.basicConfig(
    level=logging.INFO,
    format="%(message)s",
    datefmt="[%X]",
    handlers=[RichHandler(console=Console(stderr=True))]
)


logger = logging.getLogger(__name__)
class AlleleLine:
    def __init__(self, line):
        self.line = line.strip()
        line_list = line.strip().split()
        self.chrom = line_list[0]
        self.gene = line_list[1]
        self.contigs = line_list[2:]
        self.n = len(self.contigs)
    
    def __str__(self):
        return self.line 
   
    
class AlleleTable:
    """
    Allele table
    """
    def __init__(self, infile):
        self.filename = infile
        if not op.exists(self.filename):
            logger.error(f'No such file of `{self.filename}`.')
            sys.exit()

        logger.info(f'Loading AlleleTable: `{self.filename}`.')
        
    @property
    def data(self):
        _data = OrderedDict()
        with open(self.filename) as fp:
            for line in fp:
                line = AlleleLine(line)
                if line.chrom not in _data:
                    _data[line.chrom] = []
                _data[line.chrom].append(line.contigs)
        
        return _data

    @property
    def groups(self):
        _group = OrderedDict()
        with open(self.filename) as fp:
            for line in fp:
                line = AlleleLine(line)
                if line.chrom not in _group:
                    _group[line.chrom] = []
                _group[line.chrom].extend(line.contigs)
        
        for chrom in _group:
            _group[chrom] = set(_group[chrom])
        
        return _group
    
    @property
    def chromnames(self):
        return list(self.groups.keys())

def get_header(infile):
    with open(infile) as fp:
        for line in fp:
            header = line.strip().split()
            break
    return header   


def import_pairs(pairs_table):

    logger.info('Loading `{}`'.format(pairs_table))
    df = pd.read_csv(pairs_table, header=0, index_col=None,
                    sep='\t')
    df = df.astype(
        {'Contig1': 'category', 'Contig2': 'category',
        '#X': 'int32', 'Y': 'int32', 'RE1': 'int64', 
        'RE2': 'int64'}
    )

    return df

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('alleleTable', 
            help='allele table')
    pReq.add_argument('countRE', 
            help='count_RE.txt from allhic extract')
    pReq.add_argument('pairs', 
            help='pairs file from allhic extract')
    
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    at = AlleleTable(args.alleleTable)
    countRE = args.countRE
    logger.info('Loading `{}`'.format(countRE))
    chroms = [i.strip().split()[0] for i in open(countRE)
                    if i[0] != "#"]
    chrom_index = dict(zip(chroms, range(len(chroms))))
    filter_func = lambda x: x in chroms

    allelic_pairs = []
    for chrom in at.chromnames:
        for alleles in at.data[chrom]:
            if len(alleles) == 1:
                continue
            alleles = list(filter(filter_func, alleles))
            alleles = sorted(alleles, key=lambda x: chrom_index[x])
            tmp_pairs = list(combinations(alleles, 2))
            allelic_pairs.extend(tmp_pairs)
    
    allelic_pairs = list(set(allelic_pairs))

    header = get_header(args.pairs)
    pairs_df = import_pairs(args.pairs)
    pairs_df.set_index(['Contig1', 'Contig2'], inplace=True)
    
    valid_allelic_pairs = [pair for pair in allelic_pairs if pair in pairs_df.index]
   
    pairs_df = pairs_df.drop(valid_allelic_pairs, axis=0)
    pairs_df = pairs_df.reset_index()
    pairs_df = pairs_df[header]

    pairs_df.to_csv(args.output, sep='\t', header=True, index=False)
    logger.info('Done, pruned pairs table is output in `{}`'.format(args.output.name))


if __name__ == "__main__":
    main(sys.argv[1:])