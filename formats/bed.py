#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Time: 2019/7/22 8:38

"""
bed format and some util tools.
"""
from __future__ import print_function

import argparse
import os.path as op
import os
import pandas as pd
import numpy as np
import sys

from TDGP.apps.base import ActionDispatcher
from TDGP.apps.base import check_file_exists, debug
from TDGP.apps.base import listify

def main():

    actions = (
            ("randomBed", "random select several gene from bed file."),
            ("countBedByChr", "count bed file by chromosome"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

class BedpeLine(object):
    def __init__(self):
        self.chrom1 = ''
        self.start1 = 0
        self.end2 = 0
        self.chrom2 = ''
        self.start2 = 0
        self.end2 = 0
        self.name = '.'
        self.score = '.'
        self.strand1 = '.'
        self.strand2 = '.'
        self.other = []

def import_bed6(bed, sort=False):
    """
    import 6-columns bed file
    """
    names = ['chrom', 'start', 'end', 
                'name', 'score', 'strand']
    df = pd.read_csv(bed, sep='\t', header=None,
                    index_col=None, names=names)
    
    if sort:
        df = df.sort_values(by=['chrom', 'start'])
    
    return df


    

## out command
def countBedByChr(args):
    """
    %(prog)s <bed> [Options]
        count bed file by chromosome
    """

    p = argparse.ArgumentParser(prog=countBedByChr.__name__,
                        description=countBedByChr.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bed', 
            help='input bed file')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    df = import_bed6(args.bed)
    df = df.groupby('chrom').size().reset_index(name='counts')

    df.to_csv(args.output, sep='\t', header=None, index=None)

def random_select(args):
    """
    %(prog)s <gene.bed> <--target target.bed/-n nums> [Options]

        random select several genes from bed file.
    """
    p = argparse.ArgumentParser(prog=__file__,
                        description=random_select.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bed', help='source bed file')
    pOpt.add_argument('--target', help='target bed file')
    pOpt.add_argument('-n', '--nums', type=int, 
            help='select numbers (will be overlaped with `--target`)')
    pOpt.add_argument('-o', '--out', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('--include', action='store_true', default=False,
            help='if include the target gene in source bedfile [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    source_df = pd.read_csv(args.bed, sep='\t', header=None, index_col=0,
            names=('chrom', 'start', 'end', 'gene'))
    if not args.target:
        if args.nums:
            nums = args.nums
        else:
            logging.error('Must input `--nums` or `--target` options')
    else:
        target_df = pd.read_csv(args.target, sep='\t', header=None,
                index_col=0, names=('chrom', 'start', 'end', 'gene'))
        nums = len(target_df)

    if args.target and not args.include:
        source_df = source_df[~source_df['gene'].isin(target_df['gene'])]
    
    random_df = source_df.sample(nums)
    random_df = random_df.sort_values(by=['chrom', 'start'])
    random_df.to_csv(args.out, sep='\t', header=None)


if __name__ == "__main__":
    main()