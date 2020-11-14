#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
generate adjust table
"""
from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from collections import OrderedDict

def generate_adjust_list(args):
    """
    %(prog)s <haplotype.list> <adjust.tsv> [Options]

        generate an adjust list for each chromosome
    """
    p = p=argparse.ArgumentParser(prog=generate_adjust_list.__name__,
                        description=generate_adjust_list.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('table', 
            help='table of haplotype, two columns: '
            '                   example: Chr01    Chr01g1'
            '                            Chr01    Chr01g2'
            '                            Chr01    Chr02g1'
            '                            ...') 
    pReq.add_argument('list', help='list of haplotype and best assembly')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    haplotype_db = OrderedDict()
    with open(args.table) as  fp:
        for line in fp:
            hap, chrom = line.strip().split()
            if hap not in haplotype_db:
                haplotype_db[hap] = []
            haplotype_db[hap].append(chrom)
    
    reference_db = OrderedDict(i.strip().split() 
                            for i in open(args.list)
                            if i.strip())

    for hap in haplotype_db:
        for chrom in haplotype_db[hap]:
            reference = reference_db[hap]
            
            print("\t".join((chrom, reference)), file=args.output)    


if __name__ == "__main__":
    generate_adjust_list(sys.argv[1:])