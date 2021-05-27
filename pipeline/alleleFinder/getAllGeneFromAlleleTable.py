#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
get all gene from final table to a list file
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from utils import import_allele_table, alleleTable2GeneList
from utils import finalTable2AllFrame, alleleTable2AllFrame


def getAllGeneFromAlleleTable(args):
    """
    %(prog)s <final.table> [Options]
    """
    
    p = argparse.ArgumentParser(prog=getAllGeneFromAlleleTable.__name__,
                        description=getAllGeneFromAlleleTable.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('table', 
            help='allele table')
    pOpt.add_argument('--final', default=False, action='store_true',
            help='if input a final allele table with chromosome [default: %(default)s]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    allele_df = import_allele_table(args.table)
    if args.final:
        allele_df.drop(['chrom'], axis=1, inplace=True)
        allele_all_df = finalTable2AllFrame(allele_df)
    else:
        allele_all_df = alleleTable2AllFrame(allele_df)
    geneList = alleleTable2GeneList(allele_all_df)
    print("\n".join(geneList), file=args.output)


if __name__ == "__main__":
    getAllGeneFromAlleleTable(sys.argv[1:])