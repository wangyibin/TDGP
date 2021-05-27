#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
stat allele table
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from utils import *

def statFinalAlleleTable(args):
    """
    %(prog)s <final.table> [Options]
        stat final allele table
    """
    p = argparse.ArgumentParser(prog=statFinalAlleleTable.__name__,
                        description=statFinalAlleleTable.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('table', 
            help='allele table')
    pOpt.add_argument('-t', '--threads', type=int, default=8,
            help='number of program threads[default:%(default)s]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    allele_df = import_allele_table(args.table)
    res = statFinalTable(allele_df, args.threads)
    res.to_csv(args.output, sep='\t', header=True, 
                index=True, na_rep='-')
    res.to_excel(args.output.name + ".xls", header=True, index=True, 
                        na_rep="-")



if __name__ == "__main__":
    statFinalAlleleTable(sys.argv[1:])