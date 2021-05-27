#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
split_bed_by_chr.py <Chr01.bed> [Options]
split bed by different chromosome
"""
from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys
import pandas as pd
from utils import split_bed_by_chr


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bed', 
            help='bed file')
    pOpt.add_argument('-o', '--outdir', default='./',
            help='output directory [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    split_bed_by_chr(args.bed, outdir=args.outdir)

if __name__ == "__main__":
    main(sys.argv[1:])