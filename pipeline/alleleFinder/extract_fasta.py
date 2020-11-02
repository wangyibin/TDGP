#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
extract_fasta.py <in.fasta> <list> [Options]
    extract some sequences from fasta by list
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys
import gzip
from utils import extract_fasta

def main(args):
    p = p=argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('infasta', 
            help='input fasta')
    pReq.add_argument('list', help='list of sequences id')

    pOpt.add_argument('-e', '--exclude', action='store_true', 
            default=False, help='exclude these sequences [default: %(default)s]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    id_list = set([i.strip() for i in open(args.list)])
    extract_fasta(args.infasta, id_list, 
                    args.output, args.exclude)


if __name__ == "__main__":
    main(sys.argv[1:])