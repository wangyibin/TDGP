#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
rename_fasta.py <in.fasta> [Options]
    rename some sequences from fasta by database

"""

from __future__ import print_function

import argparse
import os
import os.path as op
import sys
from utils import rename_fasta, rename_fasta_by_strings

def main(args):
    p = p=argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('infasta', 
            help='input fasta')
    pOpt.add_argument('-l', '--list', default=None, 
            help='list of sequences id (two columns)')
    pOpt.add_argument('--prefix', default='', 
            help='rename fasta by add strings in prefix [default: %(default)s]')
    pOpt.add_argument('--suffix', default='', 
            help='rename fasta by add strings in suffix [default: %(default)s]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    
    if args.list:
        rename_db = dict(i.strip().split() 
                for i in open(args.list) if i.strip())
        rename_fasta(args.infasta, rename_db, args.output)
    elif args.prefix or args.suffix:
        rename_fasta_by_strings(args.infasta, args.output, 
                    prefix=args.prefix, suffix=args.suffix)
    else:    
        sys.exit(p.print_help())

if __name__ == "__main__":
    main(sys.argv[1:])