#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
rename_gff.py <in.gff> [Options]
    rename some gff gene ID by database

"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys
import gzip
from utils import rename_gff_by_strings, rename_gff_by_strings_per_hap

def main(args):
    p = p=argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('ingff', 
            help='input gff')
    # pOpt.add_argument('-l', '--list', default=None, 
    #         help='list of sequences id (two columns)')
    pOpt.add_argument('--hap_suffix', nargs='*', default=None,
            help='haplotype suffix, if refined, gff will \n'
            'rename by strings per haplotype [default: %(default)s]')
    pOpt.add_argument('--prefix', default='', 
            help='rename gff by add strings in prefix [default: %(default)s]')
    pOpt.add_argument('--suffix', default='', 
            help='rename gff by add strings in suffix [default: %(default)s]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    
    if args.hap_suffix:
        rename_gff_by_strings_per_hap(args.ingff, args.output, 
                hap_suffix_list=args.hap_suffix)
    elif args.prefix or args.suffix:
        rename_gff_by_strings(args.ingff, args.output, 
                    prefix=args.prefix, suffix=args.suffix)
    else:    
        sys.exit(p.print_help())

if __name__ == "__main__":
    main(sys.argv[1:])