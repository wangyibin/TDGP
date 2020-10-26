#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
pair_run_synteny.py -s Chr01g1 Chr01g2 Chr01g3 [Options]
    run synteny in pairs
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from utils import run_synteny, pairs



def main(args):
    p = p=argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-s', '--samples', nargs="+", 
            help='sample of pairs', required=True)
    pOpt.add_argument('-r', '--rerun', action='store_true', default=False,
            help='rerun all sample pairs [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    

    samples_pairs = pairs(args.samples)

    for pair in samples_pairs:

        if args.rerun:
            os.system('rm -rf {}'.format(pair))
        else:
            if op.exists(pair):
                continue

        run_synteny(pair)
    
if __name__ == "__main__":
    main(sys.argv[1:])