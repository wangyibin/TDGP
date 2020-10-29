#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
convert MCScanX collinearity to anchor file.
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from utils import collinearity


def collinearity2anchor(args):
    """
    %(prog)s <sample.collinearity> [Options]
        convert collinearity to anchor file

    """
    
    p = argparse.ArgumentParser(prog=collinearity2anchor.__name__,
                        description=collinearity2anchor.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('collinearity', 
            help='collinearity file')
    pOpt.add_argument('--hap_suffix', default=['g1', 'g2', 'g3', 'g4'],
            nargs="*", help='hap suffix [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    cl = collinearity(args.collinearity)
    cl.to_anchor_per_hap(hap_suffix=args.hap_suffix)

if __name__ == "__main__":
    collinearity2anchor(sys.argv[1:])
