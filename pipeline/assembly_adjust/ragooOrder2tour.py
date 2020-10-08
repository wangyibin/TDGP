#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
convert ragoo order 2 tour file
"""
from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from collections import OrderedDict


def import_order(infile):
    db = OrderedDict()
    with open(infile) as fp:
        for line in fp:
            tig, order = line.strip().split()[:2]
            db[tig] = order
    
    return db


def order2tour(order_file, output):
    order_db = import_order(order_file)

    for tig in order_db:
        print(tig + order_db[tig], end=" ", file=output)
    

def ragooOrder2tour(args):
    """
    %(prog)s <ordering.txt> <out.tour> [Options]

        convert ragoo order file to tour file.
    """
    p = p=argparse.ArgumentParser(prog=ragooOrder2tour.__name__,
                        description=ragooOrder2tour.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('order', 
            help='order file from ragoo')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    order2tour(args.order, args.output)


if __name__ == "__main__":
    ragooOrder2tour(sys.argv[1:])