#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
convert ko list to correct format for omicshare kegg pathwar enrichment
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys


def dealKOListFromEmapper(args):

    p = p=argparse.ArgumentParser(prog=dealKOListFromEmapper.__name__,
                        description=dealKOListFromEmapper.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('koList', 
            help='kolist from emapper')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    koList = args.koList
    output = args.output
    with open(koList) as fp:
        for line in fp:
            line_list = line.strip().split()
            if len(line_list) <= 1:
                print(line.strip(), file=output)
            else:
                line_list[1] = line_list[1].replace('"', "")
                for i in line_list[1].split(","):
                    print("\t".join([line_list[0], i]), file=output)
             
    

if __name__ == "__main__":
    dealKOListFromEmapper(sys.argv[1:])