#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
%(prog)s in.gff3 rename.list > out.gff3
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import gzip
import re

from collections import OrderedDict


def rename_gff_attribute(ingff, inlist):
    chrom_db = dict(i.strip().split() 
                for i in open(inlist) if i.strip())

    if ingff.endswith('.gz'):
        handle = gzip.open(ingff)
    else:
        handle = open(ingff)

    with handle as fp:
        for line in fp:
            if line.startswith('#'):
                output = line.strip()
            else:
                line_list = line.strip().split()
                attributes = OrderedDict(map(lambda x: x.split("="), line_list[8].split(";")))
               
                for attribute in attributes:
                    ID = attributes[attribute]
                    if ID in chrom_db:
                        attributes[attribute] = chrom_db[ID]
                        flag = 1
                    else:
                        logging.warning('There is not chromosome `{}` in `{}` list, will be discard'.format(ID, inlist))
                        flag = 0
                if flag == 1:      
                    attributes_string = ";".join(map(lambda x: "=".join(x), attributes.items()))
                    line_list[8] = attributes_string
                    output = "\t".join(line_list)
                    
                            
                
            print(output, file=sys.stdout)


if __name__ == "__main__":
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('ingff', help='input file of gff')
    pReq.add_argument('chromlist', help='chromosome list of two columns')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args()

    rename_gff_attribute(args.ingff, args.chromlist)
