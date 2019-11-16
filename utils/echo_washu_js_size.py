#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import sys

from collections import OrderedDict

def echo_washu(size_file):
    
    db = OrderedDict([i.strip().split() 
        for i in open(size_file) if i.strip()])
    for chrom in db.keys():
        print('    new Chromosome("{}", {}),'.format(chrom, db[chrom]))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: {} chrom.sizes".format(sys.argv[0]))
        sys.exit()
    echo_washu(sys.argv[1])
