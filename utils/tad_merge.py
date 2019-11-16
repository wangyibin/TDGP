#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import os
import os.path as op
import sys

from collections import OrderedDict
from intervaltree import IntervalTree, Interval


def tad_merge(infile):
    tree_dict = OrderedDict()
    with open(infile) as fp:
        for line in fp:
            chrom, start, end, rank = line.strip().split()
            start, end = int(start), int(end)
            if chrom not in tree_dict:
                tree_dict[chrom] = IntervalTree()
                tree_dict[chrom].addi(start, end, rank)
    
            overlap = list(tree_dict[chrom].overlap(start, end))
            if overlap:
                if (end - start) < overlap[0].length(): 
                    tree_dict[chrom].remove(overlap[0])
                    tree_dict[chrom].addi(start, end, rank)
            else:
                tree_dict[chrom].addi(start, end, rank)


    for chrom in tree_dict:
        for item in sorted(tree_dict[chrom]):
            print("\t".join((chrom, str(item.begin), str(item.end), item.data)), file=sys.stdout)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("*"*79)
        print("Usage: {} hitad.txt".format(op.basename(sys.argv[0])))
        print("*"*79)

        sys.exit()
    tad_merge(sys.argv[1])
                
            
