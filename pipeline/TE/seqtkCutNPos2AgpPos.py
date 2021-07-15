#!/usr/bin/env python
# -*- coding:utf-8 -*-


import sys
from collections import OrderedDict

def strPos2Pos(s):
    """
    convert string of "chrom:start-end" to (chrom, start, end)
    """
    chrom, pos = s.split(":")
    start, end = pos.split("-")

    return chrom, int(start), int(end)

def main(contig_list):
    contig_list = [i.strip() for i in open(contig_list) 
                    if i.strip()]
    pos_list = list(map(strPos2Pos, contig_list))
    N = 100
    for i, pos in enumerate(pos_list):
        chrom, start, end = pos
        start = start + N * i -1
        end = end + N * i 

        print(f"{chrom}:{start}-{end}", file=sys.stdout)
    



if __name__ == "__main__":
    if len(sys.argv[1:]) < 1:
        print('Usage: <contig.list>')
        sys.exit()
    contig_list = sys.argv[1]
    main(contig_list)
