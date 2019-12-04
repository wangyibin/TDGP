#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import os
import os.path as op
import sys

from collections import defaultdict


def create_bed_dict(bedfile):
    bed_dict = {}
    with open(bedfile) as fp:
        for line in fp:
            line_list = line.strip().split()
            chrom, start, end, gene = line_list[:4]

            bed_dict[gene] = (chrom, start, end)
    
    return bed_dict


def convert_anchors(bed1, bed2, anchors, outfile):
    bed1 = create_bed_dict(bed1)
    bed2 = create_bed_dict(bed2)
    
    out = open(outfile, 'w')
    with open(anchors) as fp:
        for line in fp:
            if line.startswith("#"):
                continue

            line_list = line.strip().split()
            gene1, gene2 = line_list[:2]
            
            if gene1 not in bed1 or gene2 not in bed2:
                continue
            #chrom1, start1, end1 = bed1[gene]
            #chrom2, start2, end2 = bed2[gene]
            out.write('\t'.join(['\t'.join(bed1[gene1]), '\t'.join(bed2[gene2]), gene1, gene2]) + "\n")
    
    out.close()
    os.system("sort -V {0} > .{0}".format(outfile))
    os.rename(".%s"%outfile, outfile)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: {} sample1.bed sample2.bed sample1.sample2.anchros out.bed".format(op.basename(sys.argv[0])))
        sys.exit()
    bed1, bed2, anchor, out = sys.argv[1:]
    convert_anchors(bed1, bed2, anchor, out)

                
        

        
