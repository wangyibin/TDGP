#!/usr/bin/env python
# -*- coding:utf-8 -*-


"""
%prog domain.bed [Options] > domain.bedpe

    Convert a bed file 2 a bedpe-like file
    bedfile:chrom    start    end ...
    bedpe-like:chrom1    start1    end1    chrom2    start2    end2
"""


from __future__ import print_function

import gzip

import os.path as op
import sys


def bed2bedpe(bedfile, bedpe, columns="1,2,3"):
    
    columns = list(map(lambda x: int(x) - 1, columns.split(",")))
    if not op.exists(bedfile):

        print("[Error]: No such file of %s"%bedfile, file=sys.stderr)
        sys.exit()

    if bedfile[-3] == ".gz":
        fp = gzip.open(bedfile)
    else:
        fp = open(bedfile)

    for line in fp:
        data = line.strip().split()
        chrom, start, end = [data[i] for i in columns]

        out = '\t'.join([chrom, start, end]*2)
        if bedpe is sys.stdout:
            out_bedpe = sys.stdout
        else:
            out_bedpe = open(bedpe, 'w')
        
        print(out, file=out_bedpe)


if __name__ == "__main__":
    from optparse import OptionParser

    p = OptionParser(__doc__)

    p.add_option('--columns', default='1,2,3',
                help='the columns of chrom,start,end [default: %default]')
    p.add_option('-o', '--outfile', default=sys.stdout,
            help='Output file [default: stdout]')


    opts, args = p.parse_args()

    
    if len(args) != 1:
        sys.exit(p.print_help())

    bedfile,=  args
    bed2bedpe(bedfile, opts.outfile, opts.columns)
