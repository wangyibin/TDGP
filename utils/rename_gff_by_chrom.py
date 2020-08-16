#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import gzip
import sys


def rename_gff(in_list, gfffile):
    chrom_db = dict(i.strip().split() 
            for i in open(in_list) if i.strip())

    if gfffile[-2:] == "gz":
        fp = gzip.open(gfffile)
    else:
        fp = open(gfffile)

    for line in fp:
        if line[0] == "#":
            output = line.strip()
        elif line.strip() == "":
            output = line.strip()
        else:
            line_list = line.strip().split()
            if line_list[0] in chrom_db:
                line_list[0] = chrom_db[line_list[0]]
            else:
                continue
            output = "\t".join(line_list)

        print(output, file=sys.stdout)



if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: {} in.gff3 chrom.list > out.gff3".format(sys.argv[0]))
        sys.exit()
    gfffile, in_list = sys.argv[1:]
    rename_gff(in_list, gfffile)
