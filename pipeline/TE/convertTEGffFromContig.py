#!/usr/bin/env python
# -*- coding:utf-8 -*-


import sys
import pandas as pd
from collections import OrderedDict



def strPos2Pos(s):
    """
    convert string of "chrom:start-end" to (chrom, start, end)
    """
    chrom, pos = s.split(":")
    start, end = pos.split("-")

    return chrom, int(start), int(end)

def convertPos(row):
    res_row = row.copy()
    chrom, start, end = strPos2Pos(row[0])
    res_row[0] = chrom
    res_row[3] = row[3] + start - 1
    res_row[4] = row[4] + start - 1

    return res_row

def main(inputfile):
    df = pd.read_csv(inputfile, sep='\t', header=None, index_col=None, comment="#")
    res_df = df.apply(convertPos, axis=1)

    res_df.to_csv(sys.stdout, sep='\t', header=False, index=False)
if __name__ == "__main__":
    if len(sys.argv[1:]) < 1:
        print('Usage: <input.gff>')
        sys.exit()
    inputfile = sys.argv[1]
    main(inputfile)

