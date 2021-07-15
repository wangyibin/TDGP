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

def main(inputfile, chrom_size):
    chrom_size = dict(i.strip().split() for i in open(chrom_size) 
                            if i.strip())
    i = 1
    prev_ID = ""
    with open(inputfile, 'r') as fp:
        for i, line in enumerate(fp):
            if not line.strip() or 'score' in line \
                    or 'SW    perc' in line:
                continue
            data = line.split()
            old_chrom = data[4]
            old_start, old_end = int(data[5]), int(data[6])
            old_query_left = int(data[7].strip("(").strip(")"))
            chrom, start, end = strPos2Pos(old_chrom)
            new_start = old_start + start - 1
            new_end = old_end + start - 1
            new_query_left = f"({int(chrom_size[chrom]) - new_end})"
            data[4], data[5], data[6], data[7] = chrom, new_start, new_end, new_query_left
            ID = data[14]
            if ID != prev_ID:
                i += 1
            new_ID = prev_ID = i

            try:
                asterisk = data[15]
            except IndexError:
                asterisk = ""
            print(f" {data[0]:>5} "
                    f" {data[1]:>4}"
                    f" {data[2]:>4}"
                    f" {data[3]:>4} "
                    f" {chrom:>5} "
                    f" {new_start:>11} "
                    f" {new_end:>11} "
                    f"{new_query_left:>16} "
                    f"{data[8]} "
                    f"{data[9]:<17} "
                    f" {data[10]:<17} "
                    f" {data[11]:>6} "
                    f" {data[12]:>6} "
                    f" {data[13]:>6}"
                    f" {new_ID:>6}"
                    f" {asterisk} ", file=sys.stdout)
            
                
if __name__ == "__main__":
    if len(sys.argv[1:]) < 1:
        print('Usage: <input.out> <chrom.size>')
        sys.exit()
    inputfile = sys.argv[1]
    chrom_size = sys.argv[2]
    main(inputfile, chrom_size)

