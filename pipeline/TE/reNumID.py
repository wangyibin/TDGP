#!/usr/bin/env python
# -*- coding:utf-8 -*-


import sys
import pandas as pd
from collections import OrderedDict

def main(inputfile):
    
    i = 1
    prev_ID = ""
    with open(inputfile, 'r') as fp:
        for i, line in enumerate(fp):
            if not line.strip() or 'score' in line \
                    or 'SW    perc' in line:
                continue
            data = line.split()
            try:
                ID = data[14]
            except IndexError:
                ID = ''
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
                    f" {data[4]:>5} "
                    f" {data[5]:>11} "
                    f" {data[6]:>11} "
                    f"{data[7]:>16} "
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
        print('Usage: <input.out>')
        sys.exit()
    inputfile = sys.argv[1]
    main(inputfile)

