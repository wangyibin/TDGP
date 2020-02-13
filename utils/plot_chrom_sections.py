#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
%(prog)s 
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys


import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd 


def chrom_ticks_convert(ticks):
    """
    Convert a list of  chromosome size to suitable unit.
    >>> ticks = [10000, 20000, 30000]
    >>> chrom_ticks_convert(ticks)
    ['10', '20', '30Kbp']
    """
    if ticks[-1]  - ticks[1] <= 1e3:
        labels = ["{:,.0f}".format((x)) 
                  for x in ticks] 
        labels[-1] += " bp"
    elif ticks[-1]  - ticks[1] <= 4e5:
        labels = ["{:,.0f}".format((x / 1e3)) 
                  for x in ticks]
        labels[-1] += 'Kbp'
    else:
        labels = ["{:,.1f}".format((x / 1e6)) 
                  for x in ticks]
        labels[-1] += " Mbp"
    
    return labels


def import_data(path):
    df = pd.read_csv(path, sep='\t', header=None)
    return df

def main(chrom_sizes_path, section_path, out, width=0.4):
    chrom_sizes = import_data(chrom_sizes_path)
    section = import_data(section_path)

    fig, ax = plt.subplots(figsize=(7, 5.2))
    rate = 1
    chroms = chrom_sizes[0][::-1]
    for i, (chrom, size) in enumerate(zip(chrom_sizes[0][::-1], chrom_sizes[1][::-1]), 1):
        ax.broken_barh([[0, size*rate]], (i - width/2, width), facecolor='#dcdcdc')
        tmp_df = section.loc[section[0] == chrom].sort_values(1)
        tmp_df[2] = (tmp_df[2] - tmp_df[1]) * rate
        tmp_df[1] = tmp_df[1] * rate
        ax.broken_barh(list(zip(tmp_df[1], tmp_df[2])),(i - width/2, width),facecolor='#bb4853' )
    ax.tick_params(axis='x', labelsize=12, width=1.5)
    ax.tick_params(left='off', axis='y', labelsize=12, width=0, length=0, pad=0)
    ax.spines['bottom'].set_linewidth(1.5)
    #ax.xaxis.set_visible(False)
    ax.set_xticks(np.linspace(0, max(chrom_sizes[1]), 7))
    ax.set_xticklabels(chrom_ticks_convert(np.linspace(0, max(chrom_sizes[1]), 7)))
    ax.set_yticks(list(range(1, len(chroms) + 1)))
    ax.set_yticklabels(chroms)
    sns.despine(trim=True, left=True)
    plt.savefig(out, dpi=300)



if __name__ == "__main__":
    p = p=argparse.ArgumentParser(prog=__file__,
                        description=__doc__ ,
                        conflict_handler=   'resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('chromsize', help='chromosome sizes file')
    pReq.add_argument('section', help='section file')
    pReq.add_argument('-o', '--out', help='output file', required=True)
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args()

    main(args.chromsize, args.section, args.out)