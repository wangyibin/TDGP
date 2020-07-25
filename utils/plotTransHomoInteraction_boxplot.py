#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
%(prog)s <sample_1000000_iced.cool> [Options]
 
    plotting boxplot of trans interaction within homologs per chromosome.
"""


from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import cooler
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import re
import seaborn as sns


from itertools import combinations, permutations

from TDGP.formats.hicmatrix import cool2matrix


def find_homo_chrom(homo_chroms, all_chroms):
    chrom_pairs = []
    for chrom in homo_chroms:
        regex = re.compile("({})?".format(chrom))
    
        chrom_pairs.append(list(filter(lambda x: regex.match(x).group(), all_chroms)))
    return chrom_pairs

def is_heter(chrom_pairs):
    chrom1, chrom2 = chrom_pairs
    if chrom1[:-1] == chrom2[:-1]:
        return False
    else:
        return True

def chrom_size_convert(size):
    """
    Convert the unit of chromosome size to suitable unit.
    >>> chrom_size_convert(100000)
    100 Kbp
    >>> chrom_size_convert(1000000)
    1 Mbp
    """
    if size <= 1e3:
        label = "{:,.0f}".format((size)) + " bp"
    elif size < 1e6:
        label = "{:,.0f}".format((size / 1e3)) + " Kbp"
    else:
        label = "{:,.0f}".format((size / 1e6)) + " Mbp"
    
    return label

def getChromPairsMatrix(cool, hm, chrom1, chrom2):
    idx1 = cool.bins().fetch(chrom1).index
    idx2 = cool.bins().fetch(chrom2).index
    start1 = idx1[0]
    end1 = idx1[-1]
    start2 = idx2[0]
    end2 = idx2[-1]
    return hm[start1: end1 + 1, start2: end2 + 1]


def main(args):
    p = p=argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-m', '--matrix', required=True,
            help='contact matrix with cool formats')
    pOpt.add_argument('-o', '--output', 
            help='output of picture' )
    pOpt.add_argument('--method', choices=['average', 'median'], 
            default='median',
            help='method of interaction frequency per bin [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    output = args.matrix.rsplit('.', 1)[0] + '_TransInteractionWithinHomo_boxplot.pdf' \
                if not args.output else args.output
    method_func = np.median if args.method == 'median' else np.mean
    cool = cooler.Cooler(args.matrix)
    resolution_string = chrom_size_convert(cool.binsize)
    all_chroms = cool.chromnames
    chrom_heter_pairs = list(filter(is_heter, list(combinations(all_chroms, 2))))
    homo_chroms = sorted(set(map(lambda x: x[:-1], cool.chromnames)))
    
    chrom_homo_pairs = find_homo_chrom(homo_chroms, all_chroms)

    hm = cool2matrix(cool)

    boxprops = dict(color='black', linewidth=1.5)
    medianprops=dict(color='black', linewidth=2.5)
    whiskerprops = dict(linestyle='--')
    fig, axes = plt.subplots(2, 4, figsize=(8, 4.5),sharey=True)
    plt.subplots_adjust(hspace=0.45)
    axes_list = [ax for ax_row in axes for ax in ax_row ]
    data = []
    for i, chroms in enumerate(chrom_homo_pairs):
        chr_data = []
        chrom_pairs = list(permutations(chroms, 2))
        homo_data = []
        
        for chrom in chroms:
            haps = (list(filter(lambda x: x[0]==chrom, chrom_pairs)))
            haps_data = []
            for chrom_pair in haps:
                chrom1, chrom2 = chrom_pair
                chrom_matrix = getChromPairsMatrix(cool, hm, chrom1, chrom2)
                chrom_matrix.sort(axis=1)
                haps_data.extend(method_func(chrom_matrix, axis=1))
            
            homo_data.append(haps_data)
        # homo_data.append(getChromPairsMatrix(cool, hm, chrom_pair[0], chrom_pair[1]))
        ax = axes_list[i]
        bplot = ax.boxplot(homo_data,
                showfliers=False, 
                patch_artist=True, 
                notch=True, 
                widths=0.35,
                medianprops=medianprops,
                whiskerprops=whiskerprops,
                boxprops=boxprops)
        for patch, color in zip(bplot['boxes'], ['#a83836', '#275e8c', '#df8384', '#8dc0ed']):
            patch.set_facecolor(color)
        ax.set_xticklabels(chroms, rotation=45)
    fig.text(0.02, 0.5, 
            "Interaction Frequency\n(contacts / {})".format(resolution_string), 
            rotation=90, va='center')   
    
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.savefig(output.rsplit('.', 1)[0] + '.png', dpi=300, bbox_inches='tight')


if __name__ == "__main__":
    main(sys.argv[1:])