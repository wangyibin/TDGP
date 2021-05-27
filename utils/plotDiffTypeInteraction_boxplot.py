#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
%(prog)s -m <samplt_1000000.cool> [Options]

    Plotting the different interaction type boxplot.
    
"""
from __future__ import print_function

import argparse
import cooler
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import re
import os.path as op
import os
import sys

from itertools import combinations
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
def getBetweenHomo(cool, 
                chrom_homo_pairs, 
                method='median'):
    method_func = np.median if method == 'median' else np.mean
    contact_data = []
    hm = cool2matrix(cool)
    for homo_chroms in chrom_homo_pairs:
        chrom_pairs = list(combinations(homo_chroms, 2))
        for chrom_pair in chrom_pairs:
            chrom1, chrom2 = chrom_pair
            chrom_hm = getChromPairsMatrix(cool, hm, chrom1, chrom2)
            contact_data.extend((method_func(chrom_hm, axis=1)))
    
    return contact_data

def getWithinChromosome(cool, method='median'):
    method_func = np.median if method == 'median' else np.mean
    chroms = cool.chromnames
    hm = cool2matrix(cool)
    contact_data = []
    for chrom in chroms:
        idx = cool.bins().fetch(chrom).index
        start = idx[0]
        end = idx[-1]
        contact_mean_data = method_func(hm[start: end + 1, start: end + 1], axis=1)
        contact_data.extend(contact_mean_data)
    return contact_data


def getBetweenHeter(cool, chrom_heter_pairs, method='median'):
    method_func = np.median if method == 'median' else np.mean
    contact_data = []
    hm = cool2matrix(cool)
    for chrom_pair in chrom_heter_pairs:
        chrom1, chrom2 = chrom_pair
        chrom_hm = getChromPairsMatrix(cool, hm, chrom1, chrom2)
        contact_data.extend(method_func(chrom_hm, axis=1))
    return contact_data

def plotBoxPlot(data, resolution_string, output, dpi=300):
    fig, ax = plt.subplots(figsize=(5.5, 5))
    boxprops = dict(color='black', linewidth=1.5)
    medianprops=dict(color='black', linewidth=2.5)
    whiskerprops = dict(linestyle='--')

    bplot = ax.boxplot(data, 
                       showfliers=False, 
                       patch_artist=True, 
                       notch=True, 
                       widths=0.35,
                       medianprops=medianprops,
                       whiskerprops=whiskerprops,
                       boxprops=boxprops)

    for patch, color in zip(bplot['boxes'], ['#a83836', '#df8384', '#8dc0ed']):
        patch.set_facecolor(color)
    ax.set_xticklabels(['Within\n chromosomes', 'Within\n homologs', 'Between\n homologs'], fontsize=12)
    ax.get_yaxis().set_tick_params(which='both', labelsize=12)
    ax.set_xlabel("Interaction type", fontsize=14, labelpad=12)
    ax.set_ylabel("Interaction Frequency\n(contacts / {})".format(resolution_string), fontsize=14)
    plt.savefig(output, bbox_inches='tight', dpi=dpi)
    plt.savefig(output.rsplit('.', 1)[0] + ".png", 
                bbox_inches='tight', dpi=dpi)

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-m', '--matrix', default=None, 
            help='matrix of contacts', required=True)
    pReq.add_argument('-o', '--output', required=True,
            help='output file of pictures')
    pOpt.add_argument('--method', choices=['average', 'median'], 
            default='median',
            help='method of interaction frequency per bin [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    coolfile = args.matrix
    cool = cooler.Cooler(coolfile)
    all_chroms = cool.chromnames
    homo_chroms = sorted(set(map(lambda x: x[:-1], all_chroms)))
    chrom_heter_pairs = list(filter(is_heter, list(combinations(all_chroms, 2))))
    chrom_homo_pairs = find_homo_chrom(homo_chroms, all_chroms)
    WithinChrom_data = getWithinChromosome(cool, args.method)
    BetweenHomo_data = getBetweenHomo(cool, chrom_homo_pairs, args.method)
    BetweenHeter_data = getBetweenHeter(cool, chrom_heter_pairs, args.method)
    resolution_string = chrom_size_convert(cool.binsize)
    data = [WithinChrom_data, BetweenHomo_data, BetweenHeter_data]

    plotBoxPlot(data, resolution_string, args.output)


if __name__ == "__main__":
    main(sys.argv[1:])