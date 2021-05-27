#!/usr/bin/env python
#! -*- coding:utf-8 -*-


from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys
import scipy
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import math

from collections import defaultdict, OrderedDict, Counter
from itertools import permutations, combinations

from TDGP.apps.base import listify
from TDGP.formats.bedGraph import import_bedgraph
from TDGP.graphics.ploty import plotLineRegress,savefig
from TDGP.graphics.ploty import trim_axes

def getPCValues(species, workdir='./'):
    """
    get pc values from a paired results
    
    Params:
    --------
    species: `list`
        list of species names
    workdir: `str` 
        directory of synteny compartments analysis
    
    Returns:
    --------
    data: `dict`
        results of pc dataframe
    
    Examples:
    --------
    >>> os.listdir(workdir)
    ["AP85-rice", "AP85-sorghum",...]
    >>> species = ['AP85', 'rice', 'sorghum' ...]
    >>> getPCValues(species)
    """
    data = OrderedDict()
    species_pairs = list(combinations(species, 2))
    for i, pairs in enumerate(species_pairs):
        data_dir = "{}-{}".format(*pairs)
        bg1 = "{}/{}.synteny.eigen1.nogene.bg".format(data_dir, pairs[0])
        bg2 = "{}/{}.synteny.eigen1.nogene.bg".format(data_dir, pairs[1])
        bg1 = import_bedgraph(bg1)
        bg2 = import_bedgraph(bg2)
        data[pairs] = (bg1, bg2)
    return data 


def getRvalueDf(data, species):
    """
    get a dataframe of rvalue for a correlation analysis
    
    Params:
    --------
    data: `dict`
        dict of bg data 
    species: `list`
        list of species name
    
    Returns:
    --------
    out: `DataFrame`
        `DataFrame of correlations matrix
    
    Examples:
    --------
    >>> species = ['AP85', 'rice', 'sorghum', 'foxtail', 'maize']
    >>> data = getPCValues(species)
    >>> getRvalueDf(data, species)
    """
    rvalue_list = []
    for pairs in data:
        bg1, bg2 = data[pairs]
        xdata, ydata = bg1.score, bg2.score
        _, _, rvalue, _, _ = scipy.stats.linregress(xdata, ydata)
        rvalue_list.append(rvalue)
    idx = list(combinations(range(len(species)), 2))
    matrix = np.zeros((len(species), len(species)))
    for i, r in zip(idx, rvalue_list):
        matrix[i[0], i[1]] = r
    matrix = matrix.T + matrix - np.diagonal(matrix)
    matrix[np.diag_indices_from(matrix)] = 1
    df = pd.DataFrame(matrix)
    df.columns = species
    df.index = species
    return df 


def plotClusterMap(df, output='out_heatmap', 
                   vmin=0, vmax=1, showCol=False,
                  showRow=True, cmap='Blues'):
    """
    plot clustermap for correlations
    
    Params:
    -------
    df: `pd.DataFrame`
        dataFrame of correlation matrix
    output: `str`
        output filename [default: out_heatmap.pdf]
    vmin: `int` or `float`
        minimum value of heatmap [deafult: 0]
    vmax: `int` or `float`
        maximum value of heatmap [default: 0]
    show{Col,Row}: `bool`
        show {col, row} tree strucure [default: showRow]
    
    Returns:
    -------
    out: clustermap figure
    
    Examples:
    --------
    >>> df
                AP85      rice   sorghum   foxtail     maize
        AP85     1.000000  0.012799  0.528472  0.303004  0.193328
        rice     0.012799  1.000000 -0.014729  0.018281  0.013439
        sorghum  0.528472 -0.014729  1.000000  0.475777  0.316860
        foxtail  0.303004  0.018281  0.475777  1.000000  0.122135
        maize    0.193328  0.013439  0.316860  0.122135  1.000000
    >>> plotClusterMap(df)
    """
#     plt.rcParams['xtick.labelsize'] = 18
#     plt.rcParams['ytick.labelsize'] = 18
#     plt.rcParams['xtick.alignment'] = 'center'
#     plt.rcParams['ytick.alignment'] = 'center'
    rc = {"xtick.labelsize": 18, "ytick.labelsize": 18,
                "xtick.alignment": "center", "ytick.alignment": "center"}
    sns.set(rc)
    cg = sns.clustermap(df, vmax=1, vmin=0,cmap=cmap, linewidths=1,
                      annot=True, annot_kws={"fontsize": 14}, fmt='.2f',
                   cbar_pos=(1, .25, .03, .35), tree_kws={'linewidth': 2})

    cg.ax_col_dendrogram.set_visible(showCol)
    cg.ax_row_dendrogram.set_visible(showRow)
    cax = plt.gca()
    cax.tick_params(labelsize=10)
    cax.set_ylabel('r', fontsize=16)
    savefig(output, bbox_inches='tight')


def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('species', nargs="+", 
            help='species name of plot')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    pOpt.add_argument('--cmap', default='Blues', 
            help='colormap of heatmap [default: %(default)]')
    pOpt.add_argument('-o', '--output', default="out_heatmap", 
            help='output picture name [default: out_multiSpecies]')
    args = p.parse_args(args)

    data = getPCValues(species=args.species)
    df = getRvalueDf(data, args.species)
    plotClusterMap(df, output=args.output, cmap=args.cmap)


if __name__ == "__main__":
    main(sys.argv[1:])