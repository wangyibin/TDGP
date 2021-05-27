#!/usr/bin/env python
#! -*- coding:utf-8 -*-


from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys
import scipy
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
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

def plotMultiLR(data, species, output='out_multiSpecies', workdir='./', ncols=5):
    """
    plot multi linregression scatter picture
    
    Params:
    --------
    data: `dict`
        data of pc values after synteny analysis
    output: `str`
        output figname [default: out_multiSpecies]
    workdir: `str`
        path to workdir [default: ./]
    ncols: int
        columns of picture [default: 5]
    
    Returns:
    --------
    out: picutre
    
    Examples:
    --------
    >>> data = getPCValues(species)
    >>> plotMultiLR(data)
    """
    nrows = int(math.ceil(len(data)*1.0/ncols))
    ncols = min(ncols, len(data))
    if ncols < 2:
        fig, ax = plt.subplots(nrows, ncols, figsize=(ncols*5.2, nrows*5))
        axes = listify(ax)
    else:
        fig, axes = plt.subplots(nrows, ncols , figsize=(ncols*5.2, nrows*5))
        plt.subplots_adjust(wspace=0.28)
        axes = trim_axes(axes, len(data))
    for i, pairs in enumerate(data):
        bg1, bg2 = data[pairs]
        
        ax = plotLineRegress(axes[i], bg1.score, bg2.score, 
                                 xlabel=pairs[0], ylabel=pairs[1])
    savefig("{}/{}".format(workdir, output), bbox_inches='tight')

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
    pOpt.add_argument('-n', '--ncols', type=int, default=5, 
            help='columns number of picure [default: %(default)s]')
    pOpt.add_argument('-o', '--output', default="out_multiSpecies", 
            help='output picture name [default: out_multiSpecies]')
    args = p.parse_args(args)

    data = getPCValues(species=args.species)
    plotMultiLR(data, species=args.species, output=args.output,
            ncols=args.ncols)


if __name__ == "__main__":
    main(sys.argv[1:])