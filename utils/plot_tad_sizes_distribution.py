#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
%prog infile.list [Options]
    
    To plot the sizes distribution of tad.
    Example:
        %prog infile.list
        %prog infile1,infile2,infile3
        %prog infile.list --xlim="(100,2000)"

"""


import numpy as np
import os
import os.path as op
import pandas as pd
import seaborn as sns
import sys



def import_data(infile):
    df = pd.read_csv(infile, header=None, sep='\t')
    return df

def chrom_size_convert(size):
    exts = ["", 'kb', 'Mb', 'Gb']
    i = 0
    while size >= 1000:
        size = size // 1000
        i += 1
        return "{}{}".format(int(size), exts[i])


def plot_tad_sizes(sizes_dict, xlim=(50, 800)):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    
    xlim_start, xlim_end = xlim
    plt.figure(figsize=(6,5))
    for species in sizes_dict:
        size = sizes_dict[species]
        sns.distplot(size.values//1000, hist=False, kde=True, 
                label="{} ({})".format(species, len(size)), kde_kws={'linewidth':3})
    plt.xticks(np.linspace(xlim_start, xlim_end, 11),
            list(map(int, np.linspace(xlim_start, xlim_end, 11))), rotation=45, ha='center')
    plt.xlim(xlim_start, xlim_end)
    plt.xlabel('TAD Sizes (kb)')
    plt.ylabel('Density')

    prefix = "_".join(sizes_dict.keys())
    plt.savefig('{}_tad_sizes_distribution.pdf'.format(prefix), dpi=300)


def main(infile_list, xlim):
    if op.exists(infile_list):
        infile_list = [i.strip() for i in open(infile_list) if i.strip()]
    else:
        infile_list = infile_list.split(',')

    outprefix_list = list(map(lambda x: x.split("_")[0], infile_list))
    df_list = list(map(import_data, infile_list))
    sizes_list = list(map(lambda x: x[2] - x[1], df_list))
    sizes_dict = dict(zip(outprefix_list, sizes_list))
    
    xlim = eval(xlim)
    plot_tad_sizes(sizes_dict, xlim)

if __name__ == "__main__":
    from optparse import OptionParser

    p = OptionParser(__doc__)
    p.add_option("--xlim", default="(50, 800)",
            help="the xlim of xticks [default: %default]")

    opts, args = p.parse_args()
    
    if len(args) != 1:
        sys.exit(p.print_help())
    infile_list, = args
    

    main(infile_list, opts.xlim)
