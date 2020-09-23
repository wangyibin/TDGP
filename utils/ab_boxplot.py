#!/usr/bin/env python
# -*- coding:utf-8 -*-


"""
%prog rice_eigen1_gene.bg [options]
    
    To plot boxplot of gene density of per A/B compartment.
    bg file: chrom start end eigen1 geneDensity


    ** if using xlabel or ylabel options, 
        you can also use latex to set special format (italic, up ...)
        italic: `"$\mathic{{trans}}$"`
        scientific count: `"$10 \times 6$"`
"""

from __future__ import print_function

import numpy as np
import math
import os.path as op
import pandas as pd
import sys


def wi_test(data1, data2):
    """
    Wilcoxon rank-sum tests
    return: pvalue
    """
    from scipy import stats
    wi = stats.ranksums(data1, data2)
    return wi.pvalue


def as_si(x, ndp):
    """
    Scientific count format function, use x replace e.
    """
    s = '{x:0.{ndp:d}e}'.format(x=x, ndp=ndp)
    m, e = s.split('e')
    return r'{m:s}\times 10^{{{e:d}}}'.format(m=m, e=int(e))



def ab_boxplot(a_data, b_data, chrom, ax, 
            scale=100, xlabel='', ylabel=''):
    
    # wi_test
    pvalue = wi_test(a_data[4], b_data[4])

    # property
    boxprops_r = dict(facecolor='#a83836',color='black', linewidth=1.5)
    boxprops_b = dict(facecolor='#265e8a',color='black', linewidth=1.5)
    medianprops=dict(color='black', linewidth=2.5)
    whiskerprops = dict(linestyle='--')


    #ax.set_title('Gene Density')

    a = ax.boxplot(a_data[4]/scale, showfliers=False, patch_artist=True, notch=True, widths=0.22,
               boxprops=boxprops_r, medianprops=medianprops, whiskerprops=whiskerprops)
    a_upper_extreme = [ item.get_ydata() for item in a['whiskers']][1][1]
    a_bottom_extreme = [ item.get_ydata() for item in a['whiskers']][0][1]

    b = ax.boxplot([[], b_data[4]/scale], showfliers=False, patch_artist=True, notch=True, widths=0.22,
               boxprops=boxprops_b, medianprops=medianprops, whiskerprops=whiskerprops)
    b_upper_extreme = [ item.get_ydata() for item in b['whiskers']][3][1]
    b_bottom_extreme = [ item.get_ydata() for item in b['whiskers']][2][1]

    max_upper = max(a_upper_extreme, b_upper_extreme)
    min_upper = min(a_upper_extreme, b_upper_extreme)
    min_bottom = min(a_bottom_extreme, b_bottom_extreme)
    h = max_upper/5

    if max_upper == a_upper_extreme:
        y1 = max_upper + max_upper/10
        y2 = min_upper + max_upper/10
    else:
        y1 = min_upper + max_upper/10
        y2 = max_upper + max_upper/10
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['A', 'B'], fontsize=12)
    ax.set_xlabel(r"{}".format(xlabel), fontsize=13)
    ax.set_ylabel(r"{}".format(ylabel), fontsize=13)
    ax.set_ylim((min_bottom - max_upper/5, max_upper + max_upper/2))
    ax.plot([1,1,2,2], [y1, max_upper + h, max_upper + h, y2], linewidth=1.0, color='black')
    ax.text((1 + 2)*.5, max_upper + max_upper/4.5, r'$P = {0:s}$'.format(as_si(pvalue, 2)), ha='center', va='bottom' )


def trim_axs(axs, N):
    """
    little helper to message the axs list to have correct length...
    """

    axs = axs.flat
    for ax in axs[N:]:
        ax.remove()
    
    return axs[:N]


def read_bg(bgfile):

    ab_data = pd.read_csv(bgfile, header=None, sep='\t')
    
    return ab_data


def plot(bgfile, scale=100, title="Gene Density", outfile=None, chrom_list='all', column_num=4, 
        draw_per_chrom=False, exclude=[], xlabel='', ylabel='',
        sort_func='lambda x: int(x[3:])'):
    
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    
    ab_data = read_bg(bgfile)
    
    outfile = outfile if outfile else bgfile.rsplit('.', 1)[0] + ".pdf"
    
    if not isinstance(exclude, list):
        exclude = exclude.split(',')
        

    if chrom_list == 'all':
        chrom_list = [i for i in set(ab_data[0]) if i not in exclude ]
        ab_data = ab_data[ ab_data[0].isin(chrom_list) ]
    else:
        if op.exists(chrom_list):
            chrom_list = [i.strip() for i in open(chrom_list)]
        else:
            chrom_list = chrom_list.split(',')
        chrom_list = [ i for i in chrom_list if i in set(ab_data[0]) ]
        ab_data = ab_data[ ab_data[0].isin(chrom_list) ]

    
    if not draw_per_chrom:
        fig, ax = plt.subplots(figsize=(4,4))
        chrom = xlabel
        suptitle_props = dict(fontsize=16, x=0.5, y=0.95)
        a_data = ab_data[ab_data[3] > 0]
        b_data = ab_data[ab_data[3] < 0]
        ab_boxplot(a_data, b_data, chrom, ax, scale, xlabel, ylabel)
    else:
        suptitle_props = dict(fontsize=24, x=0.5, y=1.04)
        chrom_list.sort(key=eval(sort_func))
        fig1, axs = plt.subplots(int(math.ceil(len(chrom_list)*1.0/column_num)),
                    column_num, figsize=(column_num*4 + 2,
                    int(math.ceil(len(chrom_list)*1.0/column_num))*4 + 1), 
                    constrained_layout=True)
        axs = trim_axs(axs, len(chrom_list))

        for ax, chrom in zip(axs, chrom_list):
            xlabel = r"{} ({})".format(xlabel, chrom)
            a_data = ab_data[(ab_data[0] == chrom) & (ab_data[3] > 0)]
            b_data = ab_data[(ab_data[0] == chrom) & (ab_data[3] < 0)]
            ab_boxplot(a_data, b_data, chrom, ax, scale, xlabel, ylabel)

    plt.suptitle(title, **suptitle_props)
    plt.savefig(outfile, dpi=300,bbox_inches='tight' )
    plt.savefig(outfile.rsplit('.', 1)[0] + '.png', 
                dpi=300,
                bbox_inches='tight')  



if __name__ == "__main__":
    from optparse import OptionParser

    p = OptionParser(__doc__)
    
    p.add_option('-t', '--title', default="",
            help='the title of picture')
    p.add_option('-o', '--outfile', default=None, 
                    help='outfile of plot [default: %default]')
    p.add_option('--chrom_list', default='all',
            help='list of chromosome which plot, split by comma'
            ' or a list file.[default: %default]')
    p.add_option('-n', '--column_num', default=4, type=int,
            help='the number of subplot column [default: %default]')
    p.add_option('-c', '--perchrom', default=False,
            action='store_true', help='draw boxplot per chromosome'
            '[default: %default]')
    p.add_option('--exclude', default=[],
            help='list of chromosome which exclude, split by comma'
            '[default: %default]')
    p.add_option('--scale', default=1, type=int, 
            help='the scale of yticks [default: %default]')
    p.add_option('--xlabel', default='', 
            help='xlabel of genome wide picture [default: none]')
    p.add_option('--ylabel', default='',
            help='ylabel of genome wide picture [default:none]')
    p.add_option('--sort_func', default='lambda x: int(x[3:])',
            help='chromosome name sort function [default: %default]')


    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(p.print_help())

    bgfile, = args
    title = opts.title
    outfile = opts.outfile
    chrom_list = opts.chrom_list
    column_num = opts.column_num
    draw_per_chrom = opts.perchrom
    exclude = opts.exclude
    scale = opts.scale
    xlabel = opts.xlabel
    ylabel = opts.ylabel
    sort_func = opts.sort_func
    
    plot(bgfile, scale, title, 
        outfile, chrom_list,column_num, 
        draw_per_chrom, exclude, 
        xlabel, ylabel, sort_func)
