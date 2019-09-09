#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import os
import os.path as op
import sys
import pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

def read_data(infile):
    df = pd.read_csv(infile, header=None, sep='\t')
    return df


def get_median_link(syn_data):
    df = pd.DataFrame([syn_data[0], (syn_data[2] - syn_data[1] )//2 + \
                      syn_data[1], syn_data[3], (syn_data[5] - syn_data[4] )//2 + \
                     syn_data[4], syn_data[6], syn_data[7]], index=[0,1,2,3, 4, 5]).T
    return df


def get_non_conserved_syn_gene_data(syn_data, rice_data, sorghum_data, threshold=0.005):
    i_a2b, i_b2a = [], []

    for rchrom, rstart, rend, schrom, sstart, send, gene1, gene2 \
            in zip(syn_data[0], syn_data[1], syn_data[2], syn_data[3],
                   syn_data[4], syn_data[5], syn_data[6], syn_data[7]):
        rpos = (rend - rstart) // 2 + rstart
        spos = (send - sstart) // 2 + sstart
        ivalue = rice_data[(rice_data[0] == rchrom) & (rice_data[1] <= rpos) & (rice_data[2] > rpos)]
        jvalue = sorghum_data[(sorghum_data[0] == schrom) & (sorghum_data[1] <= spos) & (sorghum_data[2] > spos)]
        if ivalue.empty or jvalue.empty:
            continue

        if float(ivalue[3]) * float(jvalue[3]) < 0:
            if float(ivalue[3]) > float(jvalue[3]) and float(ivalue[3]) - float(jvalue[3]) >= threshold:
                i_a2b.append([rchrom, rstart, rend, gene1, schrom, sstart, send, gene2])
               # j_a2b.append([schrom, sstart, send, gene2])
            elif float(ivalue[3]) < float(jvalue[3]) and float(jvalue[3]) - float(ivalue[3]) >= threshold:
                i_b2a.append([rchrom, rstart, rend, gene1, schrom, sstart, send, gene2])
                #j_b2a.append([schrom, sstart, send, gene2])
    print(i_a2b[:10], i_b2a[:10])
    return i_a2b, i_b2a


def pie_plot(i_a2b, i_b2a, syn_data, outprefix='out'):
    colors = ('#5190bb', '#dba9a9', '#ccd9d9')
    # colors = ('#a83836','#265e8a', '#ccd9d9')
    # explode = (0.1, 0.1, 0)
    sizes = [len(i_a2b), len(i_b2a), len(syn_data) - len(i_a2b) - len(i_b2a), ]
    labels = ('A2B ({})'.format(len(i_a2b)), 'B2A ({})'.format(len(i_b2a)), 'Conserved ({})'.format(sizes[2]))
    plt.figure(figsize=(5, 5))
    plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%')
    plt.savefig('{}_AB_transform_pie.pdf'.format(outprefix), dpi=300)


def bar_plot(i_a2b, i_b2a, outprefix='out'):
    plt.figure(figsize=(5.5, 5))
    plt.bar([0, 0.25, 0.75, 1], [0, len(i_a2b), len(i_b2a), 0], width=.2, color=('#a83836', '#265e8a'))
    plt.xticks([0.25, 0.75], ["A2B", "B2A"])
    plt.xlabel('Type of transform')
    plt.ylabel('Number of genes')
    plt.savefig('{}_AB_transform_bar.pdf'.format(outprefix), dpi=300)


def result_visualization(a2b_data, b2a_data, syn_data, outprefix='out'):

    bar_plot(a2b_data, b2a_data, outprefix)
    pie_plot(a2b_data, b2a_data, syn_data, outprefix)


def out_result(rice_infile, sorghum_infile, syn_file,
               outprefix='out', visual=True, threshold=0.005):
    rice_data = read_data(rice_infile)
    sorghum_data = read_data(sorghum_infile)
    syn_data = read_data(syn_file)

    a2b_data, b2a_data = get_non_conserved_syn_gene_data(syn_data,
                                        rice_data, sorghum_data, threshold=threshold)
    with open("{}_a2b.txt".format(outprefix), 'w') as out:
        for line in a2b_data:
            out.write('\t'.join(map(str, line)) + '\n')

    with open("{}_b2a.txt".format(outprefix), 'w') as out:
        for line in b2a_data:
            out.write('\t'.join(map(str, line)) + '\n')

    if visual:
        result_visualization(a2b_data, b2a_data, syn_data, outprefix)


if __name__ == "__main__":
    from optparse import OptionParser
    p = OptionParser(__doc__)

    p.add_option('--visual', default=True, action='store_false',
                 help='whether to visual result [default: %default]')
    p.add_option('--threshold', default=0.005, type=float,
                 help='the threshold of AB swith [default: %default]')

    opts, args = p.parse_args()
    if len(args) != 4:
        sys.exit(p.print_help())

    rice_file, sorghum_file, syn_file, outprefix = args
    out_result(rice_file, sorghum_file, syn_file, outprefix, opts.visual, opts.threshold)
