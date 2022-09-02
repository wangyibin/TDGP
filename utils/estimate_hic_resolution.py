#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
%(prog)s <sample_allValidPairs> <chrom_reference.sizes> [Options]
    To estimate the resolution of hic maps
"""

from __future__ import print_function

import argparse
import os
import os.path as op
import sys
import pandas as pd
import numpy as np
import multiprocessing as mp


def cacl_100bp_counts(infile, outfile=None, outdir="./"):
    outdir = op.abspath(outdir)
    if not outfile:
        outfile = outdir + "/" + op.basename(infile).rsplit('.')[0] + ".100.counts"
    
    if not os.path.exists(outfile) or not op.getsize(outfile):
        cmd = """cat {} | awk '{{a[$2][int($3/100)]++;a[$5][int($6/100)]++}}END{{for (j in a) for (i in a[j]) printf "%s\\t%s\\t%d\\n",j,i,a[j][i]}}' > {}""".format(infile, outfile)
        print(cmd)
        os.system(cmd)
    return outfile


def import_data(infile):
    df = pd.read_csv(infile, sep='\t', header=None)
    df = df.sort_values(by=[0, 1])
    df.set_index(keys=[0,1], inplace=True)
    return df


def merge_windows(args):
    df, chrom_dict, bin = args
    total = 0
    upper = 0
    for chrom in chrom_dict:
        df_chrom = df.loc[chrom]
        chrom_len = int(chrom_dict[chrom])
        df_chrom = df_chrom.reindex(range(chrom_len), fill_value=0)
        values = np.array(df_chrom[2])

        for i in range(0, chrom_len//100 + bin, bin):
            if values[i: i+bin-1].sum() > 1000:
                upper += 1
            total += 1
    return bin, upper*1.0 / total


def multi_cacl(df, chrom_dict, threads=10):    
    values = []
    task_list = [(df, chrom_dict, bin) for bin in range(10, 1001, 10)]
    pool=mp.Pool(threads)
    res = pool.map(merge_windows, task_list)
    res = dict(res)

    for i in sorted(res):
        values.append(res[i])

    return np.array(values)


def dotplot(values, outprefix, percent=0.8, outdir='./'):
    outdir = op.abspath(outdir)
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt

    values_below = values[values < percent]
    values_below_x = np.where(values < percent)[0]
    values_upper = values[values >= percent]
    values_upper_x = np.where(values >= percent)[0]
    plt.plot(values_upper_x, values_upper, '.', color='#89023E', label='Achieved resolution')
    plt.plot(values_below_x, values_below, '.', color='#b4b4b4', label='Unattained resolution')
    plt.legend(loc=0)
    plt.plot([0, 100],[percent, percent], '--', color='#BB4853') 
    plt.xticks(range(0, 101, 10), map(lambda x: str(x), range(0, 101, 10)))
    plt.xlim(0, 100)
    plt.ylim(0.0, 1.2)
    plt.xlabel("Resolution (kb)", fontsize=12)
    plt.ylabel("Percent (X100%)", fontsize=12)
    plt.title("Hi-C")
    plt.savefig('{}/{}.pdf'.format(outdir, outprefix), dpi=300)
    plt.savefig('{}/{}.png'.format(outdir, outprefix), dpi=300)

def main(infile, chrom_list, 
        outdir='./', threads=10):

    outprefix = op.basename(infile).split("_")[0]
    
    chrom_dict = dict(i.strip().split() 
                        for i in open(chrom_list) if i.strip())
    print(outprefix)
    infile = cacl_100bp_counts(infile, outdir=outdir)
    df = import_data(infile)
    print(chrom_dict)
    values = multi_cacl(df, chrom_dict, threads)
    dotplot(values, outprefix, outdir=outdir)


if __name__ == "__main__":
    p = argparse.ArgumentParser(prog=op.basename(__file__),
                        description=__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('allValidPairs', help='allValidParis file')
    pReq.add_argument('chromsize', help='chromsize file')
    pOpt.add_argument('-t', '--thread', type=int, default=10, 
            help='thread numbers of program [default: %(default)s]')
    pOpt.add_argument('-o', '--outdir', default='./',
            help='outdir of results [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args()
    

    main(args.allValidPairs, 
            args.chromsize, 
            args.outdir, 
            args.thread)
