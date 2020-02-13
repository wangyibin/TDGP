#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
%prog sample_allValidPairs chrom_reference.sizes [Options]
"""

from __future__ import print_function

import os
import sys
import pandas as pd
import numpy as np
import multiprocessing as mp


def cacl_100bp_counts(infile, outfile=None):
    if not outfile:
        outfile = infile.rsplit('.')[0] + ".100.counts"
    if not os.path.exists(outfile):
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


def multi_cacl(df,chrom_dict, threads=10):    
    values = []
    task_list = [(df, chrom_dict, bin) for bin in range(10, 1001, 10)]
    pool=mp.Pool(threads)
    res = pool.map(merge_windows, task_list)
    res = dict(res)

    for i in sorted(res):
        values.append(res[i])
    print(values)
    return np.array(values)


def dotplot(values, outprefix, percent=0.8):
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
    plt.plot([0,100],[percent, percent], '--', color='#BB4853') 
    plt.xticks(range(0,101,10), map(lambda x: str(x), range(0, 101, 10)))
    plt.xlim(0,100)
    plt.ylim(0.0, 1.2)
    plt.xlabel("Resolution (kb)", fontsize=12)
    plt.ylabel("Percent (X100%)", fontsize=12)
    plt.title("Hi-C")
    plt.savefig('{}.pdf'.format(outprefix), dpi=300)

def main(infile, chrom_list, threads=20):
    outprefix = infile.split("_")[0]
    print(outprefix)
    chrom_dict = dict(i.strip().split() for i in open(chrom_list) if i.strip())
    print(chrom_dict)
    infile = cacl_100bp_counts(infile)
    df = import_data(infile)
    values = multi_cacl(df, chrom_dict, threads)
    dotplot(values, outprefix)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("*"*79)
        print("Usage: {} sample.allValidPairs chrom.sizes".format(os.path.basename(sys.argv[0])))
        print("*"*79)

        sys.exit()

    main(sys.argv[1], sys.argv[2])
