#!/usr/bin/env python
# -*- coding:utf-8 -*-

from __future__ import print_function

import glob
import numpy as np
import multiprocessing
import pandas as pd
import os
import os.path as op
import sys



def read_matrix(matrix_file):
    df = pd.read_csv(matrix_file, sep="\t",header=None)
    return df.to_dense()

def count_num(df, values=[.25, .2, .1]):
    res = df[2].describe(values).astype(int)
    return res


def build_matrix(allpairs_file, chr_size_file, resolution, outpath,
                 soft_path='build_matrix', matrix_format='upper'):
    """
    use the script from hicpro to build matrix at the specified resolution.
    """
    if os.system(soft_path) != 0:
        print("[ERROR] command of {} not found".format(soft_path))
        sys.exit()

    outprefix = "{}_{}".format(outpath, resolution)
    cmd = '{0} --matrix-format {1} --binsize {2} --chrsizes {3}' \
          '--ifile {4} --oprefix {5}'.format(soft_path, matrix_format,
                                             resolution, chr_size_file,
                                             allpairs_file, outprefix)
    if os.system(cmd) != 0:
        print("[INFO] succeed to  build raw map of {}.matrix".format(outprefix))
    else:
        print("[ERROR] failed to build raw map of {}.matrix".format(outprefix))

    return '{}.matrix'.format(outprefix)

def get_bin_size(name):
    import re
    regex = re.compile('_(\d+).matrix')
    bin_size = regex.findall(op.basename(name))[0]

    return bin_size


def count_all_matrix(matrix_path):
    """
    count all matrix frequency

    """
    if not op.exists(matrix_path):
        print("[ERROR] no such file of  matrix_path ")
        sys.exit()

    matrix_list = glob.glob("{}/*/*matrix".format(matrix_path))
    df_list = list(map(read_matrix, matrix_list))
    bin_size_list = list(map(get_bin_size,matrix_list))
    res_list = list(map(count_num, df_list))
    res_dict = dict(zip(bin_size_list, res_list))
    df = pd.DataFrame(res_dict)
    df = df.T
    df["bin_size"] = df.index.astype('int')
    df.sort_values(['bin_size'], inplace=True)
    return df


def bin_size_to_string(bin_size):
    if bin_size < 1000000:
        return str(bin_size // 1000) + "K"
    elif bin_size >= 1000000:
        return str(bin_size // 1000000) + "M"


def plot(df):
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    bar_color = ('#48D1CA','#F5CA14','#014D9B')

    sns.set_style('white')
    sns.set_context({"figure.figsize": (6,4)})

    color75, color80, color90 = bar_color
    sns.barplot(x="bin_size", y="25%", data=df, color=color75, linewidth=0)
    sns.barplot(x="bin_size", y="20%" ,data=df, color=color80, linewidth=0)
    sns.barplot(x="bin_size", y="10%", data=df, color=color90, linewidth=0)
    sns.barplot()
    bar75 = plt.Rectangle((0,0),1,1,fc=color75, edgecolor='none')
    bar80 = plt.Rectangle((0, 0), 1, 1, fc=color80, edgecolor='none')
    bar90 = plt.Rectangle((0, 0), 1, 1, fc=color90, edgecolor='none')
    l = plt.legend([bar75, bar80, bar90],['75%', '80%', '90%'], loc=2, ncol=1, prop={'size': 10})
    l.draw_frame(False)
    plt.xlabel("Resolution")
    plt.ylabel("Depth")

    plt.savefig('t.pdf')


if __name__ == "__main__":
    from optparse import OptionParser

    p = OptionParser(__doc__)

    p.add_option('-v','--values', dest='values', default=[.90, .80, .75],
                 help='values of descrive [default: %default]')

    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(p.print_help())

    matrix_path, = args

    df = count_all_matrix(matrix_path)
    df['bin_size'] = df['bin_size'].apply(bin_size_to_string)
    plot(df)
