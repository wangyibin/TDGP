#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
%prog path_to_all_resolution_dir [options]

    calculate and plot the resolution of the HiC.
    it will produce two picture and an info_list.
    example:
        %prog hic_result/matrix/rice/raw

"""

from __future__ import print_function

import glob
import numpy as np
import pandas as pd
import os
import os.path as op
import sys

from multiprocessing import Pool, cpu_count
from scipy.sparse import coo_matrix, csr_matrix


def create_dense_matrix(matrix_file, to_dense=True):
    df = pd.read_csv(matrix_file, sep="\t",header=None)
    i, j, data = df.T.to_numpy()
    i = i.astype(int)
    j = j.astype(int)
    if min(i.min(), j.min()) == 0:
        N = max(i.max(), j.max()) + 1
        matrix = coo_matrix((data, (i, j)), shape=(N, N), dtype=float)
    else:
        N = max(i.max(), j.max())
        matrix = coo_matrix((data, (i-1, j-1)), shape=(N, N), dtype=float)

    if to_dense:
        matrix = np.array(matrix.todense())
        matrix = matrix + matrix.T
        matrix[np.diag_indices_from(matrix)] /= 2
        return matrix
    else:
        return csr_matrix(matrix)


def sum_array2df(array):
    array = np.sum(array, axis=0)
    return pd.DataFrame(array)


def count_num(df, values=[.25, .2, .1]):

    res = df[0].describe(values).astype(int)
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

def parallel(target, args, threads=cpu_count):
    p = Pool(min(len(args), threads))
    res = p.map(target, args)
    return res

def count_all_matrix(matrix_path):
    """
    count all matrix frequency
    """
    if not op.exists(matrix_path):
        print("[ERROR] no such file of  matrix_path ")
        sys.exit()

    matrix_list = glob.glob("{}/*/*.matrix".format(matrix_path))

    # array_list = parallel(create_dense_matrix, matrix_list)
    # df_list = parallel(sum_array2df, array_list)
    # res_list = parallel(count_num, df_list)
    # bin_size_list = parallel(get_bin_size, matrix_list)
    array_list = list(map(create_dense_matrix, matrix_list))
    df_list = list(map(sum_array2df, array_list))
    res_list = list(map(count_num, df_list))
    bin_size_list = list(map(get_bin_size, matrix_list))
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


def plot(df, bar_color=('#48D1CA','#F5CA14','#014D9B'), yticks=None,
         width=.2, out='resolution.pdf'):
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    color75, color80, color90 = bar_color
    
    p1 = plt.bar(df["bin_size"], df["75%"], width=width,color=color75)
    p2 = plt.bar(df['bin_size'], df["80%"], width=width,color=color80)
    p3 = plt.bar(df['bin_size'], df["90%"], width=width,color=color90)

    plt.legend((p1[0],p2[0],p3[0]), ('75%', '80%', '90%'))


    # sns.set_style('white')
    # sns.set_context({"figure.figsize": (8,5)})
    # sns.barplot(x="bin_size", y="25%", data=df, color=color75, linewidth=0,hue='floor')
    # sns.barplot(x="bin_size", y="20%" ,data=df, color=color80, linewidth=0, hue='floor')
    # sns.barplot(x="bin_size", y="10%", data=df, color=color90, linewidth=0, )
    # bar75 = plt.Rectangle((0,0),1,1,fc=color75, edgecolor='none')
    # bar80 = plt.Rectangle((0, 0), 1, 1, fc=color80, edgecolor='none')
    # bar90 = plt.Rectangle((0, 0), 1, 1, fc=color90, edgecolor='none')
    # l = plt.legend([bar75, bar80, bar90],['75%', '80%', '90%'], loc=2, ncol=1, prop={'size': 10})
    # l.draw_frame(False)
    if yticks is not None:
        yticks = eval(yticks)
        start, end, step = yticks
        plt.ylim(start, end)
        plt.yticks(range(start, end+1, step), range(start,end+1,step))

    plt.xlabel("Resolution")
    plt.ylabel("Depth")

    plt.savefig(out)
    plt.clf()

if __name__ == "__main__":
    from optparse import OptionParser

    p = OptionParser(__doc__)

    p.add_option('--yticks', dest='yticks',default='(0, 10001, 1000)',
                 help='yticks parameter (start, end, step), '
                      'e.g. (0,2000,100), [default: %default')
    p.add_option('--yticks2', dest='yticks2', default=None,
                 help='yticks parameter (start, end, step), '
                      'e.g. `(0,2000,100)`, [default: None]')
    p.add_option("--width", dest='width', default=.2, type=float,
                 help='the width of bar [default: %default]')
    p.add_option('--colors', dest='colors', default="('#3E466D','#E5E4DE','#ADCCD3')",
                 help="the colors of bar e.g. ('#48D1CA','#F5CA14','#014D9B')"
                      "[default: %default]")
    p.add_option('-t','--thread', dest='thread', default=4,type=int,
                help='the thread of script [default: %default]')
    p.add_option('--list', dest='list')
    p.add_option('-v','--values', dest='values', default=[.90, .80, .75],
                 help='values of describe [default: %default]')

    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(p.print_help())

    matrix_path, = args
    yticks = opts.yticks
    yticks2 = opts.yticks2
    width = opts.width
    colors = eval(opts.colors)

    df = count_all_matrix(matrix_path)
    df['bin_size'] = df['bin_size'].apply(bin_size_to_string)
    values = {"10%":"90%", "20%":"80%", "25%": "75%"}
    df = df.rename(columns=values)
    # out df to a file
    df.to_csv("out.info",sep="\t", index=False)

    # plot 2 figures
    plot(df, yticks=yticks, out='resolution1.pdf', width=width, bar_color=colors)
    plot(df, yticks=yticks2, out='resolution2.pdf', width=width,bar_color=colors)
