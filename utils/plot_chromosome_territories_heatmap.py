#!/usr/bin/env python3
# -*- coding:utf-8 -*-

"""
%(prog)s infile1 [infile2] -c chrom.sizes

    plot the heatmap of chromosome trans interactions to 
        show the chromosomal territories.
    also can plot two allpairs files in a heatmap. 
"""

import argparse
import numpy as np
import os.path as op
import sys
import json

from collections import defaultdict
from math import log2



def create_trans_db(allpairs):
    prefix = allpairs.split("_")[0]
    db2 = defaultdict(lambda : 0)
    db = defaultdict(lambda : db2.copy())
    with open(allpairs) as fp:
        for line in fp:
            chrom1 = line.split()[1]
            chrom2 = line.split()[4]
            if chrom1 != chrom2:
                db[chrom1][chrom2] += 1
                db[chrom2][chrom1] += 1
    

    json.dump(db, open('{}_trans.json'.format(prefix), 'w'))


def cacl_total_trans(chrom_trans_nums):
    total_trans = 0
    for chrom in chrom_trans_nums:
        total_trans += sum(list(map(int, chrom_trans_nums[chrom].values())))
    return  total_trans


def get_chrom_list(chrom_list):
    if op.exists(chrom_list):
        chrom_list = [i.strip().split()[0] 
                for i in open(chrom_list) if i.strip()]
    else:
        chrom_list = chrom_list.strip(",").split(",")
    
    return chrom_list



def get_whole_obs_exp_matrix(allpairs, chrom_list):
    prefix = allpairs.split("_")[0]

    if not op.exists('{}_trans.json'.format(prefix)):
        create_trans_db(allpairs)
    else:
        print('[Warning]: There has already existed'
                ' `{}_trans.json`, using old json file.'.format(prefix))
    chrom_trans_nums = json.load(open('{}_trans.json'.format(prefix)))
    
    whole_obs_exp_matrix = np.zeros((len(chrom_list), len(chrom_list)))
    total_trans = cacl_total_trans(chrom_trans_nums)

    for i, chrom1 in enumerate(chrom_list):
        chrom1_trans = sum(list(map(int, chrom_trans_nums[chrom1].values())))
        for j, chrom2 in enumerate(chrom_list):
            if i==j:
                continue
            chrom2_trans = sum(list(map(int, chrom_trans_nums[chrom2].values())))
            obs = float(chrom_trans_nums[chrom1][chrom2])
            exp = (((chrom1_trans / total_trans) * (chrom2_trans / (total_trans - chrom1_trans)) + (chrom2_trans / total_trans) * (chrom1_trans / (total_trans - chrom2_trans))) * (total_trans / 2))
            whole_obs_exp_matrix[i][j] = log2(obs / exp)


    return whole_obs_exp_matrix


def merge_matrix(matrix1, matrix2):
    """
    merge two matrix to a matrix spearate by up and down triangle
    """
    down = np.tril(matrix1, k=-1)
    up = np.triu(matrix2, k=1)
    merge = up + down
    return merge



def plot_heatmap(whole_obs_exp_matrix, chrom_list, prefix,
        color='coolwarm', valfmt='{x: .3f}', figsize=(10, 10), 
        vmin=0.2, vmax=-0.2):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    
    # set xtickslabel move to the top
    plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
    plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True

    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(whole_obs_exp_matrix, cmap=color, vmax=vmax, vmin=vmin)
    ax.figure.colorbar(im, ax=ax,fraction=0.045 )

    ax.set_xticks(np.arange(len(chrom_list)))
    ax.set_xticklabels(chrom_list, fontsize=14, rotation=45, ha='left')
    ax.set_yticks(np.arange(len(chrom_list)))
    ax.set_yticklabels(chrom_list, fontsize=14)

    valfmt = mpl.ticker.StrMethodFormatter(valfmt)

    for i in range(len(chrom_list)):
        for j in range(len(chrom_list)):
            if i == j:
                continue
            text = ax.text(j, i, valfmt(whole_obs_exp_matrix[i, j], None),
                      ha='center', va='center')

    plt.savefig('{}_whole_chromosome_positions.pdf'.format(prefix), dpi=300)




if __name__ == "__main__":
    p=argparse.ArgumentParser(prog=op.basename(__file__),
                        description=__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('infile', nargs="+", help='allValidpairs file')
    pReq.add_argument('-c', '--chrom', help='the chromsome size file')
    pOpt.add_argument('--vmin', type=float, default=-0.2, 
            help="min value of heatmap [default: %(default)s]")
    pOpt.add_argument('--vmax', type=float, default=0.2,
            help='max value of heatmap [default: %(default)s]')
    pOpt.add_argument('--figsize', type=str, default='(10, 10)', 
            help='figsize of picture [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args()
    chrom_list = get_chrom_list(args.chrom)
    if len(args.infile) == 1: 
        allpairs, = args.infile
        matrix = get_whole_obs_exp_matrix(allpairs, chrom_list)
        prefix = allpairs.split(".")[0]
    elif len(args.infile) == 2:
        allpairs1, allpairs2 = args.infile
        matrix1 = get_whole_obs_exp_matrix(allpairs1, chrom_list)
        matrix2 = get_whole_obs_exp_matrix(allpairs2, chrom_list)
        matrix = merge_matrix(matrix1, matrix2)
        prefix = allpairs1.split(".")[0] + "-" + allpairs2.split("_")[0]
    else:
        print('[Error]: infiles number should be 1 or 2')
        sys.exit()
    plot_heatmap(matrix, chrom_list, prefix, vmax=args.vmax, 
            vmin=args.vmin, figsize=eval(args.figsize))
