#!/usr/bin/env python
# -*- coding:utf-8 -*-


from __future__ import print_function

import argparse

import os.path as op
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns



def import_ab_data(bg5file):
    df = pd.read_csv(bg5file, sep='\t', 
                    header=None, 
                    names=['chrom', 'start', 
                            'end', 'pca1', 
                            'value'])
    
    return df


def plot(ax, df, 
        output,
        xlabel='PCA1', 
        ylabel='Count',
        title='gene',
        scale=1,
        yMin=None,
        yMax=None,
        log1p=True,
        topA=0,
        topB=0
        ):
    plt.rc('font', size=14)
    df['value'] = df['value'] / scale
    if log1p:
        df['value'] = np.log1p(df['value'])
        ylabel = "log({} + 1)".format(ylabel)
    
    a_data = df[df['pca1'] > 0]
    b_data = df[df['pca1'] < 0]  

    sns.scatterplot(x='pca1', y='value', data=a_data,
                    color='#a83836', axes=ax, label='A')
    sns.scatterplot(x='pca1', y='value', data=b_data, 
                    color='#265e8a', axes=ax, label='B')
    
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.set_title(title, fontsize=15)

    if yMax or yMin:
        ax.set_ylim(yMin, yMax)

    if topA != 0:
        a_data_sorted = a_data.sort_values(by=['value'], ascending=False)
        top_data = a_data_sorted.iloc[:topA, :]
        topA_output = output.rsplit('.', 1)[0] + "_topA{}.tsv".format(topA)
        top_data.to_csv(topA_output, sep='\t', header=None, index=None)
        ax.hlines(top_data['value'].iloc[-1], 0, a_data_sorted['pca1'].max(),
                linestyles='--', color='#bcbcbc')
    if topB != 0:
        b_data_sorted = b_data.sort_values(by=['value'], ascending=False)
        top_data = b_data_sorted.iloc[:topB, :]
        topB_output = output.rsplit('.', 1)[0] + "_topB{}.tsv".format(topB)
        top_data.to_csv(topB_output, sep='\t', header=None, index=None)
        ax.hlines(top_data['value'].iloc[-1], b_data_sorted['pca1'].min(), 0,
                linestyles='--', color='#bcbcbc')
    
    
    plt.legend()
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.savefig(output.rsplit('.', 1)[0] + '.png', dpi=300, bbox_inches='tight')

def main(args):
    """
    Plotting the scatter of A/B compartment by pca1 in x-axis,
        values in y-axis
    """
    p = argparse.ArgumentParser(prog=main.__name__,
                            description=main.__doc__,
                            conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bg5file', help='bg file of eigen1 and value')
    
    pOpt.add_argument('-o', '--output', default=None,
            help='output file [default: inputfile.pdf]')
    pOpt.add_argument('--xlabel', default='PCA1',
            help='xLabel of picture [default: %(default)s]')
    pOpt.add_argument('--ylabel', default='Count',
            help='yLabel of picture [default: %(default)s]')
    pOpt.add_argument('--title', default='gene', 
            help='title of picture [default: %(default)s]')
    pOpt.add_argument('--scale', default=1, type=int,
            help='scale of yticks [default: (default)s]')
    pOpt.add_argument('--log1p', default=False, action='store_true',
            help='Plot the log1p of the values.')
    pOpt.add_argument('--topA', default=0, type=int,
            help='output topA number of high expression compartments'
                'and plotting hlines in pictures [default: %(default)s]')
    pOpt.add_argument('--topB', default=0, type=int,
            help='output topB number of high expression compartments'
                'and plotting hlines in pictures [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    fig, ax = plt.subplots(figsize=(5, 5))
    df = import_ab_data(args.bg5file)
    output = args.output if args.output \
                else op.basename(args.bg5file).split('.')[0] + '_dotplot.pdf'
    plot(ax, df, output, 
            xlabel=args.xlabel, 
            ylabel=args.ylabel,
            title=args.title, 
            scale=args.scale,
            log1p=args.log1p,
            topA=args.topA,
            topB=args.topB)
    


if __name__ == "__main__":
    main(sys.argv[1:])