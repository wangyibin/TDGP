#!/usr/bin/env python
# -*- coding:utf-8 -*-


from __future__ import print_function

import argparse

import os.path as op
import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
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
        scale=1,
        yMin=0,
        yMax=0.6
        ):
    plt.rc('font', size=14)
    df['value'] = df['value'] / scale
    sns.scatterplot(x='pca1', y='value', data=df[df['pca1'] > 0], 
                    color='#a83836', axes=ax, label='A')
    sns.scatterplot(x='pca1', y='value', data=df[df['pca1'] < 0], 
                    color='#265e8a', axes=ax, label='B')
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(yMin, yMax)
    plt.legend()
    plt.savefig(output, bbox_inches='tight')

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
    pOpt.add_argument('--xLabel', default='PCA1',
            help='xLabel of picture [default: %(default)s]')
    pOpt.add_argument('--yLabel', default='Count',
            help='yLabel of picture [default: %(default)s]')
    pOpt.add_argument('--scale', default=1, type=int,
            help='scale of yticks [default: (default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    fig, ax = plt.subplots(figsize=(5, 5))
    df = import_ab_data(args.bg5file)
    output = args.output if args.output \
                else op.basename(args.bg5file).split('.')[0] + '_dotplot.pdf'
    plot(ax, df, output, 
            xlabel=args.xLabel, 
            ylabel=args.yLabel, 
            scale=args.scale)
    


if __name__ == "__main__":
    main(sys.argv[1:])