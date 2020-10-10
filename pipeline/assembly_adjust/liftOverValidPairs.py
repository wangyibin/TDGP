#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
liftover allValidPairs from hicpro
"""

from __future__ import print_function

import argparse
import logging
import gc
import os
import os.path as op
import sys
import multiprocessing as mp
import pandas as pd
import numpy as np
from collections import OrderedDict


ValidPairsHeader = ["read", "chr1", "pos1", "strand1",
                    "chr2", "pos2", "strand2", "size",
                    "site1", "site2"]

AGP_NAMES_tig = ['chrom', 'start', 'end', 'number',
                 'type', 'id', 'tig_start', 'tig_end', 'orientation']
AGP_NAMES_tig_dtype = {'chrom': '|S', 'start': np.int64, 
                    'end': np.int64, 'number': np.int32,
                 'type': '|S', 'id': '|S', 'tig_start': np.int32, 
                 'tig_end': np.int32, 'orientation': '|S'}
AGP_NAMES_gap = ['chrom', 'start', 'end', 'number',
                 'type', 'length', 'type', 'linkage', 'evidence']


def import_agp(agpfile, split=True):
    """
    import agp file and return a dataframe
    """
    df = pd.read_csv(agpfile, sep='\t',
                     header=None, index_col=None)
    if split:
        tig_df = df[df[4] == 'W']
        gap_df = df[df[4] == 'U']
        tig_df.columns = AGP_NAMES_tig
        gap_df.columns = AGP_NAMES_gap
        tig_df.set_index('chrom', inplace=True)
        gap_df.set_index('chrom', inplace=True)
        tig_df.loc[:, 'tig_end'] = tig_df['tig_end'].astype('int32')
        tig_df.loc[:, 'tig_start'] = tig_df['tig_start'].astype('int32')
        
        logging.debug('Load `{}`'.format(agpfile))
        return tig_df, gap_df
    else:
        return df

def import_validpairs(infile, chunksize=1000):
    df = pd.read_csv(infile, usecols=range(10), header=None,
                    sep='\t', index_col=None, 
                    names=ValidPairsHeader, 
                    chunksize=chunksize)
    return df


def get_new_coord(chrom, pos, old_agp_df, new_agp_df):
    """
    get a new coordinate according agp files.

    Params:
    --------
    chrom: `str` chromosome from validpairs
    pos: `int` position from validpairs
    old_agp_df: `dataframe` old agp dataframe 
    new_agp_df: `dataframe` new agp dataframe

    Returns:
    --------
    chrom: `str` chromosome
    new_pos: `int` new position

    Examples:
    --------
    >>> get_new_coord("Chr01g1", 1334, old_agp_df, new_agp_df)
    ('Chr01g1', 1334)
    """
    old_chrom_df = old_agp_df.loc[chrom]
    old_tig_res = old_chrom_df[(old_chrom_df['start'] <= pos) & 
                            (old_chrom_df['end'] >= pos) ]
    
    try:
        if len(old_tig_res) == 0:
            return chrom, np.nan
    except:
        return chrom, np.nan
    
    old_tig, = old_tig_res['id']
    old_tig_end, = old_tig_res['tig_end']
    old_tig_orientation, = old_tig_res['orientation']
    old_chrom_start, = old_tig_res['start']
    if old_tig_orientation == "+":
        old_tig_pos = pos - old_chrom_start + 1
    else:
        old_tig_pos = old_tig_end - (pos - old_chrom_start + 1)
    new_chrom_df = new_agp_df.loc[chrom]
    try:
        new_tig_res = new_chrom_df.loc[old_tig]
    except KeyError:
        return chrom, np.nan
    new_chrom_start = new_tig_res['start']
    new_chrom_end = new_tig_res['end']
    new_tig_orientation, = new_tig_res['orientation']

    if new_tig_orientation == '+':
        new_chrom_pos = new_chrom_start + old_tig_pos - 1
    else:
        new_chrom_pos = new_chrom_end - old_tig_pos
    del old_chrom_df, new_chrom_df
    gc.collect()
    return chrom, new_chrom_pos


def liftOverChunk(old_agp_df, new_agp_df, chunk_df ):
    
    chunk_df.chr1, chunk_df.pos1 = zip(*chunk_df.apply(lambda x: get_new_coord(x.chr1, 
                                            x.pos1, old_agp_df, new_agp_df), axis=1))
    chunk_df.chr2, chunk_df.pos2 = zip(*chunk_df.apply(lambda x: get_new_coord(x.chr2, 
                                            x.pos2, old_agp_df, new_agp_df), axis=1))
    chunk_df.dropna(inplace=True)
    chunk_df.loc[:, 'pos2'] = chunk_df['pos2'].astype('int32')
    return chunk_df


def liftOverValidPairs(args):
    """
    %(prog)s <validpairs> <old.agp> <new.agp> [Options]
        
        liftover validparis file to a new coordinate file. 

    """
    p = p=argparse.ArgumentParser(prog=liftOverValidPairs.__name__,
                        description=liftOverValidPairs.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('validpairs', 
            help='validpairs from hicpro')
    pReq.add_argument('old_agp', help='old agp file')
    pReq.add_argument('new_agp', help='new agp file')
    
    pOpt.add_argument('-t', '--threads', type=int, default=4,
            help='number of threads [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    old_agp_df, _ = import_agp(args.old_agp)
    new_agp_df, _ = import_agp(args.new_agp)
    new_agp_df = new_agp_df.reset_index().set_index(['chrom', 'id'])
    chunk_df_iter = import_validpairs(args.validpairs)
   
    with mp.Pool(args.threads) as pool:
        task_list = []
        for df in chunk_df_iter:
            f = pool.apply_async(liftOverChunk, 
                    [old_agp_df, new_agp_df, df])
            task_list.append(f)
        
        results = []
        for f in task_list:
            results.append(f.get())
        
        result_df = pd.concat(results)
        result_df.to_csv('out.tsv', sep='\t', header=None, index=None)
        


if __name__ == "__main__":
    liftOverValidPairs(sys.argv[1:])
