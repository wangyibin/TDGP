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
import warnings

from collections import OrderedDict
from pandarallel import pandarallel

warnings.filterwarnings('ignore')

ValidPairsHeader = ["read", "chr1", "pos1", "strand1",
                    "chr2", "pos2", "strand2", "size",
                    "site1", "site2"]
ValidPairsHeader_dtype = {'pos1': np.int32, 'pos2': np.int32,
                        'size': np.int32}
AGP_NAMES_tig = ['chrom', 'start', 'end', 'number',
                 'type', 'id', 'tig_start', 'tig_end', 'orientation']
AGP_NAMES_tig_dtype = {'start': np.int32, 'end': np.int32, 
                'number': np.int32}
AGP_NAMES_gap = ['chrom', 'start', 'end', 'number',
                 'type', 'length', 'type', 'linkage', 'evidence']

def debug(level=logging.DEBUG):
    """
    Basic config logging format
    """
    from TDGP.apps.font import magenta, green, yellow, white
    formats = white("%(asctime)s") 
    formats += magenta(" <%(module)s:%(funcName)s>")
    formats += white(" [%(levelname)s]")
    formats += yellow(" %(message)s")
    logging.basicConfig(level=level, format=formats, datefmt="[%Y-%m-%d %H:%M:%S]")


def import_agp(agpfile, split=True):
    """
    import agp file and return a dataframe
    """
    df = pd.read_csv(agpfile, sep='\t',
                     header=None, index_col=None,
                     names=AGP_NAMES_tig, 
                     dtype=AGP_NAMES_tig_dtype)
    if split:
        tig_df = df[df['type'] == 'W']
        gap_df = df[df['type'] == 'U']
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

def import_validpairs(infile, chunksize=10000):
    """
    import validpairs which from hicpro.

    Params:
    --------
    infile: `str` input validpairs file
    chunksize: `int` chunksize for pandas reader [default: 10000].

    Returns:
    --------
    out: `pandas.DataFrame.chunk.iter`

    Examples:
    --------
    >>> import_validpairs("All.ValidPairs")
    """
    df = pd.read_csv(infile, usecols=range(10), header=None,
                    sep='\t', index_col=None, chunksize=chunksize,
                    names=ValidPairsHeader,  
                    dtype=ValidPairsHeader_dtype)
    logging.debug("Load `{}` with chunksize `{}`".format(infile, chunksize))
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
    # del old_chrom_df, new_chrom_df, old_tig_res
    # gc.collect()
    return chrom, int(new_chrom_pos)


def liftOverChunk(ns, chunk_df, idx):
    """
    liftover coordinate of validpairs in a chunk.

    Params:
    --------
    ns: multiprocessing.Manager.NameSpace with `old_agp_df`, `new_agp_df`, and `tmp`.
    chunk_df: `dataframe` dataframe of validpairs.
    idx: `int` number of chunk.

    Returns:
    --------
    chunk_df or name of results

    Examples:
    --------
    >>> liftOverChunk(old_agp_df, new_agp_df, chunk_df, 1, "tmp")
    "tmp/1.tmp.ValidPairs

    """
    
    chunk_df.chr1, chunk_df.pos1 = zip(*chunk_df.apply(lambda x: get_new_coord(x.chr1, 
                                            x.pos1, ns.old_agp_df, ns.new_agp_df), axis=1))
    chunk_df.chr2, chunk_df.pos2 = zip(*chunk_df.apply(lambda x: get_new_coord(x.chr2, 
                                            x.pos2, ns.old_agp_df, ns.new_agp_df), axis=1))
    chunk_df.dropna(inplace=True)
    chunk_df.loc[:, 'pos1'] = chunk_df['pos1'].astype('int32')
    chunk_df.loc[:, 'pos2'] = chunk_df['pos2'].astype('int32')
    if tmp:
        res = "{}/{}.tmp.ValidPairs".format(ns.tmp, idx)
        chunk_df.to_csv(res, sep='\t', header=None, index=None)
        # del chunk_df
        # gc.collect()
        logging.debug("Successful converted `{}` chunk".format(idx))
        return res
    else:
        return chunk_df

def liftOverParallel(chunk_df_iter, old_agp_df, 
            new_agp_df, tmp, threads=4):
    """
    liftover coordinate of validpairs in parallel.

    Params:
    --------
    chunk_df_iter: `dataframe` dataframe of validpairs.
    old_agp_df: `dataframe` old agp dataframe.
    new_agp_df: `dataframe` new agp dataframe.
    tmp: `str` tempoary direcory of output.
    threads: `int` threads number of program [default: 4]

    Returns:
    --------
    chunk_df or name of results

    Examples:
    --------
    >>> liftOverParallel(chunk_df_iter, "tmp", threads=4)
    ["tmp/1.tmp.ValidPairs]

    """
    results = []
    pandarallel.initialize(nb_workers=threads, verbose=1)
    for i, chunk_df in enumerate(chunk_df_iter):
        chunk_df.chr1, chunk_df.pos1 = \
            zip(*chunk_df.parallel_apply(lambda x: get_new_coord(x.chr1, 
                                        x.pos1, old_agp_df, new_agp_df), axis=1))
        chunk_df.chr2, chunk_df.pos2 = \
            zip(*chunk_df.parallel_apply(lambda x: get_new_coord(x.chr2, 
                                        x.pos2, old_agp_df, new_agp_df), axis=1))
        chunk_df.dropna(inplace=True)
        if tmp:
            res = "{}/{}.{}.ValidPairs".format(tmp, i, os.getpid())
            chunk_df.to_csv(res, sep='\t', header=None, index=None, 
                        float_format='%.0f')
        
            logging.debug("Successful converted `{}` chunk".format(i))
            results.append(res)
        else:
            results.append(chunk_df)
    
    return results

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
    pOpt.add_argument('-o', '--output', default="All.New.ValidPairs", 
            help='output file [default: All.New.ValudPairs]')
    pOpt.add_argument('-t', '--threads', type=int, default=4,
            help='number of threads [default: %(default)s]')
    pOpt.add_argument('-c', '--chunksize', type=int, default=10000,
            help='chunksize of  pandas read [default: %(default)s]')
    pOpt.add_argument('-T', dest='tmp', default='tmp', 
            help='tempoary directory of output, if set to None \n'
            'that will not split write to disk, direct write to \n'
            'a file. [default: %(default)s]')
    pOpt.add_argument('--verbose', default=2, type=int,
            choices=[0, 1, 2],
            help='The verbosity level\n'
                '0 - Don\'t display any logs\n'
                '1 - Display only warning logs\n'
                '2 - Display all logs\n')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    # verbose
    if args.verbose == 0:
        debug(level=logging.ERROR)
    elif args.verbose == 1:
        debug(level=logging.WARNING)
    elif args.verbose == 2:
        debug(level=logging.DEBUG)

    old_agp_df, _ = import_agp(args.old_agp)
    new_agp_df, _ = import_agp(args.new_agp)
    new_agp_df = new_agp_df.reset_index().set_index(['chrom', 'id'])
    
    chunk_df_iter = import_validpairs(args.validpairs, args.chunksize)
    if args.tmp:
        try:
            os.makedirs(args.tmp)
        except:
            pass
    results = liftOverParallel(chunk_df_iter, old_agp_df, new_agp_df, 
                    'tmp', args.threads)
    # with mp.Pool(args.threads, maxtasksperchild=args.threads) as pool:
    #     task_list = []
    #     for i, df in enumerate(chunk_df_iter):
    #         f = pool.apply_async(liftOverChunk, 
    #                 [ns, df, i])
    #         task_list.append(f)
        
    #     results = []
    #     for f in task_list:
    #         results.append(f.get())
    if not args.tmp:
        result_df = pd.concat(results)
        result_df.to_csv(args.output, sep='\t', header=None, index=None)
        logging.debug("Done, results is writed in `{}`".format(args.output))
    else:
        flag = os.system("cat {} > {}".format(" ".join(results), args.output))
        os.system("rm -rf {}".format(args.tmp))
        if not flag:
            logging.debug("Done, results is writed in `{}`".format(args.output))
        else:
            logging.debug("Error, something wrong happend in cat")


if __name__ == "__main__":
    liftOverValidPairs(sys.argv[1:])
