#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Some utils tools for parsing hicmatrix file (cool, hicpro ...)
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import cooler
import numpy as np 
import pandas as pd

from scipy.sparse import csr_matrix
from TDGP.apps.base import ActionDispatcher
from TDGP.apps.base import check_file_exists, debug
from TDGP.apps.base import listify


def main():

    actions = (
            ("getCisCounts", "get cis counts from hicmatrix"),
            ("getTransCounts", "get trans counts from hicmatrix")
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


class myCool(object):
    """
    my cool file parser library
    """
    def __init__(self):
        pass

def cool2matrix(cool, symmetric=True):
    """
    convert cool file to matrix.
    this function is copy from `github.com/deeptools/HiCMatrix`.
    Params:
    -------
    coolfile: `str` coolfile path
    symmetric: `bool` if return symmetric matrics

    Returns:
    --------
    hm: `numpy.array`

    Examples:
    ---------
    >>> cool = cooler.Cooler('sample.cool')
    >>> hicmatrix = cool2matrix(cool)
    """
    matrixDataFrame = cool.matrix(balance=False,
                                    sparse=False, 
                                    as_pixels=True)

    used_dtype = np.int32
    if np.iinfo(np.int32).max < cool.info['nbins']:
        used_dtype = np.int64
    count_dtype = matrixDataFrame[0]['count'].dtype
    data = np.empty(cool.info['nnz'], dtype=count_dtype)
    instances = np.empty(cool.info['nnz'], dtype=used_dtype)
    features = np.empty(cool.info['nnz'], dtype=used_dtype)
    i = 0
    size = cool.info['nbins'] // 32
    if size == 0:
        size = 1
    start_pos = 0
    while i < cool.info['nbins']:
        matrixDataFrameChunk = matrixDataFrame[i:i + size]
        _data = matrixDataFrameChunk['count'].values.astype(count_dtype)
        _instances = matrixDataFrameChunk['bin1_id'].values.astype(used_dtype)
        _features = matrixDataFrameChunk['bin2_id'].values.astype(used_dtype)

        data[start_pos:start_pos + len(_data)] = _data
        instances[start_pos:start_pos + len(_instances)] = _instances
        features[start_pos:start_pos + len(_features)] = _features
        start_pos += len(_features)
        i += size
        del _data
        del _instances
        del _features

    matrix = csr_matrix((data, (instances, features)), 
                        shape=(np.int(cool.info['nbins']), 
                                np.int(cool.info['nbins'])),
                        dtype=count_dtype)
    hm = matrix.toarray()
    ## whether obtain symmetric matrices or not up_triangle
    if symmetric:
        hm += hm.T - np.diag(hm.diagonal())

    return hm


## out command 

def getCisCounts(args):
    """
    %(prog)s <sample.cool> [Options]

        To get cis interaction count from coolfile per bins.
    
    """
    p = p=argparse.ArgumentParser(prog=getCisCounts.__name__,
                        description=getCisCounts.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('coolfile', 
            help='input file of cool')
    pOpt.add_argument('-o', '--out', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    cool = cooler.Cooler(args.coolfile)
    hm = cool2matrix(cool)
    
    counts = np.zeros(cool.info['nbins'])
    bins = cool.bins()[:].copy()
    for chrom in cool.chromnames:
        idx = cool.bins().fetch(chrom).index
        counts[idx] = hm[idx][:, idx].sum(axis=1)
    
    bins['counts'] = counts
    bins.to_csv(args.out, sep='\t', header=None, index=None)


def getTransCounts(args):
    """
    %(prog)s <sample.cool> [Options]

        To get trans interaction count from coolfile per bins.
    
    """
    p = p=argparse.ArgumentParser(prog=getTransCounts.__name__,
                        description=getTransCounts.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('coolfile', 
            help='input file of cool')
    pOpt.add_argument('-o', '--out', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    cool = cooler.Cooler(args.coolfile)
    hm = cool2matrix(cool)

    counts = np.zeros(cool.info['nbins'])
    bins = cool.bins()[:].copy()
    for chrom in cool.chromnames:
        idx = cool.bins().fetch(chrom).index
        start = idx[0]
        end = idx[-1]
        hm[start: end + 1, start: end + 1] = 0
        counts[idx] = hm[idx].sum(axis=1)

    
    bins['counts'] = counts 
    bins.to_csv(args.out, sep='\t', header=None, index=None)


if __name__ == "__main__":
    main()