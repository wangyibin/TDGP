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
        ("getTransCounts", "get trans counts from hicmatrix"),
        ("getInteractionScore", "get interaction score of haplotypes")
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
    # whether obtain symmetric matrices or not up_triangle
    if symmetric:
        hm += hm.T - np.diag(hm.diagonal())

    return hm


class CalObsExp(object):
    """
    Calculate the Observed/Expected normalization of hic matrix.

    Examples:
    --------
    >>> oe = CalObsExp(`sample_10000_iced.cool`)
    >>> oe.matrix_list 
    # [chrom1_oe_matrix, chrom2_oe_matrix ...]
    >>> oe.getInteractionScore()  ## get interaction score in each bins

    """

    def __init__(self, coolfile):

        self.coolfile = coolfile
        self.cool = cooler.Cooler(self.coolfile)
        self.matrix_list = []
        self.chromnames = self.cool.chromnames
        for c in self.chromnames:
            raw = self.cool.matrix(balance=False, sparse=False).fetch(c)
            raw[np.isnan(raw)] = 0
            expected = self.expected_matrix(raw)
            with np.errstate(invalid="ignore"):
                obs_exp = raw / expected
            obs_exp[expected == 0] = 0
            self.matrix_list.append(obs_exp)

    def expected_matrix(self, raw):
        """
        Calculate expected of hic matrix

        Params:
        --------
        raw: `matrix`  hic matrix

        Returns:
        --------
        out: `matrix` expected matrix

        Examples:
        --------
        >>> expected = self.expected_matrix(raw)

        """
        tmp = raw.sum(axis=0) != 0  # valid rows or columns
        n = raw.shape[0]
        expected = np.zeros_like(raw)
        idx = np.arange(n)
        for i in idx:
            if i > 0:
                valid = tmp[:-i] * tmp[i:]
            else:
                valid = tmp
            current = raw.diagonal(i)[valid]
            if current.size > 0:
                v = current.mean()
                if i > 0:
                    expected[idx[:-i], idx[i:]] = v
                    expected[idx[i:], idx[:-i]] = v
                else:
                    expected[idx, idx] = v
        return expected

    def getInteractionScore(self, output=""):
        """
        get the interaction score which were calculated as the 
        average of the diatance-normalized contacts in each bin.

        Params:
        --------
        None

        Returns:
        --------
        out: `pd.DataFrames` 

        Examples:
        ---------
        >>> self.getInteractionScore()
            chrom start end score
         0   Chr1 0 10000 1.2
         1   Chr1 10000 20000 1.1
        """
        self.binsDf = self.cool.bins()[:].copy(deep=True)
        try:
            self.binsDf.drop('DIs', axis=1, inplace=True)
        except KeyError:
            pass
        df_lists = []
        for matrix, chrom in zip(self.matrix_list, self.chromnames):
            tmp_df = self.binsDf[self.binsDf['chrom'] == chrom].copy(deep=True)
            idx = tmp_df.index.tolist()
            tmp_df['score'] = matrix.mean(axis=1)
            df_lists.append(tmp_df)

        df = pd.concat(df_lists, axis=0)
        if output:
            df.to_csv(output, sep='\t', header=None, index=None)
            logging.debug('Done, results is in `{}`'.format(output))

        return df

    @staticmethod
    def getHapInteractionScore(chromlists, wrkdir="./", 
                coolfile='AP85_10000_iced.cool', outdir=""):
        """
        calculate interaction score in each haplotype hic matrix

        Params:
        -------
        chromlists: `list` list of haplotype names
        wrkdir: `string` workdir [default: "./"]
        coolfile: `string` name of coolfile [default: "AP85_10000_iced.cool"]

        Returns:
        --------
        out: `pd.DataFrame` all of interaction results

        Examples:
        ---------
        >>> hap_list = ['ChrHA', 'ChrHB', 'ChrHC', 'ChrHD']
        >>> CalObsExp.getHapInteractionScore(hap_list)

        """
        res = []
        for hap in chromlists:
            cool = "{}/{}/{}".format(wrkdir, hap, coolfile)
            
            if outdir:
                output = hap + "_" + \
                    op.basename(coolfile).rsplit(".", 1)[
                        0] + "_interactionScore.bg"
            else:
                output = ""
            hap_df = CalObsExp(cool).getInteractionScore(output)
        
            res.append(hap_df)

        df = pd.concat(res, axis=0)
        df = df.reindex()
        return df



# out command

def getInteractionScore(args):
    """
    %(prog)s <sample_10000_iced.cool> [Options]

        Calculate the interaction scores base O/E matrix
    """
    p = argparse.ArgumentParser(prog=getInteractionScore.__name__,
                        description=getInteractionScore.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    #pReq.add_argument('coolfile', 
    #        help='coolfile of hicmatrix')
    pOpt.add_argument('--haplists', default=['ChrHA', 'ChrHB', 'ChrHC', 'ChrHD'],
            nargs="*", help='haplotype names [default: %(default)s]')
    pOpt.add_argument('--wrkdir', default='./', 
            help='directory of haplotypes [default: %(default)s] ')
    pOpt.add_argument('--coolfile', default='AP85_10000_iced.cool',
            help='name of coolfile [default: %(default)s]')
    pOpt.add_argument('-o', '--outdir', default='',
            help='directory of output each haplotype results [default: None]')
    pOpt.add_argument('-o', '--out', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    df = CalObsExp.getHapInteractionScore(args.haplists, 
                                        args.wrkdir, 
                                        args.coolfile,
                                        outdir=args.outdir)
    
    df.to_csv(args.out, sep='\t', header=None, index=False)
    logging.debug('Successful, results is in `{}`'.format(args.out.name))

    


def getCisCounts(args):
    """
    %(prog)s <sample.cool> [Options]

        To get cis interaction count from coolfile per bins.

    """
    p = p = argparse.ArgumentParser(prog=getCisCounts.__name__,
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
    logging.debug(
        "Successful outputting cis counts bedgraph to `{}`".format(args.out.name))

def getTransCounts(args):
    """
    %(prog)s <sample.cool> [Options]

        To get trans interaction count from coolfile per bins.

    """
    p = p = argparse.ArgumentParser(prog=getTransCounts.__name__,
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
    logging.debug("Successful outputting trans counts bedgraph to `{}`".format(args.out.name))


if __name__ == "__main__":
    main()
