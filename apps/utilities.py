#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
utility libraries.
"""

from __future__ import print_function

import cooler
import os
import os.path as op
import sys

import numpy as np

from collections import defaultdict


def chrRangeID(args, axis=0):
    """
    Chrom range transformation.
    Examples:
    --------
    >>> args = ["Chr1", 100, 200]
    >>> chrRangeID(args)
    "Chr1:100-200"
    >>> args = "Chr1:100-200"
    >>> chrRangeID(args, axis=1)
    ("Chr1", "100", "200")
    """
    if axis == 0:
        chrom, start, end = map(str, args)
        return "{}:{}-{}".format(chrom, start, end)
    elif axis == 1:
        chrom, ranges = args.split(':')
        start, end = ranges.split('-')
        return chrom, start, end
    else:
        return 


def dictCounter(inDict):
    """
    To count dict list value num.
    Examples:
    -------
    >>> d = {1: [1,3], 2: [1]}
    >>> dictCounter(d)
    {1: 2, 2: 1}
    """
    for item in inDict:
        inDict[item] = len(inDict[item])
    return inDict



def isCooler(filename):
    """
    judge a file if a cool
    """
    try:
        c = cooler.Cooler(filename)
    except IOError:
        return False
    else:
        if isinstance(c, cooler.Cooler):
            return True
        else:
            return False


def chrom_size_convert(size):
    """
    Convert the unit of chromosome size to suitable unit.
    >>> chrom_size_convert(100000)
    100 Kbp
    >>> chrom_size_convert(1000000)
    1 Mbp
    """
    if size <= 1e3:
        label = "{:,.0f}".format((size)) + " bp"
    elif size <= 4e5:
        label = "{:,.0f}".format((size / 1e3)) + " Kbp"
    else:
        label = "{:,.1f}".format((size / 1e6)) + " Mbp"
    
    return label


def chrom_ticks_convert(ticks):
    """
    Convert a list of  chromosome size to suitable unit.
    >>> ticks = [10000, 20000, 30000]
    >>> chrom_ticks_convert(ticks)
    ['10', '20', '30Kbp']
    """
    if ticks[-1]  - ticks[1] <= 1e3:
        labels = ["{:,.0f}".format((x)) 
                  for x in ticks] 
        labels[-1] += " bp"
    elif ticks[-1]  - ticks[1] <= 4e5:
        labels = ["{:,.0f}".format((x / 1e3)) 
                  for x in ticks]
        labels[-1] += 'Kbp'
    else:
        labels = ["{:,.1f}".format((x / 1e6)) 
                  for x in ticks]
        labels[-1] += " Mbp"
    
    return labels

def chrom_ticks_convert(ticks):
    """
    Convert a list of  chromosome size to suitable unit.
    >>> ticks = [10000, 20000, 30000]
    >>> chrom_ticks_convert(ticks)
    ['10', '20', '30Kbp']
    """
    if ticks[-1]  - ticks[1] <= 1e3:
        labels = ["{:,.0f}".format((x)) 
                  for x in ticks] 
        labels[-1] += " bp"
    elif ticks[-1]  - ticks[1] <= 4e5:
        labels = ["{:,.0f}".format((x / 1e3)) 
                  for x in ticks]
        labels[-1] += 'Kbp'
    else:
        labels = ["{:,.1f}".format((x / 1e6)) 
                  for x in ticks]
        labels[-1] += " Mbp"
    
    return labels


def makeChromWindows(chrom_size, window=5e6):
    """
    To segment chromosome interval into certain windows.
    """

    chrom_size_db = dict(i.strip().split() 
                for i in open(chrom_size) if i.strip())
    
    interval_db = defaultdict(list)
    window = int(window)
    for chrom in chrom_size_db:
        size = int(chrom_size_db[chrom])
        for i in range(0, size, window):
            interval_db[chrom].append([i, i+window])
    
        interval_db[chrom][-1][1] = size

    return interval_db



def wi_test(data1, data2):
    """
    Wilcoxon rank-sum tests
    return: pvalue
    """
    from scipy import stats
    wi = stats.ranksums(data1, data2)
    return wi.pvalue


def bezier(o, c1, c2, d, step=0.05):
    """
    Return cubic beizer curve points array: 
    <http://www.moshplant.com/direct-or/bezier/math.html>
    o: origin, c1, c2: control, d: detination

    Returns:
    --------
    out: `array` xt, yt

    Examples:
    --------
    >>> beizer((0 ,4,), ( 2.5, 4), (2.5, 0), (5, 0))
    (array([0.      , 0.356875, 0.68    , 0.973125, 1.24    , 1.484375,
        1.71    , 1.920625, 2.12    , 2.311875, 2.5     , 2.688125,
        2.88    , 3.079375, 3.29    , 3.515625, 3.76    , 4.026875,
        4.32    , 4.643125, 5.      ]),
    array([4.   , 3.971, 3.888, 3.757, 3.584, 3.375, 3.136, 2.873, 2.592,
        2.299, 2.   , 1.701, 1.408, 1.127, 0.864, 0.625, 0.416, 0.243,
        0.112, 0.029, 0.   ]))

    """
    t = np.arange(0, 1 + step, step)
    pts = (o, c1, c2, d)
    px, py = zip(*pts)
    def get_array(pts):
        # calculation
        o, c1, c2, d = pts
        c = 3 * (c1 - o)
        b = 3 * (c2 - c1) - c
        a = d - o - c - b

        tsquared = t ** 2
        tcubic = tsquared * t

        return a * tcubic + b * tsquared + c * t + o
    
    return get_array(px), get_array(py)
