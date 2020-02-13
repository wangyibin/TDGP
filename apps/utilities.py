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
