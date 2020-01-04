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