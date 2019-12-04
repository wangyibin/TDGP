#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import os
import os.path as op
import sys

from collections import OrderedDict, defaultdict
from intervaltree import Interval, IntervalTree
from TDGP.apps.base import BaseFile
from TDGP.apps.base import debug


class BedGraphLine(object):
    """
    The object of bedGraph line

    Params:
    ---------
    infile: str

    Returns:
    ---------
    out: object of BedGraphLine

    Examples:
    ---------
    >>> bgl = BedGraphLine(line)
    >>> bgl.chrom
    Chr1
    >>> bgl.score 
    0.520
    """
    def __init__(self, line):
        self.line = line
        self.line_list = line.strip().split()

        self.chrom, self.start, self.end, \
        self.score = self.line_list
        self.start, self.end = int(self.start), int(self.end)
        self.score = float(self.score)
    
    @property
    def range(self):
        return (self.chrom, self.start, self.end, self.score)

    def __str__(self):
        return self.line

    __repr__ = __str__


class BedGraph(BaseFile):
    """
    The object of bedGraph file.

    Params:
    --------
    infile: str

    Returns:
    --------
    out: BedGraph

    Examples:
    ---------
    >>> bg = BedGraph('sample.bg')
    >>> bg.get()

    """

    def __init__(self, infile):
        super(BedGraph, self).__init__(infile)

        self.infile = infile
        self.get()
        self.getChromSizes()
        self.totalChromSize = sum(self.chromSizes.values())
        self.chromList = list(self.bedGraphDict.keys())
        _binSize, = list(self.bedGraphDict.values())[0][0]
        self.binSize = _binSize.length()
        
    def getLine(self):
        with open(self.infile) as fp:
            for line in fp:
                yield BedGraphLine(line)

    def __iter__(self):
        for bgl in self.getLine():
            yield bgl


    def get(self):
        self.bedGraphDict = OrderedDict()
        for bgl in self.getLine():
            if bgl.chrom not in self.bedGraphDict:
                self.bedGraphDict[bgl.chrom] = IntervalTree()

            self.bedGraphDict[bgl.chrom].add(Interval(bgl.start, 
                bgl.end, bgl.score))
        return self.bedGraphDict

        
    def getChromSizes(self):
        self.chromSizes = OrderedDict()
        for chrom in self.bedGraphDict:
            if chrom not in self.chromSizes:
                self.chromSizes[chrom] = 0
            for item in self.bedGraphDict[chrom]:

                if self.chromSizes[chrom] < item[1]:
                    self.chromSizes[chrom] = item[1]
        
        return self.chromSizes
    

    def getChromBinRanges(self, chrom):
        binRanges = []
        for interval_object in sorted(self.bedGraphDict[chrom]):
            binRanges.append((interval_object.begin, interval_object.end))

        return binRanges
    
    
    def getBinRanges(self):
        self.binRanges = OrderedDict()
        for chrom in self.bedGraphDict:
            self.binRanges[chrom] = self.getChromBinRanges(chrom)

        return self.binRanges 
    

