#!/usr/bin/env python
# -*- coding:utf-8 -*-


import logging
import os
import os.path as op

from collections import OrderedDict
from intervaltree import Interval, IntervalTree
from TDGP.apps.base import BaseFile, Line

from TDGP.apps.base import debug

debug()


class LinkLine(Line):
    """
    Object of anchor file line.
    """

    def __init__(self, line):
        super(LinkLine, self).__init__(line)

        assert len(self.line_list) >= 6, \
            "This Line is not a normal linkfile line. Please check file."
        (self.chrom1, self.start1, self.end1,
         self.chrom2, self.start2, self.end2) = self.line_list[:6]

        self.start1, self.end1, self.start2, self.end2 = list(
            map(int, (self.start1, self.end1, self.start2, self.end2)))

        if  6 < len(self.line_list) <= 8:
            self.gene1, self.gene2 = self.line_list[6: 9]


    def length1(self):
        return self.end1 - self.start1 + 1

    def length2(self):
        return self.end2 - self.start2 + 1

    def interval1(self):
        return Interval(self.start1, self.end1, self.gene1)

    def interval2(self):
        return  Interval(self.start2, self.end2, self.gene2)

    def reverse(self):
        self.chrom1, self.chrom2 = self.chrom2, self.chrom1
        self.start1, self.start2 = self.start2, self.start1
        self.end1, self.end2 = self.end2, self.end1
        self.gene2, self.gene2 = self.gene2, self.gene1
        self.interval1, self.interval2 = self.interval2, self.interval1
        self.length1, self.length2 = self.length2, self.length1
        self.line_list = list(map(str, [self.chrom1, self.start1, self.end1,
                         self.chrom2, self.start2, self.end2,
                         self.gene1, self.gene2]))


class LinkFile(BaseFile):
    """
    Object of linkfile.
    """
    def __init__(self, infile, axis=0):
        super(LinkFile, self).__init__(infile)
        self.infile = infile
        self.getLinkDict(axis)

    def getLine(self):
        with open(self.infile) as fp:
            for line in fp:
                yield LinkLine(line)

    def __iter__(self):
        for ll in self.getLine():
            yield ll

    def getLinkDict(self, axis=0):
        """
        get the dictionary of link.

        Params:
        --------
        axis: `int` the axis of dict key.
                if axis=0 key is chrom1
                else axis=1 key is chrom2

        Returns:
        ---------
        out: `OrderedDict` of link

        Examples:
        ---------
        >>> lf = LinkFile(sample.link)
        >>> lf.LinkDict


        """
        assert axis == 0 or axis == 1, \
            "The axis must bewtein 0 and 1."
        self.LinkDict = OrderedDict()
        for ll in self.getLine():
            if axis == 1:
                ll.reverse()

            if ll.chrom1 not in self.LinkDict:
                self.LinkDict[ll.chrom1] = IntervalTree()
            self.LinkDict[ll.chrom1].addi(ll.start1, ll.end1, ll)

        return self.LinkDict

