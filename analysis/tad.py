#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
TADs analysis libraries.
"""
from __future__ import print_function

import logging
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import os.path as op
import os
import sys

from collections import OrderedDict, defaultdict
from intervaltree import Interval, IntervalTree
from optparse import OptionParser
from TDGP.apps.base import ActionDispatcher
from TDGP.apps.base import debug, listify, check_file_exists
from TDGP.apps.base import BaseFile, Line
from TDGP.apps.utilities import chrRangeID


debug()

def main():

    actions = (
        ('plotSizeDist', 'plot the tad size distribution a list of samples'),
        ('getBottom', 'get bottom tads from hitad results'),
        ('stat', 'stat TADs informations'),
        ('getBoundaryBed', 'get tad boundary'), 
        ('testPipe', 'testPipe'),
        ('annotate', 'annotate tad'),
        ('whichTAD', 'find gene location in TADs')
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

class TADLine(Line):
    """
    The object of TAD Line, `chrom start end [level]`.

    Params:
    ---------
    line: `str` line of TAD.

    Returns:
    ---------
    out: `class`

    Examples:
    ----------
    >>> line = "Chr1    0    10000    2"
    >>> tl = TADLine(line)
    >>> tl.chrom
    Chr1
    >>> tl.interval
    Interval(0, 10000, "Chr1")
    """
    def __init__(self, line):
        super(TADLine, self).__init__(line)
        self.chrom, self.start, self.end = self.line_list[:3]
        self.start, self.end = int(self.start), int(self.end)
        self.interval = Interval(self.start, self.end, self.chrom)

        if len(self.line_list) > 3:
            self.level = int(self.line_list[3])


class TADFile(BaseFile):
    """
    The object of TAD File.

    Params:
    --------
    infile: `str`

    Returns:
    ---------
    out: `class`

    Examples:
    ---------
    >>> tdf = TADFile('sample.tad.txt')
    >>> tdf
    """
    def __init__(self, infile):
        super(TADFile, self).__init__(infile)
        self.infile = infile
        self.getTAD()
        self.getBottomDict()
        self.getBottom()
        self.getBottomSizes()
        self.getSizes()

    def getLine(self):
        with open(self.infile) as fp:
            for line in fp:
                yield TADLine(line)

    def __iter__(self):
        pass

    def getTAD(self):
        self.TADDict = defaultdict(list)
        for tl in self.getLine():
            self.TADDict[tl.chrom].append(
                Interval(tl.start, tl.end, tl.chrom))
    
    def getSizes(self):
        """
        Get all tads sizes
        
        Returns:
        ---------
        out: `list`
        
        Examples:
        ---------
        >>> tf.getSizes()
        >>> tf.sizes
        [1000, 20000, ...]
        """
        self.sizes = []
        for value in self.TADDict.values():
            for interval in value:
                self.sizes.append(interval.length())
        
        return self.sizes

    def getBottomDict(self):
        self.bottomDict = defaultdict(lambda :IntervalTree())
        for chrom in self.TADDict:
            for interval in self.TADDict[chrom]:
                overlaps = list(self.bottomDict[chrom].overlap(
                    interval.begin, interval.end))
                if overlaps:
                    for overlap in overlaps:
                        if (interval.length()) < overlap.length():
                            self.bottomDict[chrom].remove(overlap)
                self.bottomDict[chrom].add(interval)

    def getBottom(self):
        self.bottom = []
        for chrom in self.bottomDict:
            for interval in sorted(self.bottomDict[chrom]):
                self.bottom.append((chrom, interval.begin,
                                   interval.end))

    def getBottomSizes(self):
        """
        Get a list of bottom sizes.
        """
        self.bottomSizes = []
        for item in self.bottom:
            size = item[2] - item[1]
            self.bottomSizes.append(size)
    
    def getBoundaryDict(self):
        """
        Get the TADs boundary
        Params:
        -------
        None

        Returns:
        --------
        out: `OrderedDict`

        Examples:
        --------
        >>> tf = TADFile('sample.txt')
        >>> tf.getBottomDict()
        >>> tf.getBoundaryDict()
        >>> tf.boundaryDict
        OrderedDict(('Chr1', {1000, 20000 ...} ...))
        """
        self.boundaryDict = OrderedDict()
        for chrom in self.bottomDict:
            self.boundaryDict[chrom] = set()
            for interval in self.bottomDict[chrom]:
                self.boundaryDict[chrom].add(interval.begin)
                self.boundaryDict[chrom].add(interval.end)

        return self.boundaryDict


    def getBoundary(self):
        """
        list of TADs boundary

        Returns:
        --------
        out: `list`
        """
        self.boundary = []
        for chrom in self.BoundaryDict:
            for boundary in sorted(self.BoundaryDict[chrom]):
                self.boundary.append((chrom, boundary))
        return self.boundary


    @classmethod
    def getBoundaryBed(self, boundaryDict, chromSize, 
            updistance=2e4, downdistance=2e4):
        """
        Get the TADs boundary bed with .

        Params:
        --------
        boundaryDict: `dict` the dict of boundary
        chromSize: `dict` the dict of chrom sizes
        [Options]
        updistance: `int` the upstream distance of boundary[default: 2e4]
        downdistance: `int` the downstream distance of boundary[default: 2e4]

        Returns:
        --------
        out: `list` a list of boundary bed

        Examples:
        --------
        >>>tf.getBoundaryBed(boundaryDict, chromSize)

        """
        self.boundaryBed = []
        for chrom in boundaryDict:
            for boundary in boundaryDict[chrom]:
                upstream = boundary - updistance\
                    if (boundary - updistance) >= 0 \
                    else 0
                downstream = boundary + downdistance \
                    if (boundary + downdistance <= chromSize[chrom]) \
                    else chromSize[chrom]

                res = (chrom, upstream, downstream, boundary)
                self.boundaryBed.append(res)

        return self.boundaryBed

    def getTADbyLevel(self):
        """
        Get TAD Dict use the level as the keys.
        """
        self.LevelDict = defaultdict(list)
        for tl in self.getLine():
            level = tl.level
            self.LevelDict[level].append(tl.interval)

        return self.LevelDict

    def getSizeDictPerLevel(self):
        self.sizeDict = OrderedDict()
        for level in sorted(self.LevelDict.keys()):
            for interval in self.LevelDict[level]:
                if level not in self.sizeDict:
                    self.sizeDict[level] = []
                self.sizeDict[level].append(interval.length())


    def plotSizeDistPerLevel(self, ax, out, exclude=[], scale=1000, xmin=0,
                     xmax=1000, step=200):
        scale_units = {1: 'bp', 1000: 'kb', 1e6: 'Mb'}


        for level in self.sizeDict:
            data = np.array(self.sizeDict[level]) / scale
            if level in exclude:
                continue
            sns.distplot(data, hist=False, kde=True, ax=ax,
                         label="level %d (%d)"%(level, len(data)))

        ax.set_xlim(xmin, xmax)
        ax.set_xticks(range(xmin, xmax, step))
        ax.set_xlabel('TAD Size ({})'.format(scale_units[scale]))
        ax.set_ylabel('Frequency')
        ax.set_title('TAD Size Distribution Per Level')
        plt.savefig(out, dpi=300)



    @classmethod
    def plotSizeDist(self,  ax, data, out, label='Sample', scale=1000,
                     xmin=0, xmax=800, step=100):

        """
        Plot
        """
        scale_units = {1: 'bp', 1000: 'kb', 1e6: 'Mb'}

        data = np.array(data) / scale
        sns.distplot(data, hist=False, kde=True, ax=ax,
                     label="{} ({})".format(label, len(data)))
        ax.set_xlim(xmin, xmax)
        ax.set_xticks(range(xmin, xmax + 1, step))
        ax.set_xlabel('TAD Size ({})'.format(scale_units[scale]))
        ax.set_ylabel('Frequency')
        ax.set_title('TAD Size Distributions')
        plt.savefig(out, dpi=300)

    def plotSizeDistMulti(self, ax, data, label, scale=1000):
        scale_units = {1: 'bp', 1000: 'kb', 1e6: 'Mb'}
        data = np.array(data) / scale
        ax = sns.distplot(data, hist=False, kde=True, ax=ax,
                     label="{} ({})".format(label, len(data)))

        return ax


class TADConservedLine(Line):
    """
    Object of TADConservedLine.

    Params:
    --------

    Returns:
    --------
    
    Examples:
    --------
    """
    def __init__(self, line):
        super(TADConservedLine, self).__init__(line)
        self.chrom, self.start, self.end = self.line_list[:3]
        self.genes = self.line_list[3]
        self.start, self.end = int(self.start), int(self.end)

    @classmethod
    def from_list(self, input):
        """
        function of get object from a list.
        list must contain `chrom`, `start`, `end`, `gene1,gene2`
        """
        assert isinstance(input, list), "input is not a list"

        self.line_list = input
        self.chrom, self.start, self.end = self.line_list[:3]
        self.genes = self.line_list[3]
        self.start, self.end = int(self.start), int(self.end)
    



class TADConserved(object):
    """
    Object of TAD conservation analysis.

    Params:
    --------

    Returns:
    --------

    Examples:
    --------
    
    """
    def __init__(self):
        pass
    
    def fetchSyntenyGene(self, tad, gene_tree):
        if not isinstance(tad, Interval):
            tad = Interval
        result = gene_tree.overlap(*tad)
        sorted_result = sorted(result)

        return result
    
    @staticmethod
    def getGene(tads, genes, fraction=0.7):
        check_file_exists(tads)
        check_file_exists(genes)
        if 0 > float(fraction) > 1:
            logging.error('The option `-F` must set in '
                'range [0, 1], and you set {}'.format(fraction))
            sys.exit()

        bedtools_cmd = "bedtools intersect -a {} -b {} -wao -F {} |"
                        " cut -f 1-3,7 ".format(tads, genes, fraction)
        db = OrderedDict()
            
        for line in os.popen(bedtools_cmd):
            line_list = line.strip().split()
            ID = chrRangeID(line_list[:3])
            gene = line_list[3]
            if ID not in db:
                db[ID] = []
            db[ID].append(gene)
        
        for ID in db:
            gene_list = db[ID]
        
        return db

### outsite command start ###
def getBottom(args):
    """
    %prog getBottom <tad.txt> [Options]

    To get the bottom tads of hitad results.
    """
    p = OptionParser(getBottom.__doc__)

    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(p.print_help())

    tad, = args
    tf = TADFile(tad)
    for bottom in tf.bottom:
        print("\t".join(map(str, bottom)))


def getBoundaryBed(args):
    """
    %prog <tad.bed> <chrom.sizes> [Options]

    get a bed file of the tad boundary.
    """
    p = OptionParser(getBoundaryBed.__doc__)
    p.add_option('-a', '--up', type=int, default=250,
            help='the upstrean distance of boundary '
            '[default: %default]')
    p.add_option('-b', '--down', type=int, default=250,
            help='the downstream distance of boundary '
            '[default: %default]')
    
    opts, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(p.print_help())

    tadFile, chromSize = args
    up, down = opts.up, opts.down
    if not op.exists(tadFile) or not \
            op.exists(chromSize):
        logging.error('The input file is not exists')
    chrom_dict = dict(i.strip().split() \
            for i in open(chromSize) if i.strip())

    tf = TADFile(tadFile)
    tf.getBoundaryDict()
    boundaryBed = tf.getBoundaryBed(tf.boundaryDict, chrom_dict, 
            up, down)
    for item in sorted(boundaryBed):
        print("\t".join(map(str, item[:3])))





def plotSizeDist(args):
    """
    %prog plotSizeDist <tad1.bed> [tad2.bed ...] [Options]

    Given some tad bed file to plot their sizes distributions.
    """
    scale_units = {1: 'bp', 1000: 'kb', '1e6': 'Mb'}
    p = OptionParser(plotSizeDist.__doc__)
    p.add_option('-o', '--out', default='tad_sizes_dist.pdf',
                help='out of plot [default: %default]')
    p.add_option('--all', default=False, action='store_true',
                help='plot all levels of tads [default: %default]')
    p.add_option('-s', '--scale', default=1000, type=int,
                help='the scale of xticks [default: %default]')
    p.add_option('--xmin', default=0, type=int,
                help='min value of xticks [default: %default]')
    p.add_option('--xmax', default=800, type=int,
                help='max value of xticks [default: %default]')
    p.add_option('--step', default=100, type=int,
                help='the step of xticks [default: %default]')

    opts, args = p.parse_args(args)
    if len(args) < 1:
        sys.exit(p.print_help())
    out, scale, xmin, xmax, step = opts.out, opts.scale, \
                                    opts.xmin, opts.xmax, \
                                    opts.step
    
    logging.debug('Plotting ...')
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    for tad in args:
        label = tad.rsplit(".")[0]
        tf = TADFile(tad)
        data = tf.Sizes if opts.all else tf.bottomSizes
        ax = tf.plotSizeDistMulti(ax, data, label=label,
                                 scale=scale)
    ax.set_xlim(xmin, xmax)
    ax.set_xticks(range(xmin, xmax + 1, step))
    ax.set_xlabel('TAD Size ({})'.format(scale_units[scale]))
    ax.set_ylabel('Frequency')
    ax.set_title('TAD Size Distributions')
    plt.savefig(out, dpi=300, bbox_inches='tight')
    logging.debug('Success file is save as {}'.format(out))


def stat(args):
    """
    tads informations stat
        total_num total_size genome_size percentage
    """
    p = OptionParser(stat.__doc__)
    p.add_option('-g','--genome', type=int,
            help='the genome size of species')
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())
    if not opts.genome:
        logging.error('Must input genome size '
                'with `-g` option')
        sys.exit()
    tad, = args
    tf = TADFile(tad)
    total_num, total_size = 0, 0
    for size in tf.bottomSizes:
        total_num += 1
        total_size += size

    print('#Total number\tTotal size\tGenome size\tPercentage')
    print("{}\t{}\t{}\t{:.2%}".format(total_num, total_size, 
        opts.genome, total_size * 1.0 / opts.genome))



def testPipe(args):
    """
    test pipe for tad annotate.
    """
    db = OrderedDict()
    if not sys.stdin.isatty():
        handle = sys.stdin
    else:
        pass

    for line in sys.stdin:
        line_list = line.strip().split()
        ID = chrRangeID(line_list[:3])
        gene = line_list[3]
        if ID not in db:
            db[ID] = []
        db[ID].append(gene)
    
    for ID in db:
        print(ID + "\t" + ",".join(db[ID]))


def annotate(args):
    """
    %prog tad.bed gene.bed [Options]
    Annotate tads with gene.
    """

    p = OptionParser(annotate.__doc__)
    p.add_option('-F', dest='fraction', default='0.7', 
                    help='the fraction of gene overlap of tads')

    opts, args = p.parse_args(args)

    if len(args) != 2 :
        sys.exit(p.print_help())
    
    tads, genes = args
    fraction = opts.fraction
    check_file_exists(tads)
    check_file_exists(genes)
    if 0 > float(fraction) > 1:
        logging.error('The option `-F` must set in '
            'range [0, 1], and you set {}'.format(fraction))
        sys.exit()

    bedtools_cmd = "bedtools intersect -a {} -b {} -wao -F {} | cut -f 1-3,7 ".format(
        tads, genes, fraction)
    db = OrderedDict()
        
    for line in os.popen(bedtools_cmd):
        line_list = line.strip().split()
        ID = chrRangeID(line_list[:3])
        gene = line_list[3]
        if ID not in db:
            db[ID] = []
        db[ID].append(gene)
    
    for ID in db:
        gene_list = db[ID]
        length = len(gene_list) if "." not in gene_list else 0
        print("\t".join(chrRangeID(ID, axis=1)) + "\t" + \
            ",".join(gene_list) + "\t" + \
            str(length), file=sys.stdout)
    


def whichTAD(args):
    """
    %prog gene.bed tad.bed Options
    
    find gene location in tads
    """
    p = OptionParser(annotate.__doc__)
    p.add_option('-f', dest='fraction', default='0.7', 
                    help='the fraction of gene overlap of tads')

    opts, args = p.parse_args(args)

    if len(args) != 2 :
        sys.exit(p.print_help())
    
    genes, tads = args
    fraction = opts.fraction
    check_file_exists(tads)
    check_file_exists(genes)
    if 0 > float(fraction) > 1:
        logging.error('The option `-f` must set in '
            'range [0, 1], and you set {}'.format(fraction))
        sys.exit()

    bedtools_cmd = "bedtools intersect -a {} -b {} -wao -f {} | cut -f 4-8 ".format(
        genes, tads, fraction)
    for line in os.popen(bedtools_cmd):
        print(line.strip())
    """
    db = OrderedDict()
        
    for line in os.popen(bedtools_cmd):
        line_list = line.strip().split()
        gene = line_list[0]
        ID = chrRangeID(line_list[1:4])
        
        if ID not in db:
            db[ID] = []
        db[ID].append(gene)
    """
    """
    for ID in db:
        gene_list = db[ID]
        length = len(gene_list) if "." not in gene_list else 0
        print("\t".join(chrRangeID(ID, axis=1)) + "\t" + \
            ",".join(gene_list) + "\t" + \
            str(length), file=sys.stdout)
    """
if __name__ == "__main__":
    main()
