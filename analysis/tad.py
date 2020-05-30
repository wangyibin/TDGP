#!/usr/bin/env python
# -*- coding:utf-8 -*-
"""
TADs analysis libraries.
"""
from __future__ import print_function

import argparse 
import cooler
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

    actions = (('plotSizeDist',
                'plot the tad size distribution a list of samples'),
               ('getBottom', 'get bottom tads from hitad results'),
               ('stat', 'stat TADs informations'), ('getBoundaryBed',
                                                    'get tad boundary'),
               ('testPipe', 'testPipe'), ('annotate', 'annotate tad'),
               ('whichTAD',
                'find gene location in TADs'), ('getSyntenicTADs',
                                                'get syntenic tads'),
               ('plotBoundary',
                'plot omics data density in boundary'), ('test', 'test'),
                ('quickPlotTAD',
                'quick plot picture to tads result visualization'),)

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
            self.TADDict[tl.chrom].append(Interval(tl.start, tl.end, tl.chrom))

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
        self.bottomDict = defaultdict(lambda: IntervalTree())
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
                self.bottom.append((chrom, interval.begin, interval.end))

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
    def getBoundaryBed(self,
                       boundaryDict,
                       chromSize,
                       updistance=2e4,
                       downdistance=2e4):
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
                    if (boundary + downdistance <= int(chromSize[chrom])) \
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

    def plotSizeDistPerLevel(self,
                             ax,
                             out,
                             exclude=[],
                             scale=1000,
                             xmin=0,
                             xmax=1000,
                             step=200):
        scale_units = {1: 'bp', 1000: 'kb', 1e6: 'Mb'}

        for level in self.sizeDict:
            data = np.array(self.sizeDict[level]) / scale
            if level in exclude:
                continue
            sns.distplot(data,
                         hist=False,
                         kde=True,
                         ax=ax,
                         label="level %d (%d)" % (level, len(data)))

        ax.set_xlim(xmin, xmax)
        ax.set_xticks(range(xmin, xmax, step))
        ax.set_xlabel('TAD Size ({})'.format(scale_units[scale]))
        ax.set_ylabel('Frequency')
        ax.set_title('TAD Size Distribution Per Level')
        plt.savefig(out, dpi=300)
        plt.savefig(out.rsplit(".")[0] + ".png", dpi=300)


    @classmethod
    def plotSizeDist(self,
                     ax,
                     data,
                     out,
                     label='Sample',
                     scale=1000,
                     xmin=0,
                     xmax=800,
                     step=100):
        """
        Plot
        """
        scale_units = {1: 'bp', 1000: 'kb', 1e6: 'Mb'}

        data = np.array(data) / scale
        sns.distplot(data,
                     hist=False,
                     kde=True,
                     ax=ax,
                     label="{} ({})".format(label, len(data)))
        ax.set_xlim(xmin, xmax)
        ax.tick_params(labelsize=12)
        ax.set_xticks(range(xmin, xmax + 1, step))
        
        ax.set_xlabel('TAD Size ({})'.format(scale_units[scale]), fontsize=14)
        ax.set_ylabel('Frequency', fontsize=14)
        ax.set_title('TAD Size Distributions', fontsize=16)
        plt.savefig(out, dpi=300)
        plt.savefig(out.rsplit(".")[0] + ".png", dpi=300)
        
    def plotSizeDistMulti(self, ax, data, label, scale=1000):
        scale_units = {1: 'bp', 1000: 'kb', 1e6: 'Mb'}
        data = np.array(data) / scale
        ax = sns.distplot(data,
                          hist=False,
                          kde=True,
                          ax=ax,
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
    def getGene(tads, genes, fraction=0.7, isnum=False, isPlot=False):
        """
        Annotate tads with genes and return as dict.

        Params:
        --------
        tads: `str` bed3 file of tad
        genes: `str` bed4 file of gene
        fraction: `str` or `float` fraction of gene 
                    overlap with tads [default: 0.7]
        isnum: `bool` if set output the gene number instead 
                    of gene list. [default: False]
        isPlot: `bool` if plot the gene number per TADs 
                    distribution. [default: False]

        Returns:
        --------
        out: `dict` dictionary of TADs annotation

        Examples:
        --------
        >>> db = TADConserved().getGene("tad.bed", "gene.bed")
        """
        check_file_exists(tads)
        check_file_exists(genes)
        if 0 > float(fraction) > 1:
            logging.error('The option `-F` must set in '
                          'range [0, 1], and you set {}'.format(fraction))
            sys.exit()

        bedtools_cmd = "bedtools intersect -a {} -b {} -wao -F {} | \
                         cut -f 1-3,7 ".format(tads, genes, fraction)
        db = OrderedDict()

        for line in os.popen(bedtools_cmd):
            line_list = line.strip().split()
            ID = chrRangeID(line_list[:3])
            gene = line_list[3]
            if ID not in db:
                db[ID] = set()
            db[ID].add(gene)

        if isnum:
            for ID in db:
                db[ID] = len(db[ID])

        if isPlot:
            assert isnum, 'isnum must specify as True'
            fig, ax = plt.subplots(figsize=(5, 5))
            sns.distplot(db.values(), hist=False, kde=True, ax=ax)
            ax.set_xticks(range(0, 41, 5))
            ax.set_xlim(0, 40)
            ax.set_xlabel('Gene number')
            ax.set_ylabel('Frequence')
            ax.set_title('Gene Number Distribution ({:,})'.format(
                sum(db.values())))
            plt.savefig('{}.gene_num_dist.pdf'.format(genes.rsplit('.', 1)[0]),
                        dpi=300)
            logging.debug('Successful to plot gene number distribution '
                          '`{}.gene_num_dist.pdf`.'.format(
                              genes.rsplit('.', 1)[0]))

        return db

    @staticmethod
    def genePairTAD(genes, tads, fraction=0.7):
        """
        Annotate tads with genes and return as dict.

        Params:
        --------
        genes: `str` bed4 file of gene
        tads: `str` bed3 file of tad
        fraction: `str` or `float` fraction of gene 
                    overlap with tads [default: 0.7]
        
        Returns:
        --------
        out: `dict` dictionary of TADs annotation

        Examples:
        --------
        >>> db = TADConserved().getGene("gene.bed", "tad.bed")

        """

        check_file_exists(tads)
        check_file_exists(genes)
        if 0 > float(fraction) > 1:
            logging.error('The option `-f` must set in '
                          'range [0, 1], and you set {}'.format(fraction))
            sys.exit()

        bedtools_cmd = "bedtools intersect -a {} -b {} -wao -f {} | cut -f 4-7 ".format(
            genes, tads, fraction)
        db = OrderedDict()
        for line in os.popen(bedtools_cmd):
            line_list = line.strip().split()
            gene, chrom, start, end = line_list
            ID = chrRangeID([chrom, start, end]) if chrom != "." \
                else "."
            if ID == ".":
                continue

            db[gene] = ID

        return db

    @staticmethod
    def getConserved(tad1,
                     tad2,
                     syngene1,
                     syngene2,
                     gene1,
                     gene2,
                     anchors,
                     fraction=0.7,
                     threshold=0,
                     gene_num=0,
                     synthre=0):
        """
        Get all syntenic TADs between two species.
        
        out: tad1 tad2 geneNum1 geneNum2 synNum1 synNum2 \
            genePer1 genePer2 synPer1 synPer2 geneList1 geneList2
        
        >>> tc = TADConserved
        >>> tc.getConserved(tad1, tad2, syngene1, syngene2, gene1, gene2, anchor)
        ...
        """
        logging.debug('Start ...')
        check_file_exists(anchors)
        tc = TADConserved()
        tadSynGeneNum1 = tc.getGene(tad1, syngene1, fraction, isnum=True)
        tadSynGeneNum2 = tc.getGene(tad2, syngene2, fraction, isnum=True)
        tadGeneNum1 = tc.getGene(tad1, gene1, fraction, isnum=True)
        tadGeneNum2 = tc.getGene(tad2, gene2, fraction, isnum=True)
        geneTAD1 = tc.genePairTAD(syngene1, tad1, fraction)
        geneTAD2 = tc.genePairTAD(syngene2, tad2, fraction)

        db = OrderedDict()
        with open(anchors, 'r') as fp:
            for line in fp:
                if line[0] == "#":
                    continue
                gene1, gene2, length = line.strip().split()

                try:
                    anchor1 = geneTAD1[gene1]
                    anchor2 = geneTAD2[gene2]
                except KeyError:
                    continue
                if anchor1 not in db:
                    db[anchor1] = OrderedDict()
                if anchor2 not in db[anchor1]:
                    db[anchor1][anchor2] = []
                db[anchor1][anchor2].append((gene1, gene2))
        header = ('#tad1', 'tad2', 'total_gene_num1', 'total_gene_num2',
                  'syn_gene_num1', 'syn_gene_num2', 'gene_per1', 'gene_per2',
                  'syngene_per1', 'syngene_per2', 'gene_list1', 'gene_list2')
        print("\t".join(header), file=sys.stdout)
        for anchor1 in db:
            for anchor2 in db[anchor1]:
                tmp = np.array(db[anchor1][anchor2])
                geneNum1 = tadGeneNum1[anchor1]
                geneNum2 = tadGeneNum2[anchor2]
                synGeneNum1 = tadSynGeneNum1[anchor1]
                synGeneNum2 = tadSynGeneNum2[anchor2]
                genePer1 = len(tmp[:, 0]) * 1.0 / geneNum1
                genePer2 = len(tmp[:, 1]) * 1.0 / geneNum2
                synGenePer1 = len(tmp[:, 0]) * 1.0 / synGeneNum1
                synGenePer2 = len(tmp[:, 1]) * 1.0 / synGeneNum2
                if genePer1 >= threshold and genePer2 >= threshold and \
                        geneNum1 >= gene_num and geneNum2 >= gene_num and \
                            synGeneNum1 >= synthre and synGeneNum2 >= synthre:

                    print("\t".join(
                        map(str, (anchor1, anchor2, geneNum1, geneNum2,
                                  synGeneNum1, synGeneNum2, genePer1, genePer2,
                                  synGenePer1, synGenePer2, ",".join(
                                      tmp[:, 0]), ",".join(tmp[:, 1])))),
                          file=sys.stdout)
        logging.debug('Done')
    
    @staticmethod
    def randomTAD(parameter_list):
        pass


class CalObsExp(object):
    """
    From https://www.nature.com/articles/s41477-019-0479-8#Sec23
    To Calculate the observed/expected 
    Parameters
    ----------
    object : [type]
        [description]
    
    Returns
    -------
    [type]
        [description]
    """
    def __init__(self, matrix_source, outpre):
        lib = cooler.Cooler(matrix_source)
        for c in lib.chromnames:
            raw = lib.matrix(balance=False, sparse=False).fetch(c)
            raw[np.isnan(raw)] = 0
            expected = self.expected_matrix(raw)
            obs_exp = raw / expected
            obs_exp[expected==0] = 0
            outfil = outpre + '.{0}.npy'.format(c)
            np.save(outfil, obs_exp)

    def expected_matrix(self, raw):

        tmp = raw.sum(axis=0)!=0 # valid rows or columns
        n = raw.shape[0]
        expected = np.zeros_like(raw)
        idx = np.arange(n)
        for i in idx:
            if i > 0:
                valid = tmp[:-i] * tmp[i:]
            else:
                valid = tmp
                current = raw.diagonal(i )[valid]
            if current.size > 0:
                v = current.mean()
                if i > 0:
                    expected[idx[:-i], idx[i:]] = v
                    expected[idx[i:], idx[:-i]] = v
                else:
                    expected[idx, idx] = v
        return expected


def test(args):
    tad1, tad2, syngene1, syngene2, gene1, gene2, anchor = args
    TADConserved.getConserved(tad1, tad2, syngene1, syngene2, gene1, gene2,
                              anchor)


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
    p.add_option('-a',
                 '--up',
                 type=int,
                 default=0,
                 help='the upstrean distance of boundary '
                 '[default: %default]')
    p.add_option('-b',
                 '--down',
                 type=int,
                 default=1,
                 help='the downstream distance of boundary '
                 '[default: %default]')

    opts, args = p.parse_args(args)
    if len(args) < 2:
        sys.exit(p.print_help())

    tadFile, chromSize = args
    check_file_exists(chromSize)
    up, down = opts.up, opts.down
    if not op.exists(tadFile) or not \
            op.exists(chromSize):
        logging.error('The input file is not exists')
    chrom_dict = dict(i.strip().split() \
            for i in open(chromSize) if i.strip())

    tf = TADFile(tadFile)
    tf.getBoundaryDict()
    boundaryBed = tf.getBoundaryBed(tf.boundaryDict, chrom_dict, up, down)

    for item in sorted(boundaryBed):
        print("\t".join(map(str, item[:3])))
    logging.debug('Successful output boundary bed')


def plotSizeDist(args):
    """
    %prog plotSizeDist <tad1.bed> [tad2.bed ...] [Options]

    Given some tad bed file to plot their sizes distributions.
    """
    scale_units = {1: 'bp', 1000: 'kb', '1e6': 'Mb'}
    p = OptionParser(plotSizeDist.__doc__)
    p.add_option('-o',
                 '--out',
                 default='tad_sizes_dist.pdf',
                 help='out of plot [default: %default]')
    p.add_option('--all',
                 default=False,
                 action='store_true',
                 help='plot all levels of tads [default: %default]')
    p.add_option('-s',
                 '--scale',
                 default=1000,
                 type=int,
                 help='the scale of xticks [default: %default]')
    p.add_option('--xmin',
                 default=0,
                 type=int,
                 help='min value of xticks [default: %default]')
    p.add_option('--xmax',
                 default=800,
                 type=int,
                 help='max value of xticks [default: %default]')
    p.add_option('--step',
                 default=100,
                 type=int,
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
        label = tad.rsplit("_")[0]
        tf = TADFile(tad)
        data = tf.Sizes if opts.all else tf.bottomSizes
        ax = tf.plotSizeDistMulti(ax, data, label=label, scale=scale)
    ax.set_xlim(xmin, xmax)
    ax.tick_params(labelsize=11)
    ax.set_xticks(range(xmin, xmax + 1, step))
    ax.set_xticklabels(list(range(xmin, xmax+1, step)), rotation=45, ha='right')
    ax.set_xlabel('TADs Size ({})'.format(scale_units[scale]), fontsize=13)
    ax.set_ylabel('Frequency', fontsize=13)
    ax.set_title('TADs Size Distribution', fontsize=14)
    plt.savefig(out, dpi=300, bbox_inches='tight')
    logging.debug('Success file is save as {}'.format(out))


def stat(args):
    """
    tads informations stat
        total_num total_size genome_size percentage
    """
    p = OptionParser(stat.__doc__)
    p.add_option('-g', '--genome', type=int, help='the genome size of species')
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())
    if not opts.genome:
        logging.error('Must input genome size ' 'with `-g` option')
        sys.exit()
    tad, = args
    tf = TADFile(tad)
    total_num, total_size = 0, 0
    for size in tf.bottomSizes:
        total_num += 1
        total_size += size

    print('#Total number\tTotal size\tGenome size\tPercentage')
    print("{}\t{}\t{}\t{:.2%}".format(total_num, total_size, opts.genome,
                                      total_size * 1.0 / opts.genome))


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
    p.add_option('-F',
                 dest='fraction',
                 default='0.7',
                 help='the fraction of gene overlap of tads'
                    ' [default: %default]')
    p.add_option('--isnum',
                 default=False,
                 action='store_true',
                 help='if output the gene number [default: %default]')
    p.add_option('--plot',
                 default=False,
                 action='store_true',
                 help='if plot the gene number '
                 'distribution [default: %default]')
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    tads, genes = args
    fraction = opts.fraction
    db = TADConserved().getGene(tads, genes, fraction, opts.isnum, opts.plot)
    if opts.isnum:
        return
    for ID in db:
        gene_list = sorted(db[ID])
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
    p.add_option('-f',
                 dest='fraction',
                 default='0.7',
                 help='the fraction of gene overlap of tads')

    opts, args = p.parse_args(args)

    if len(args) != 2:
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


def getSyntenicTADs(args):
    """
    %prog tad1.bed tad2.bed synteny1.bed synteny2.bed 
        gene1.bed gene2.bed 1.2.anchors [Options]
    %prog species1 species2
        Cautions: if input species, these file should exists
                `species1.tad.bed`,  tad bed3
                `species2.tad.bed`,
                `species1.synteny.bed`, synteny gene bed4
                `species2.synteny.bed`,
                `species1.bed`, all gene bed4
                `species2.bed`,
                `species1.species2.anchors`. synteny anchors file

    To get syntenic TADs table, default is not filter ouput all result.
    """
    p = OptionParser(getSyntenicTADs.__doc__)
    p.add_option('--fraction',
                 default='0.7',
                 help='fraction of gene overlap with '
                 'TADs [defalut: %default]')
    p.add_option('--threshold',
                 default=0,
                 type=float,
                 help='the threshold of non-change syn-gene / total gene num'
                 '[default: %default]')
    p.add_option('--gene_num',
                 default=0,
                 type=int,
                 help='the least gene number of TAD [default: %default]')
    p.add_option('--synthre',
                 default=0,
                 type=float,
                 help='the threshold of non-change syn-gene / syn-gene num')

    opts, args = p.parse_args(args)
    if len(args) == 2:
        logging.debug('less args mode, Input two species prefix')
        species1, species2 = args
        tad1 = species1 + ".tad.bed"
        tad2 = species2 + ".tad.bed"
        syngene1 = species1 + '.synteny.bed'
        syngene2 = species2 + '.synteny.bed'
        gene1 = species1 + '.bed'
        gene2 = species2 + '.bed'
        anchor = species1 + '.' + species2 + '.anchors'

    elif len(args) == 7:
        tad1, tad2, syngene1, syngene2, gene1, gene2, anchor = args
    else:
        sys.exit(p.print_help())

    TADConserved.getConserved(tad1, tad2, syngene1, syngene2, gene1, gene2,
                              anchor, opts.fraction, opts.threshold,
                              opts.gene_num, opts.synthre)


def plotBoundary(args):
    """
    %prog boundary.bed data.bw samplelabel [options]
        To plot omics data density in tads boundary.
    """
    p = OptionParser(plotBoundary.__doc__)
    p.add_option('-b',
                 dest='up',
                 default=50000,
                 type=int,
                 help='upstream distance of boundary [default: %default]')
    p.add_option('-a',
                 dest='down',
                 default=50000,
                 type=int,
                 help='downstream distance of boundary [default: %default]')
    p.add_option('--binSize',
                 default=1000,
                 type=int,
                 help='calculate binSize [default: %default]')
    p.add_option('-p',
                 '--process',
                 default=8,
                 type=int,
                 help='process of program [default:%default]')
    p.add_option('-o',
                 '--output',
                 default=None,
                 help='the plot output prefix [default: tadprefix_label]')
    opts, args = p.parse_args(args)
    if len(args) != 3:
        sys.exit(p.print_help())

    boundary, data, label = args
    check_file_exists(boundary)
    check_file_exists(data)
    up = opts.up
    down = opts.down
    binSize = opts.binSize
    process = opts.process

    prefix = op.basename(boundary).replace('.bed', '') if not opts.output \
            else opts.output
    compute_cmd = """
    computeMatrix reference-point -S {data} -R {boundary} \\
        --referencePoint center -b {up} -a {down} --binSize {binSize} \\
            --samplesLabel {label} -p {process} --missingDataAsZero \\
                --skipZeros -o {prefix}_{label}.matrix.gz\\
                    --outFileSortedRegions {prefix}_{label}.bed
    """.format(data=data,
               boundary=boundary,
               binSize=binSize,
               up=up,
               down=down,
               label=label,
               prefix=prefix,
               process=process)

    plot_cmd = """
     plotProfile -m {prefix}_{label}.matrix.gz --refPointLabel Boundary \\
         -out {prefix}_{label}.pdf --plotHeight 10 --plotWidth 12 
    """.format(prefix=prefix, label=label)

    with open('run_{}_{}.sh'.format(prefix, label), 'w') as out:
        out.write(compute_cmd + "\n")
        out.write(plot_cmd)
    logging.debug('Starting plot {} density in boundary'.format(label))
    os.system('sh run_{}_{}.sh'.format(prefix, label))
    logging.debug('Done, picture is `{prefix}_{label}.pdf`'.format(
        prefix=prefix, label=label))




def quickPlotTAD(args):
    """
    %(prog)s <sample_4000.iced.cool> <sample.domain> [Options]

        Quick plot picture to view all TADs results.

    """

    p=argparse.ArgumentParser(prog=quickPlotTAD.__name__,
                        description=quickPlotTAD.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cool', help='cool file of hic matrix')
    pReq.add_argument('domain', help='domain file of TAD, '
                        'three columns(chrom start end).')
    pReq.add_argument('chrom_size', help='chromosome sizes file')
    pOpt.add_argument('--bg', default=None, 
            help='Direction index bedGraph file [default: %(default)s]')

    pOpt.add_argument('--min_value', default=3, type=str, 
            help='min value of hic matrix [default: %(default)s]')
    pOpt.add_argument('-d', '--depth', default=None, type=str,
            help='hicmatrix depth [default: window*1.2].')
    pOpt.add_argument('-w', '--window', type=float, default=5e6,
            help='window of chromosome sizes [default: %(default)s]')
    pOpt.add_argument('-o', '--outdir', default='quickPlotTAD_result',
            help='outdir of result [default: %(default)s]')
    pOpt.add_argument('--pdf', default=False, action='store_true',
            help='if output pdf format [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    if not args.depth:
        depth = 1.2 * args.window
    else:
        depth = args.depth
    

    import configparser
    from TDGP.apps.utilities import makeChromWindows
    cf = configparser.ConfigParser()

    cf.add_section('hic matrix')
    cf.set('hic matrix', 'file', args.cool)
    cf.set('hic matrix', 'title', 'Hi-C')
    cf.set('hic matrix', 'depth', str(int(depth)))
    cf.set('hic matrix', 'min_value', str(args.min_value))
    cf.set('hic matrix', 'transform', 'log1p')
    cf.set('hic matrix', 'file_type', 'hic_matrix')

    cf.add_section('tad')
    cf.set('tad', 'file', args.domain)
    cf.set('tad', 'display', 'triangles')
    cf.set('tad', 'border_color', 'black')
    cf.set('tad', 'color', 'none')
    cf.set('tad', 'overlay_previous', 'share-y')

    
    if args.bg:
        cf.add_section("spacer")
        cf.add_section("DI_bg")
        cf.set('DI_bg', 'file', args.bg)
        cf.set('DI_bg', 'height', '4')
        cf.set('DI_bg', 'title', 'DI')
        cf.set('DI_bg', 'negative_color', "#0B1D51")
        cf.set('DI_bg', 'color', '#787596')

    
    cf.add_section('x-axis')
    cf.set('x-axis', 'where', 'bottom')
    if not op.exists(args.outdir):
        os.makedirs(args.outdir)
    with open('{}/quickPlotTAD.tad.ini'.format(args.outdir), 'w+') as f:
        cf.write(f)
    
    chrom_windows_db = makeChromWindows(args.chrom_size, args.window)
    plot_cmd_formatter = "pyGenomeTracks --tracks {1}/quickPlotTAD.tad.ini -o {1}/{0}.{2} --region {0}"
    
    ext = 'pdf' if args.pdf else 'png'
    
    for chrom in chrom_windows_db:
        for (start, end) in chrom_windows_db[chrom]:
            region = '{}:{}-{}'.format(chrom, start, end)
            print(plot_cmd_formatter.format(region, args.outdir, ext),
                    file=sys.stdout)
    

if __name__ == "__main__":
    main()
