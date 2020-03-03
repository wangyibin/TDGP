#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Hi-C analysis quality control, such as IDE ...
"""

from __future__ import print_function

import argparse
import logging
import numpy as np


import cooler
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
import os
import os.path as op
import sys

from collections import OrderedDict, defaultdict
from itertools import chain
from joblib import Parallel, delayed, Memory
from optparse import OptionParser
from TDGP.apps.base import ActionDispatcher
from TDGP.apps.base import check_file_exists, debug
from TDGP.apps.base import listify
from TDGP.apps.base import BaseFile
from TDGP.formats.hicmatrix import cool2matrix


debug()

def main():

    actions = (
            ('plotCisTrans', 'plot the barplot of cis and trans interactions'),
            ("plotDistDensity", "Plot the IDE"),
            ("plotIDEMulti", 'plot multi sample IDE'),
            ('statFrag', 'stat the reality fragments')
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


class ValidPairsLine(object):
    """
    The HiCPro ValidPairs line object.
    Returns:
    ---------
    Examples:
    ----------
    >>> vp = ValidParisLine(line)
    >>> vp.chr1
    Chr1
    >>> vp.pos2 
    2341
    >>> print(vp)
    Chr1 2342 + Chr2 324 - ....
    """
    def __init__(self, line):
        self.line = line
        if not self.line:
            logging.error('Nothing item in this line.')
            exit(1)
        self.line_list = self.line.strip().split()
        (self.read, self.chr1, self.pos1, self.strand1, 
         self.chr2, self.pos2, self.strand2, self.size, 
         self.site1, self.site2, _, _) = self.line_list
        self.pos1, self.pos2 = int(self.pos1), int(self.pos2)
        self.size = int(self.size)
        #self.distance = self.pos2 - self.pos1
    def isCis(self):
        """
        If pairs is Cis return True.
        
        Returns:
        ---------
        out: bool
        True or False of pairs whether is Cis.
        
        Examples:
        ----------
        >>> vpl.isCrans()
        True
        """
        if self.chr1 == self.chr2:
            return True
        else:
            return False
    
    def isTrans(self):
        """
        If pairs is Trans return True.
        
        Returns:
        ---------
        out: bool
        True or False of pairs whether is trans.
        
        Examples:
        ----------
        >>> vpl.isTrans()
        True
        """
        if self.chr1 != self.chr2:
            return True
        else:
            return False
    
    def getCisDistance(self):
        """
        Calculate the distance bewteen pairs.
        Returns:
        ---------
        out: int
        
        Examples:
        ---------
        >>>vpl.getCisDistance()
        23424
        """
        if self.isCis():
            if self.pos1 < self.pos2:
                distance = self.pos2 - self.pos1
            else:
                distance = self.pos1 - self.pos2
        else:
            distance = 0
        return distance
    
    def __str__(self):
        return self.line.strip()
    __repr__ = __str__
    
    
class ValidPairs(BaseFile):
    """
    Class to handle validPairs from HiCPro result.
    Also contain some util function to deal validPairs.
    
    Examples:
    ----------
    >>> vp = ValidPairs('sample.validPairs')
    >>> vp.getLine():
    
    """
    def __init__(self, infile, mem_cache='.'):
        super(ValidPairs, self).__init__(infile)
        self.infile = infile
        self.mem_cache = mem_cache
        memory = Memory(self.mem_cache)
        self.getCisDistance = memory.cache(self._getCisDistance)
        
    def getLine(self):
        """
        Get Line of validPairs.
        
        Returns:
        ----------
        out:  object of `ValidParisLine`
        
        Examples:
        ----------
        """
        with open(self.infile) as fp:
            for line in fp:
                yield ValidPairsLine(line)
    
    def __iter__(self):
        for vpl in self.getLine():
            yield vpl
            
    def getCisLine(self):
        """
        Get all cis pairs line.
        
        Returns:
        ---------
        out: ValidPairsLine
        
        Examples:
        ---------
        >>> next(vp.getCisLine)
        read1 Chr1 2343 ...
        """
        for vpl in self.getLine():
            if vpl.isCis():
                yield vpl
    
    def getTransLine(self):
        """
        Get all trans pairs line.
        
        Returns:
        ---------
        out: ValidPairsLine
        
        Examples:
        ---------
        >>> next(vp.getTransLine)
        read1 Chr1 2343 ...
        """
        for vpl in self.getLine():
            if vpl.isTrans():
                yield vpl
    
    def _getCisDistance(self, chrom=None):
        """
        Get all chromosome cis distance.
        
        Returns:
        ----------
        out: dict
        
        Examples:
        ----------
        >>> vp.getCisDistance()
        {'Chr1': [32, 4434, 23223, ...], 'Chr2': [2342, ...]}
        """
        cis_dist_db = OrderedDict()
        
        for vpl in self.getCisLine():
            if chrom:
                chrom = listify(chrom)
                if vpl.chr1 in chrom:
                    if vpl.chr1 not in cis_dist_db:
                        cis_dist_db[vpl.chr1] = []
                    cis_dist_db[vpl.chr1].append(vpl.getCisDistance())
            else:
                if vpl.chr1 not in cis_dist_db:
                    cis_dist_db[vpl.chr1] = []
                cis_dist_db[vpl.chr1].append(vpl.getCisDistance())
        self.cis_dist_db = cis_dist_db
        return self.cis_dist_db
    
    def getRealFrags(self):
        """
        Obtain the reality fragments, which is the mapped restrict
            enzyme fragments.
        
        Returns:
        --------
        out: `set` 

        Examples:
        --------
        >>> vp = ValidPairs('allvalidpairs')
        >>> vp.getRealFrags()
        set('HiC_1', 'HiC_2' ...)
        """
        self.realFrags = set()
        for vpl in self.getLine():
            self.realFrags.add(vpl.site1)
            self.realFrags.add(vpl.site2)
        return self.realFrags 


    @classmethod
    def plotDistDensity(self, distance_db, out, perchrom=True, scale=100000,
            xmin=1e5, xmax=1e8,):# color=[]):
        """
        Plot the density of contact distance per chromosome or whole chromosome
        
        Params:
        --------
        distance_db: `dict` or per chromosome distance
        perchrom: `bool` default=True
        scale: `int` default=100000
        
        Returns:
        --------
        out: figure or distance density.
        
        Examples:
        ---------
        >>> vp = ValidPairs('all.validpairs')
        >>> distance_db = vp.getCisDistance()
        >>> out = 'ide.pdf'
        >>> plotDistDensity(distance_db, out)
        """
        """
        if color:
            color_pallete = listify(color)
            single_color = listify(color)[0]
        else:
            single_color = '#209093'
            if len(distance_db) <= 8:
                color_pallete = sns.color_palette('Set2')
            else:
                color_pallete = sns.color_palette('hls', len(distance_db))
        """
        from scipy.stats import linregress
        plt.figure(figsize=(5, 5))
        if perchrom:
            for i, chrom in enumerate(distance_db): 
                #c = color_pallete[i % len(color_pallete)]
                data = np.array(distance_db[chrom]) // scale * scale
                data = data[(data >= xmin) & (data <= xmax)]
                unique, counts = np.unique(data, return_counts=True)
                db = OrderedDict(zip(unique, counts))
                slope = linregress(np.log10(unique), np.log10(counts)).slope
                label = "{} ({:.2f})".format(chrom, slope)
                plt.plot(list(db.keys()), list(db.values()), label=label,)# color=c)
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        else:
            data = list(chain(*list(distance_db.values())))
            data = np.array(data) // scale * scale
            data = data[(data >= xmin) & (data <= xmax)]
            unique, counts = np.unique(data, return_counts=True)
            db = OrderedDict(zip(unique, counts))
            slope = linregress(np.log10(unique), np.log10(counts)).slope
            label = 'all'
            label = "{} ({:.2f})".format(label, slope)
            plt.plot(list(db.keys()), list(db.values()), 
                label=label)#, color=single_color)
            plt.legend(loc='best', fontsize=14)
        #plt.xlim(xmin, xmax)
        plt.ylabel('Contact probability', fontsize=14)
        plt.xlabel('Distance (bp)', fontsize=14)
        plt.yscale('log')
        plt.xscale('log')
        plt.savefig(out, dpi=300, bbox_inches='tight')
        plt.savefig(out.rsplit(".", 1)[0] + '.png', 
                    dpi=300, bbox_inches='tight')
        logging.debug('Successful, picture is in `{}`'.format(out))


## outside command ##
def plotDistDensity(args):
    """
    %prog all.validpairs [Options]
        Plot the IDE of all genome or per chromosome
    """
    p = OptionParser(plotDistDensity.__doc__)
    p.add_option('--chrom', default=None, 
            help='plot chrom list')
    p.add_option('--perchr', default=True, action='store_false',
            help='whether to plot per chromosome [default: %default]')
    p.add_option('-s', '--scale', default=100000, type=int,
            help='the scale of data [default: %default]')
    p.add_option('--xmin', default=1e5, type=float,
            help='min value of xtick [default: %default]')
    p.add_option('--xmax', default=1e8, type=float, 
            help='max value of xtick [default: %default]')
    p.add_option('--color', default='', 
            help='color palette')
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())
    
    pairsFile, out = args
    if opts.chrom:
        if op.exists(opts.chrom):
            chrom = [i.strip().split()[0] for i in open(opts.chrom) if i.strip()]
        else:
            chrom = opts.chrom.split(',')
    else:
        chrom = opts.chrom
    """
    if opts.color:
        if op.exists(opts.color):
            color = [i.strip() for i in open(opts.color) if i.strip()]
        else:
            color = opts.color.split(',')
    else:
        color = opts.color
    """
    vp = ValidPairs(pairsFile)
    distance_db = vp.getCisDistance(chrom=chrom)
    vp.plotDistDensity(distance_db, out, perchrom=opts.perchr, scale=opts.scale,
            xmin=opts.xmin, xmax=opts.xmax) #, color=color)


def plotIDEMulti(args):
    """
    %(prog) 1.ValidPairs 2.ValidPairs ... [Options]
        To multi sample IDE in a picture.
    """
    p = p=argparse.ArgumentParser(prog=plotIDEMulti.__name__,
                        description=plotIDEMulti.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('validpairs', nargs="+", 
            help='validpairs file')
    pReq.add_argument('--labels', nargs='+', required=True,
            help='lable for legend')
    pReq.add_argument('-o', '--out', required=True,
            help='output file')
    pOpt.add_argument('--chrom', default=None, help='plot chrom list')
    pOpt.add_argument('--scale', default=100000, type=int, metavar='int',
            help='the scale of data [default: %(default)]')
    p.add_argument('--xmin', default=1e5, type=float, metavar='float',
            help='min value of xtick [default: %(default)]')
    p.add_argument('--xmax', default=2e7, type=float, metavar='float',
            help='max value of xtick [default: %(default)]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    from scipy.stats import linregress
    from matplotlib.lines import Line2D
    scale = args.scale
    xmin = args.xmin
    xmax = args.xmax
    out = args.out
    fig, ax = plt.subplots(figsize=(5, 5))

    if args.chrom:
        if op.exists(args.chrom):
            chrom = [i.strip().split()[0] 
                        for i in open(args.chrom) if i.strip()]
        else:
            chrom = args.chrom.split(',')
    else:
        chrom = args.chrom
    for i in args.validpairs:
        check_file_exists(i)
    assert len(args.validpairs) == len(args.labels), \
        'input validpair file must equal to labels'
    i = 0
    for validpair, label in zip(args.validpairs, args.labels):
        vp = ValidPairs(validpair)
        distance_db = vp.getCisDistance(chrom=chrom)
        data = list(chain(*list(distance_db.values())))
        data = np.array(data) // scale * scale
        
        data = data[(data >= xmin) & (data <= xmax)]
       
        unique, counts = np.unique(data, return_counts=True)
        db = OrderedDict(zip(unique, counts))
        slope = linregress(np.log10(unique), np.log10(counts)).slope
        label = "{} ({:.2f})".format(label, slope)
        #sns.regplot(list(db.keys()), list(db.values()), label=label, 
         #   marker=Line2D.filled_markers[i], ci=0, truncate=True,
        #    )
        #i += 1
        plt.plot(list(db.keys()), list(db.values()), label=label)
    plt.legend(loc='best', fontsize=13)
    #plt.xlim(xmin, xmax)
    plt.ylabel('Contact probability', dict(size=14))
    plt.xlabel('Distance (bp)', dict(size=14))
    plt.yscale('log')
    plt.xscale('log')
    #sns.despine(trim=True)
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.savefig(out.rsplit(".", 1)[0] + '.png', 
                    dpi=300, bbox_inches='tight')
    logging.debug('Successful, picture is in `{}`'.format(out))


def plotCisTrans(args):
    """
    %(prog)s <coolfile> [coolfile ...] [Options]

        calculate the cis and trans interaction, and plot the barplot.
    
    """
    p = p=argparse.ArgumentParser(prog=plotCisTrans.__name__,
                        description=plotCisTrans.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cool', nargs='+', 
            help='cool file of hicmatrix')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    cis_list = []
    trans_list = []
    for coolfile in args.cool:
        cool = cooler.Cooler(coolfile)
        hm = cool2matrix(cool)
        # calculate cis counts
        counts = np.zeros(cool.info['nbins'])
        for chrom in cool.chromnames:
            idx = cool.bins().fetch(chrom).index
            counts[idx] = np.triu(hm[idx][:, idx]).sum(axis=1)
        cis_list.append(int(counts.sum()))

        ## calculate trans counts
        counts = np.zeros(cool.info['nbins'])
        for chrom in cool.chromnames:
            idx = cool.bins().fetch(chrom).index
            start = idx[0]
            end = idx[-1]
            hm[start: end + 1, start: end + 1] = 0
            counts[idx] = hm[idx].sum(axis=1)
        counts = counts / 2
        trans_list.append(int(counts.sum()))
        #print(counts.sum())

    
    print('\t'.join(map(str, cis_list)))
    print('\t'.join(map(str, trans_list)))

def statFrag(args):
    """
    %(prog)s allValidParis [Options]
        stat the Ratio of theoretically digested genomic 
        fragments covered by valid paired Hi-C reads.
    """
    p = p=argparse.ArgumentParser(prog=statFrag.__name__,
                        description=statFrag.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('validpairs',  help='Validparis file')
    pReq.add_argument('enzyme', help='restriction enzyme site bed file')
    pOpt.add_argument('-o', '--out', type=argparse.FileType('w'), default=sys.stdout,
            help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    theo_num = 0
    reality_num = 0
    check_file_exists(args.enzyme)

    with open(args.enzyme, 'r') as fp:
        for line in fp:
            if line.strip():
                theo_num += 1
    vp = ValidPairs(args.validpairs)
    realFrags = vp.getRealFrags()
    reality_num = len(realFrags)
    
    print("Theoretical Fragments\t{}".format(theo_num), file=args.out)
    print("Reality Fragments\t{}".format(reality_num), file=args.out)
    print("Reality Fragments Ratio (%)\t{:.2%}".format(reality_num*1.0/theo_num), 
            file=args.out)


if __name__ == "__main__":
    main()