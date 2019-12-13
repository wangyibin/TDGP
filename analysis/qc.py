#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Hi-C analysis quality control, such as IDE ...
"""

from __future__ import print_function

import logging
import numpy as np
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import os.path as op
import sys

from collections import OrderedDict, defaultdict
from itertools import chain
from optparse import OptionParser
from TDGP.apps.base import ActionDispatcher
from TDGP.apps.base import check_file_exists, debug
from TDGP.apps.base import listify
from TDGP.apps.base import BaseFile


debug()

def main():

    actions = (
            ("plotDistDensity", "Plot the IDE"),
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
    def __init__(self, infile):
        super(ValidPairs, self).__init__(infile)
        self.infile = infile
        
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
    
    def getCisDistance(self, chrom=None):
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
        cis_dist_db = defaultdict(list)
        
        for vpl in self.getCisLine():
            if chrom:
                chrom = listify(chrom)
                if vpl.chr1 in chrom:
                    cis_dist_db[vpl.chr1].append(vpl.getCisDistance())
            else:
                cis_dist_db[vpl.chr1].append(vpl.getCisDistance())
        self.cis_dist_db = cis_dist_db
        return self.cis_dist_db
    
    @classmethod
    def plotDistDensity(self, distance_db, out, perchrom=True, scale=100000):
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

        plt.figure(figsize=(5, 5))
        if perchrom:
            for chrom in distance_db: 
                data = np.array(distance_db[chrom]) // scale * scale
                data = data[(data >= 1e5) & (data <= 2e7)]
                unique, counts = np.unique(data, return_counts=True)
                db = OrderedDict(zip(unique, counts))
                plt.plot(db.keys(), db.values(), label=chrom)
                plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        else:
            data = list(chain(*list(distance_db.values())))
            data = np.array(data) // scale * scale
            data = data[(data >= 1e5) & (data <= 2e7)]
            unique, counts = np.unique(data, return_counts=True)
            db = OrderedDict(zip(unique, counts))
            plt.plot(db.keys(), db.values(), label='all chromosome')
            plt.legend(loc='best')
        plt.ylabel('$\log_{10}$ (contact frequency)')
        plt.xlabel( '(distance)')
        plt.yscale('log')
        plt.xscale('log')
        plt.savefig(out, dpi=300, bbox_inches='tight')

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
    
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())
    
    pairsFile, out = args
    if opts.chrom:
        chrom = opts.chrom.split(',')
    else:
        chrom = opts.chrom
    vp = ValidPairs(pairsFile)
    distance_db = vp.getCisDistance(chrom=chrom)
    print(opts.perchr)
    vp.plotDistDensity(distance_db, out, perchrom=opts.perchr, scale=opts.scale)

if __name__ == "__main__":
    main()