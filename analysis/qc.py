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
import glob
import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy
import seaborn as sns
import os
import os.path as op
import pandas as pd
import sys

from argparse import Namespace
from collections import OrderedDict, defaultdict
from itertools import chain
from joblib import Parallel, delayed, Memory
from optparse import OptionParser
from TDGP.apps.base import ActionDispatcher
from TDGP.apps.base import check_file_exists, debug
from TDGP.apps.base import listify
from TDGP.apps.base import BaseFile
from TDGP.formats.bedGraph import import_bedgraph
from TDGP.formats.hicmatrix import cool2matrix, getCisCounts, getTransCounts


debug()

def main():

    actions = (
            ('validStat', 'stat hicpro valid data'),
            ('plotCisTrans', 'plot the barplot of cis and trans interactions'),
            ("plotDistDensity", "Plot the IDE"),
            ("plotIDEMulti", 'plot multi sample IDE'),
            ('plotFragmentSizeHist', 'plot histogram of fragment size'),
            ('getRealityBed', 'get reality fragments bed file'),
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
         self.site1, self.site2) = self.line_list[:10]

        self.pos1, self.pos2 = int(self.pos1), int(self.pos2)
        self.size = int(self.size)
        #self.distance = self.pos2 - self.pos1
    def isCis(self):
        """
        If pairs is Cis return True.
        
        Returns:
        ---------
        out: `bool`,
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
        out: `bool`,
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
            xmin=1e5, xmax=1e8, plotSlope=False, slopeRange='500000-7000000'):# color=[]):
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
        slope_left, slope_right = list(map(int, slopeRange.split("-")))
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
            slope_data = data[(data >= slope_left) & (data <= slope_right)]
            slope_unique, slope_counts = np.unique(slope_data, return_counts=True)
            slope_unique_log, slope_counts_log = np.log10(slope_unique), np.log10(slope_counts)
            slope, intercept, rvalue, pvalue, stderr = linregress(slope_unique_log, 
                                                                slope_counts_log)
            label = 'Intrachromosomal'
            if plotSlope:
                label = "{}".format(label)
            else:
                label = "{} ({:.2f})".format(label, slope)
            
            plt.plot(list(db.keys()), list(db.values()), 
                label=label)#, color=single_color)
            
            if plotSlope:
                slope_x_values = slope_unique
                slope_y_values = 10 ** (slope_unique_log * slope + intercept)
                slope_label = "{:.2f}".format(slope)
                plt.plot(slope_x_values, slope_y_values, "--", label=slope_label, lw=2)
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

def getCisTransData(cool, sample='sample'):
    """
    get cis/trans ratio per chromosome
    
    Params:
    --------
    cool: `str`
            path of cool file

    sample: `str`, optional
            label of cool file, also show in barplot
    
    Returns:
    --------
    data: `DataFrame`
    
    Examples:
    ---------
    >>> cool = '/path/to/sample.cool'
    >>> sample = 'sample'
    >>> data = getCisTransData(cool, sample)
    """
    prefix = op.splitext(cool)[0] 
    cis_bg_file = prefix + ".cis.bg"
    trans_bg_file = prefix + ".trans.bg"
    getCisCounts([cool, '-o', cis_bg_file])
    cis_bg = import_bedgraph(cis_bg_file)
    getTransCounts([cool, '-o', trans_bg_file])
    trans_bg = import_bedgraph(trans_bg_file)
    data = cis_bg.groupby('chrom').score.sum() / trans_bg.groupby('chrom').score.sum()
    data = data.to_frame()
    data.reset_index(inplace=True)
    data['sample'] = sample
    return data


def concatCisTransData(cools, samples):
    """
    concat cis/trans data from cools file
    
    Params:
    --------
    cools: `list`
            cool file path list

    sampels: `list`
            label of barplot
    
    Returns:
    --------
    data: `DataFrame`
            data of cis/trans ratio per chromosome
    
    Examples:
    ---------
    >>> cools = ['sample.cool', 'sample2.cool']
    >>> samples = ['sample', 'sample2']
    >>> data = concatCisTransData(cools, samples)
    
    """
    cool_dict = dict(zip(samples, cools))
    data_list = []
    for sample in cool_dict:
        data_list.append(getCisTransData(cool_dict[sample], sample))
    data = pd.concat(data_list, axis=0)
    return data

def plotCisTransPerChrom(data, rotation=True, 
                            legend=True, colors=None):
    """
    barplot of cis/trans ratio.
    
    Params:
    --------
    data: `DataFrame`
            from getCisTransData or concatCisTransData

    rotation: `bool`, optional
            wether rotation xticks [default: True]

    legend: `bool`, optional
            wether show legend [default: True]

    colors: `list`, optional
            color list for par plot per sample [default: None]
    
    Returns:
    ---------
    ax: `axes`
    
    Examples:
    ---------
    >>> data = getCisTransData(cool)
    >>> plotCisTransPerChrom(data)
    
    """
    rotation = 45 if rotation else 0
    if colors:
        ax = sns.barplot(x='chrom', y='score', data=data, 
                         hue='sample', pallete=colors)
    else:
        ax = sns.barplot(x='chrom', y='score', data=data, 
                         hue='sample')
    
    ax.tick_params(axis='x', rotation=0)
    ax.set_ylim(0, data.score.max() * 1.4)
    ax.set_xlabel('')
    ax.set_ylabel('$\mathit{{cis}}$ / $\mathit{{trans}}$', 
                  fontsize=14)

    ax.tick_params(which='major', axis='y', labelsize=13, width=2)
    ax.tick_params(which='major', axis='x', labelsize=13, width=0)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(0)
    sns.despine(trim=True)
    if legend:
        plt.legend()
    else:
        plt.legend([],[], frameon=False)

    return ax

## outside command ##

def validStat(args):
    """
    %(prog)s <statpath> <outfile>

        Stat the hicpro valid pais dat to a table
    """
    p = argparse.ArgumentParser(prog=validStat.__name__,
                        description=validStat.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('statpath', 
            help='the path of stat file')
    pReq.add_argument('outfile', 
            help='output file')
    pOpt.add_argument('-f', '--format', choices=['1', '2'], 
            default='1',
            help='the format of table, {1|multi line; 2|one line} [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
        
    
    args = p.parse_args(args)

    def import_stat(statfile):
        logging.debug('Read `{}`'.format(statfile))
        df = pd.read_csv(statfile, 
                            sep='\t',
                            header=None, 
                            index_col = 0,
                            names=['item', 'count'],
                            usecols=[0, 1],
                            comment="#")
        return df


    def get_statfile(path):
        files = glob.glob(path + "/*stat")
        db = {}
        for stat in files:
            name = op.basename(stat).rsplit('.', 1)[-1]
            if name not in db:
                db[name] = import_stat(stat)
            else:
                db[name + "2"] = import_stat(stat)
        
        return Namespace(**db)

    def main(path, out, fmt='1'):
        statfiles = get_statfile(path)
        mRS = statfiles.mRSstat
        mpair = statfiles.mpairstat
        mmap = statfiles.mmapstat
        mmap = statfiles.mmapstat2
        merge = statfiles.mergestat
        
        db = OrderedDict()
        if fmt == '1':

            db["Statistics of mapping"] = ""
            db["Clean Paired-end Reads"] = "{:,}".format(mpair.loc['Total_pairs_processed']['count'])
            db["Unmapped Paired-end Reads"] = "{:,}".format(mpair.loc['Unmapped_pairs']['count'])
            db['Unmapped Paired-end Reads Rate (%)'] = "{:.2f}".format(1.0 * mpair.loc['Unmapped_pairs']['count']
                                                                    / mpair.loc['Total_pairs_processed']['count']
                                                                    * 100.0)
            db['Paired-end Reads with Singleton'] = "{:,}".format(mpair.loc['Pairs_with_singleton']['count'])
            db['Paired-end Reads with Singleton Rate (%)'] = "{:.2f}".format(1.0 * mpair.loc['Pairs_with_singleton']['count']
                                                                    / mpair.loc['Total_pairs_processed']['count'] * 100
                                                                )
            db['Multi Mapped Paired-end Reads'] = "{:,}".format(mpair.loc['Multiple_pairs_alignments']['count'])
            db['Multi Mapped Rate (%)'] = "{:.2f}".format(1.0 * mpair.loc['Multiple_pairs_alignments']['count']
                                                                    / mpair.loc['Total_pairs_processed']['count'] * 100
                                                                )
            db['Low Mapped Quality Reads'] = "{:,}".format(mpair.loc['Low_qual_pairs']['count'])
            db['Low Quality Mapped Rate (%)'] = "{:.2f}".format(1.0 * mpair.loc['Low_qual_pairs']['count']
                                                                    / mpair.loc['Total_pairs_processed']['count'] * 100)
            db['Unique Mapped Paired-end Reads'] = "{:,}".format(mpair.loc['Unique_paired_alignments']['count'])
            db['Unique Mapped Rate (%)'] =  "{:.2f}".format(1.0 * mpair.loc['Unique_paired_alignments']['count']
                                                                    / mpair.loc['Total_pairs_processed']['count'] * 100)
            db['Statistics of valid reads'] = ""
            db['Unique Mapped Paired-end Reads'] = "{:,}".format(mpair.loc['Unique_paired_alignments']['count'])
            db['Dangling End Paired-end Reads'] = "{:,}".format(mRS.loc['Dangling_end_pairs']['count'])
            db['Dangling End Rate (%)'] = "{:.2f}".format(1.0 * mRS.loc['Dangling_end_pairs']['count'] /
                                                        mpair.loc['Unique_paired_alignments']['count'] * 100)
            db['Self Circle Paired-end Reads'] = "{:,}".format(mRS.loc['Self_Cycle_pairs']['count'])
            db['Self Circle Rate (%)'] = "{:.2f}".format(1.0 * mRS.loc['Self_Cycle_pairs']['count'] /
                                                        mpair.loc['Unique_paired_alignments']['count'] * 100)
            db['Dumped Paired-end Reads'] = "{:,}".format(mRS.loc['Dumped_pairs']['count'])
            db['Dumped Rate (%)'] = "{:.2f}".format(1.0 * mRS.loc['Dumped_pairs']['count'] /
                                                        mpair.loc['Unique_paired_alignments']['count'] * 100)
            
            db['Interaction Paired-end Reads'] = '{:,}'.format(merge.loc['valid_interaction']['count'])
            db['Interaction Rate (%)'] = '{:.2f}'.format(1.0 * merge.loc['valid_interaction']['count'] /
                                                        mpair.loc['Unique_paired_alignments']['count'] * 100)
            db['Lib Valid Paired-end Reads'] = '{:,}'.format(merge.loc['valid_interaction_rmdup']['count'])
            db['Lib Valid Rate (%)']= '{:.2f}'.format(1.0 * merge.loc['valid_interaction_rmdup']['count'] /
                                                        merge.loc['valid_interaction']['count'] * 100)
            db['Lib Dup Rate (%)'] = '{:.2f}'.format(100- (1.0 * merge.loc['valid_interaction_rmdup']['count'] /
                                                        merge.loc['valid_interaction']['count'] * 100))
        elif fmt == '2':
            db['Raw reads'] = '{:,}'.format(mpair.loc['Total_pairs_processed']['count'])
            db['Mapped pairs'] = '{:,}'.format(mpair.loc['Total_pairs_processed']['count'] - 
                                                    mpair.loc['Unmapped_pairs']['count'])
            db['Unique pairs'] = '{:,}'.format(mpair.loc['Unique_paired_alignments']['count'])
            db['Self-circle'] = '{:,}'.format(mRS.loc['Self_Cycle_pairs']['count'])
            db['Dangling'] = '{:,}'.format(mRS.loc['Dangling_end_pairs']['count'])
            db['PCR duplicate'] = '{:,}'.format(merge.loc['valid_interaction']['count'] -
                                                merge.loc['valid_interaction_rmdup']['count'])
            db['Valid contact'] = '{:,}'.format(merge.loc['valid_interaction_rmdup']['count'])


        df = pd.DataFrame([db])
        df = df.T if fmt == '1' else df
        header = None if fmt == '1' else True
        index = True if fmt == '1' else None 
        df.to_excel(out + ".xls", header=header, index=index)
        df.to_csv(out, sep='\t', header=header, index=index)
        logging.debug('Output file to `{}`'.format(out))


        
    main(args.statpath, args.outfile, args.format)
    
    
    

def plotDistDensity(args):
    """
    %progs all.validpairs out [Options]
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
    p.add_option('--plotSlope', action='store_true', default=False,
            help='plotting slope line in picture [default: %default]')
    p.add_option('--slopeRange', default='500000-7000000', 
            help='slope range [default: %default]')
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
            xmin=opts.xmin, xmax=opts.xmax, plotSlope=opts.plotSlope, 
            slopeRange=opts.slopeRange) #, color=color)
def plotIDEMulti(args):
    """
    %(prog)s 1.ValidPairs 2.ValidPairs ... [Options]
        To multi sample IDE in a picture.
    """
    p = argparse.ArgumentParser(prog=plotIDEMulti.__name__,
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
            help='the scale of data [default: %(default)s]')
    pOpt.add_argument('--xmin', default=1e5, type=float, metavar='float',
            help='min value of xtick [default: %(default)s]')
    pOpt.add_argument('--xmax', default=2e7, type=float, metavar='float',
            help='max value of xtick [default: %(default)s]')
    pOpt.add_argument('--plotSlope', action='store_true', default=False,
            help='plotting slope line in picture [default: %(default)s]')
    pOpt.add_argument('--slopeRange', default='500000-7000000', 
            help='slope range [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    from scipy.stats import linregress
    from matplotlib.lines import Line2D
    scale = args.scale
    xmin = args.xmin
    xmax = args.xmax
    plotSlope = args.plotSlope
    slope_range = args.slopeRange
    slope_left, slope_right = list(map(int, slope_range.split("-")))
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
   
        
        slope_data = data[(data >= slope_left) & (data <= slope_right)]
        slope_unique, slope_counts = np.unique(slope_data, return_counts=True)
        slope_unique_log, slope_counts_log = np.log10(slope_unique), np.log10(slope_counts)
        slope, intercept, rvalue, pvalue, stderr = linregress(slope_unique_log, 
                                                                slope_counts_log)
        if plotSlope:
            label = "{}".format(label)
        else:
            label = "{} ({:.2f})".format(label, slope)
        plt.plot(list(db.keys()), list(db.values()), label=label)
        if plotSlope:
            slope_x_values = slope_unique
            slope_y_values = 10 ** (slope_unique_log * slope + intercept)
            slope_label = "{:.2f}".format(slope)
            plt.plot(slope_x_values, slope_y_values, "--", label=slope_label, lw=2)
        #sns.regplot(list(db.keys()), list(db.values()), label=label, 
         #   marker=Line2D.filled_markers[i], ci=0, truncate=True,
        #    )
        #i += 1
        

        
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

def plotIDEMultiv1(args):
    """
    %(prog) 1.ValidPairs 2.ValidPairs ... [Options]
        To multi sample IDE in a picture.
    """
    p = argparse.ArgumentParser(prog=plotIDEMulti.__name__,
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
    p = argparse.ArgumentParser(prog=plotCisTrans.__name__,
                        description=plotCisTrans.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cool', nargs='+', 
            help='cool file of hicmatrix')
    pOpt.add_argument('--sample', nargs='+',
            help='sample name of each cool file [default: coolprefix]', )
    pOpt.add_argument('--perchr', action='store_true', default=False,
            help='plot cis/trans ration per chromosome [default: %(default)]')
    pOpt.add_argument('--hideLegend', action='store_true', default=False, 
            help='hide legend label [default: %(default)s]')
    pOpt.add_argument('--colors', nargs="*", default=None,
            help='colors of bar [default: %(default)s]')
    pOpt.add_argument('-o', '--outprefix', default='out',
            help='out prefix of output picture and result [default: %(default)]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    outprefix = args.outprefix
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
    if len(cis_list) == 1 and len(trans_list) == 1:
        cis, trans = cis_list[0], trans_list[0]
        fig, ax = plt.subplots(figsize=(2.8, 4))
        ax.bar([1.1, 2], [cis, trans], width=.6, align='center',
            color=['#ffa040', '#55a0fb'], edgecolor='#606060',
            linewidth=2)
        ax.set_xticks([1, 2])
        ax.set_xlim(0.5, 2.5)
        ax.set_ylim(0, max([cis, trans]) * 1.2)
        ax.set_xticklabels(['${cis}$', '${trans}$'], fontsize=14)
        ax.set_ylabel('Contact probability', fontsize=15)
        ax.tick_params(axis='both', which='major', labelsize=14,)
        ax.tick_params(axis='x', which='major', width=0)
        ax.tick_params(axis='y', which='major', width=2, length=5)
        ax.spines['left'].set_linewidth(2)
        ax.spines['bottom'].set_linewidth(0)
        ax.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
        sns.despine()
        outpdf = outprefix + "_cis_trans.pdf"
        outpng = outprefix + "_cis_trans.png"
        plt.savefig(outpdf, bbox_inches='tight', dpi=300)
        plt.savefig(outpng, bbox_inches='tight', dpi=300)
        logging.debug("Successful plotting total cis trans barplot in `{}`".format(outpdf))
    output = outprefix + "_cis_trans.tsv"
    ## output total cis and trans counts 
    with open(output, 'w') as out:
        print("Cis", end='\t', file=out)
        print('\t'.join(map(str, cis_list)), file=out)
        print("Trans", end='\t', file=out)
        print('\t'.join(map(str, trans_list)), file=out)
        logging.debug("Output total cis and trans count to `{}`".format(output))
    plt.close() ## close plt
    ## plotting per chromosome cis/trans ratio barplot
    if args.perchr:
        if not args.sample:
            samples = list(map(lambda x: x.split(".")[0], args.cool))
        else:
            samples = args.sample
            if len(samples) != len(args.cool):
                logging.error('cool must equal sample')
        legend = False if args.hideLegend else True
        colors = args.colors if args.colors else None
        data = concatCisTransData(args.cool, samples)
        outtsv = outprefix + "_cis_trans_ratio_perchr.tsv"
        data.to_csv(outtsv, sep='\t', header=None, index=None)
        logging.debug("Successful outputting cis/trans ratio data to `{}`".format(outtsv))
        ax = plotCisTransPerChrom(data, legend=legend)
        outpdf = outprefix + "_cis_trans_ratio_perchr.pdf"
        outpng = outprefix + "_cis_trans_ratio_perchr.png"
        plt.savefig(outpdf, bbox_inches='tight', dpi=300)
        plt.savefig(outpng, bbox_inches='tight', dpi=300)
        logging.debug("Successful plotting cis/trans ratio per chrom barplot in `{}`".format(outpdf))


def statFrag(args):
    """
    %(prog)s allValidParis [Options]
        stat the Ratio of theoretically digested genomic 
        fragments covered by valid paired Hi-C reads.
    """
    p = argparse.ArgumentParser(prog=statFrag.__name__,
                        description=statFrag.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('validpairs',  help='Validparis file')
    pReq.add_argument('enzyme', help='restriction enzyme site bed file')
    pOpt.add_argument('--unmap', default=False, action='store_true',
            help='output the unmap fragments bed [default: %(default)s]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
            help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    theo_num = 0
    reality_num = 0
    check_file_exists(args.enzyme)

    enzyme_df = pd.read_csv(args.enzyme, sep='\t', 
                    header=None, index_col='name',
                    names=['chrom', 'start', 'end',
                            'name', 'flag', 'strand'])
    theoFrags = set(enzyme_df.index)
    theo_num = len(theoFrags)
    vp = ValidPairs(args.validpairs)
    realFrags = vp.getRealFrags()
    reality_num = len(realFrags)
    
    if args.unmap:
        unmapFrags = theoFrags - realFrags 
        unmap_df = enzyme_df.loc[unmapFrags]
        unmap_df.reset_index(inplace=True)
        df = unmap_df[['chrom', 'start', 'end',
                            'name', 'flag', 'strand']]
        df = df.sort_values(by=['chrom', 'start'])
        df.to_csv(args.output.name + ".unmap", sep='\t', header=None, index=None)

    print("Theoretical Fragments\t{}".format(theo_num), file=args.output)
    print("Reality Fragments\t{}".format(reality_num), file=args.output)
    print("Reality Fragments Ratio (%)\t{:.2%}".format(reality_num*1.0/theo_num), 
            file=args.output)

def getRealityBed(args):
    """
    %(prog)s <allValidParis> <theoretical.bed> [Options] > <reality.bed> 
        get reality fragment bed file from valid pairs file and 
            theoretical bed.
    """
    p = argparse.ArgumentParser(prog=getRealityBed.__name__,
                        description=getRealityBed.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('validpairs',  help='Validparis file')
    pReq.add_argument('enzyme', help='restriction enzyme site bed file')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'), default=sys.stdout,
            help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    check_file_exists(args.enzyme)

    enzyme_df = pd.read_csv(args.enzyme, sep='\t', 
                    header=None, index_col='name',
                    names=['chrom', 'start', 'end',
                            'name', 'flag', 'strand'])

    vp = ValidPairs(args.validpairs)
    realFrags = vp.getRealFrags()
    df = enzyme_df.loc[realFrags]
    df.reset_index(inplace=True, )
    
    df = df[['chrom', 'start', 'end',
                            'name', 'flag', 'strand']]
    df = df.sort_values(by=['chrom', 'start'])
    df.to_csv(args.output, sep='\t', header=None, index=None)


def plotFragmentSizeHist(args):
    """
    %(prog)s <fragment1 ...> [Options]
        plot histogram of hic validpairs fragment size
    """

    p = argparse.ArgumentParser(prog=plotFragmentSizeHist.__name__,
                        description=plotFragmentSizeHist.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fragsize',  nargs="+",
            help='fragment size')
    pOpt.add_argument('--xMax', default=1500, type=int,
            help='max value of x-axis [default: %(default)s]')
    pOpt.add_argument('--bin', default=25, type=int,
            help='bin size of histogram [default: %(default)s]')
    pOpt.add_argument('-o', '--output', default='plotFragmentSizeHist.png',
            help='output of picture')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    def import_fragment_size(sizefile, sample):
        df = pd.read_csv(sizefile, header=None, index_col=None,
                        sep='\t', names=['size'], dtype='int32')
        df['size'] = df['size'].astype('int32')
        df['sample'] = sample
        df['sample'] = df['sample'].astype('category')

        return df

    samples = list(map(lambda x: x.rsplit(".")[0], args.fragsize))
    df_list = []
    for fragsize, sample in list(zip(args.fragsize, samples)):
        tmp_df = import_fragment_size(fragsize, sample)
        df_list.append(tmp_df)
    df = pd.concat(df_list, axis=0)

    ## plot
    sns.displot(df, x='size', hue='sample', stat='probability',#stat='density', 
            bins=range(0, args.xMax, args.bin), common_norm=False)
    ax = plt.gca()
    ax.set_xlabel('fragment size (bp)', fontsize=14)
    ax.set_ylabel('probability', fontsize=14)

    sns.despine(trim=True)

    plt.savefig(args.output, dpi=300)
    

    

if __name__ == "__main__":
    main()
