#!/usr/bin/env python2
# -*- coding:utf-8 -*-

"""
A/B compartments analysis libraries
"""

from __future__ import print_function

import argparse
import logging
import sys
import os
import os.path as op
import numpy as np
import pandas as pd

import cooler
import math
import scipy
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from collections import defaultdict, OrderedDict
from joblib import delayed, Memory, Parallel
from TDGP.analysis.genome import Genome
from TDGP.apps.mathutils import completeIC, observedOverExpected
from TDGP.apps.base import debug, check_file_exists, listify
from TDGP.apps.base import ActionDispatcher
from TDGP.apps.utilities import isCooler
from TDGP.formats.bedGraph import BedGraph
from TDGP.formats.hicmatrix import cool2matrix
from TDGP.graphics.ploty import change_width
debug()
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def main():

    actions = (
            ("test", "test"),
            ("plotLineRegress", "Plot two species pca1 lineregress"),
            ('plotMultiLineRegress', "plot two species pca1 lineregress per switch type"),
            ('plotLRPerChrom', 'plot two samples pca1 linregress per chromosome'),
            ('plotEnrichment', 'plot compartment strength enrichment'),
            ('plotStrength', 'plot compartment strength in multi samples'),
            ('plotSwitchPie', 'plot two samples switch type pie picture'),
            ('plotBoxPerChrom', 'plot boxplot of some data per chromosomes'),
            ('plotBoxMultiSamples', 'plot boxplot of some data per samples on one chromosome'),
            ('quickPlot', 'quick plot A/B compartments to show pc1 and gene density'),
            ('getSyntenyGenePca', "get the synteny gene pairs pca value"),
            ('annotateSwitchType', "annotate swithch type for synteny gene pairs"),
            ('annotateType', 'annotate the compartment type for a pca bg file'),
            ('statAB', 'stat A/B compartments informations'),
            ('buscoGeneDist', 'to analysis busco gene distribution between A and B'),
            ('getSwitchLink', 'to get links between sample1 and sample2 to plot circos')
            
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def test(args):
    """
    %(prog)s 
    """
    import argparse
    from argparse import ArgumentParser

    p = ArgumentParser(prog=test.__name__, description=test.__doc__,
            conflict_handler='resolve', formatter_class=argparse.RawDescriptionHelpFormatter)
    req = p.add_argument_group('Required arguments')
    req.add_argument('infile', help='infile', nargs=2)
    opt = p.add_argument_group('Optional arguments')
    
    opt.add_argument('-p', '--pi', default=1, nargs=2, required=False)
    opt.add_argument('--help', action='help')
    args = p.parse_args(args)
    print(args.infile, args.pi)


class Compartment(object):
    """
    Compartment analysis library

    Params:
    --------

    Returns:
    --------

    Examples:
    --------
    """
    def __init__(self, mem_cache='.'):
        memory = Memory(mem_cache, verbose=0)
        self.getStrengthError = memory.cache(self._getStrengthError)
    

    @staticmethod
    def calStrength(strength, permutted):
        def substrength(i):
            # comment return below to make kind of weird but realistic 50/50 compartment split 
            # rather than default 20% compartments
            # ------------!!!!!!!!!!!!!!!!!!!!!!-------------------
            return np.log(i[0,0] * i[-1, -1] / i[0, -1]**2)
            t1 = i[:3, :3].sum() + i[:2, :2].sum()
            t2 = i[2:, 2:].sum() + i[3:, 3:].sum()
            t3 = i[:3, 2:].sum() + i[:2, 3:].sum()
            return np.log(t1 * t2 / t3**2)
        s1 = substrength(strength)
        sp = [substrength(i) for i in permutted]
        return s1, np.std(sp)

    
    
    def _getStrengthError(self, data, eig, genome, resolution, chr=[], 
            correct=False, iterCorrect=True, thread=24):
        """
        Calculate compartments strength.
            reference: Ilya M.Flyamer et.al., Nature, 2019.
            code refrence: https://bitbucket.org/mirnylab/hiclib/src/default 
                    /examples/singleCellScripts/singleShared.py
        """
        
        
        chr = chr if chr else genome.chromLabels
        if eig == 'GC':
                _eig = genome.getGCBin(window=resolution, 
                    chr=chr, correct=True, thread=thread)
        else:
            _eig = eig

        strengthes = []
        permutted = []
        strength = np.zeros((5, 5), dtype=float)
        for i in range(100):
            permutted.append(np.zeros((5, 5), dtype=float))
        
        for n, chrom in enumerate(chr):
            idx = genome.label2idx[chrom]
            if correct:
                data[chrom] = completeIC(data[chrom])
            start, end = genome.chromBinsDict[chrom]    
            
            cmatrix = data[chrom]
            #cmatrix[np.where(cmatrix == 0)] = 1.0 # zero bin mask with 1.0
            np.savetxt('completeIC.{}.txt'.format(chrom), cmatrix, delimiter='\t')
            cmatrix = observedOverExpected(cmatrix)
            np.savetxt('obs_exp.{}.txt'.format(chrom), cmatrix, delimiter='\t')  
            mask = np.sum(cmatrix, axis=0) > 0
            cmatrix = cmatrix[mask]
            cmatrix = cmatrix[:, mask]

            
            gc = np.array(_eig[genome.label2idx[chrom]]) if eig == 'GC' \
                else np.array(_eig[chrom])
            gc = gc[mask]
            if len(gc) > 5:
                for i in range(5):
                    for j in range(5):
                        g1, g2 = np.percentile(gc, [20*i, 20*i + 20])
                        mask1 = (gc > g1) * (gc < g2)
                        g1, g2 = np.percentile(gc, [20*j, 20*j + 20])
                        mask2 = (gc > g1) * (gc < g2)

                        addition = cmatrix[np.ix_(mask1, mask2)]
                        if iterCorrect:
                            addition = np.reshape(addition, (-1))
                            for k in range(100):
                                reshampled = np.random.choice(addition, 
                                    len(addition), replace=True)
                                permutted[k][i, j] += reshampled.mean()
                        
                        strength[i, j] += addition.mean()
            if not iterCorrect:
                strengthes.append(strength)
                strength = np.zeros((5,5), dtype = float)
        
        strengthes = strength if iterCorrect else strengthes
        return strengthes, permutted   
    
class ABComparision(object):
    """
    A/B compartments comparision across same species.
    
    Parameters
    ----------
    object : [type]
        [description]
    """
    def __init__(self):
        self.stype = ('AA', 'BB', 'BA', 'AB')

    @staticmethod
    def annotateSwitchType(infile, out, 
            stype=('AA', 'BB', 'BA', 'AB')):
        """
        To annotate switch type for compartments pairs 

            ** infile format:
                chrom1 start1 end1 value1 chrom2 start2 end2 value2
        """
        check_file_exists(infile)
        with open(infile, 'r') as fp:
            for line in fp:
                line_list = line.strip().split()
                v1 = line_list[3]
                v2 = line_list[7]
                type_ = switch_type(float(v1), float(v2))
                if type_ in listify(stype):
                    line_list.append(type_)
                    print("\t".join(line_list), file=out)
        logging.debug('Successful ... result is in `{}`'.format(out.name)) 

class ABComparisionSpecies(object):
    """
    A/B compartments comparision across two species.
    
    Parameters
    ----------
    object : [type]
        [description]
    """
    def __init__(self):
        self.stype = ('AA', 'BB', 'BA', 'AB')
    
    @staticmethod
    def getSyntenyGenePca(bg1, bg2, bed1, bed2, plot=True):
        """
        To get ab PC1 value of per synteny gene.

        Params:
        ---------
        bg1, bg2, bed1, bed2: `str` infile

        Returns:
        ---------
        out: some output file

        Examples:
        ---------
        >>> ABC = ABComparisionSpecies()
        >>> ABC.getSyntenyGenePca(bg1, bg2, bed1, bed2)
        """
        list(map(check_file_exists, (bg1, bg2, bed1, bed2)))
        bedtools_formatter = "bedtools intersect -a {} -b {} -wao -f 0.5  |cut -f 1-4,8 > {}\n"
        out_formatter = "{}.synteny.eigen1.bg"
        cut_pca_formatter = "cut -f 5 {0}.synteny.eigen1.bg > {0}.pca1 \n"
        prefix1 = bed1.split('.')[0]
        prefix2 = bed2.split('.')[0]

        bedtools_cmd1 = bedtools_formatter.format(bed1, bg1, 
                            out_formatter.format(prefix1))
        bedtools_cmd2 = bedtools_formatter.format(bed2, bg2, 
                            out_formatter.format(prefix2))
        cut_cmd1 = cut_pca_formatter.format(prefix1)
        cut_cmd2 = cut_pca_formatter.format(prefix2)
        paste_cmd = "paste {0}.synteny.eigen1.bg {1}.synteny.eigen1.bg > {0}-{1}.synteny.eigen1.bg\n".format(prefix1, 
                                                                            prefix2)
        list(map(os.system, (bedtools_cmd1, bedtools_cmd2, cut_cmd1, cut_cmd2, paste_cmd)))
        logging.debug('Successful to getSyntenyGenePca')
        
        synteny_pairs = '{0}-{1}.synteny.eigen1.bg'.format(prefix1, prefix2)

        ABC = ABComparisionSpecies()
        for type_ in ABC.stype:
            out = '{0}-{1}.{2}.synteny.eigen1.bg'.format(prefix1, prefix2, type_)
            annotateSwitchType([synteny_pairs,
                                "--stype", type_,
                                "-o", out])
            cmd1 = 'cut -f 5 {} > {}.{}.pca1'.format(out, prefix1, type_)
            cmd2 = 'cut -f 10 {} > {}.{}.pca1'.format(out, prefix2, type_)
            os.system(cmd1)
            os.system(cmd2)

            if plot:
                plotLineRegress(["{}.{}.pca1".format(prefix1, type_),
                            "{}.{}.pca1".format(prefix2, type_),
                            "--xlabel", prefix1,
                            "--ylabel", prefix2,
                            "-o", 
                            "{}-{}.{}.lineregress.pdf".format(
                                prefix1, prefix2, type_)])
        if plot:
            
            plotLineRegress(["{}.pca1".format(prefix1),
                            "{}.pca1".format(prefix2),
                            "--xlabel", prefix1,
                            "--ylabel", prefix2,])
            

    @staticmethod
    def annotateSwitchType(infile, out, 
            stype=('AA', 'BB', 'BA', 'AB')):
        """
        To annotate switch type for synteny gene pairs 
             ** infile format:
                chrom1 start1 end1 value1 chrom2 start2 end2 value2
        """
        check_file_exists(infile)
        with open(infile, 'r') as fp:
            for line in fp:
                line_list = line.strip().split()
                v1 = line_list[4]
                v2 = line_list[9]
                type_ = switch_type(float(v1), float(v2))
                if type_ in listify(stype):
                    line_list.append(type_)
                    print("\t".join(line_list), file=out)
        logging.debug('Successful ... result is in `{}`'.format(out.name))    

    

    @staticmethod
    def plotLineRegress(a, b, out=None, xlabel='species1', ylabel='species2'):
        """
        To plot the line regress of two list pca1.

        Params:
        --------
        a: `array-like`
        b: `array-like`

        Returns:
        --------
        picture

        Examples:
        ---------
        >>> a = [-0.1, -0.2, 0.3 ...]
        >>> b = [0.2, 0.5, -0.1 ...]
        >>> ABC = ABComparisionSpecies()
        >>> ABC.plotLineRegress(a, b)
        
        """
        from matplotlib.lines import Line2D
        from scipy.stats import linregress
        from seaborn import regplot
        assert len(a) == len(b), "`a` and `b` length must euqual"
        if not out:
            out = "{}-{}_ab_lineregress.pdf".format(xlabel, ylabel)
        scatter_params = dict(color='#209093', s=2)
        line_params = dict(color='#032F49', lw=2)
        slope, intercept, rvalue, pvalue, stderr = linregress(a, b)
        #r2 = rvalue ** 2 if rvalue > 0 else -1 * rvalue ** 2
        label = [r"r = {:.2f}  $\mathit{{P}}$ = {:.2e}".format(rvalue, pvalue)]
        fig, ax = plt.subplots(figsize=(5, 5))
        regplot(a, b, ax=ax, truncate=True, 
                scatter_kws=scatter_params, line_kws=line_params)
        legend_elements = [Line2D([0], [0], **line_params)]
        #ax.set_title('The regression of PC1')
        ax.set_xlabel("{}".format(xlabel), fontsize=12)
        ax.set_ylabel("{}".format(ylabel), fontsize=12)
        plt.legend(legend_elements, label, loc='best')
                #bbox_to_anchor=(1, 0.5))
        plt.savefig(out, dpi=300, bbox_inches='tight')
        logging.debug('Successful, picture is in `{}`'.format(out))

    @staticmethod
    def plotMultiLineRegress(annotated_file, out=None, xlabel='species1', 
            ylabel='species2'):
        """
        To plot four switch type line regression plot in a picture
            through seaborn lmplot.
        """
        from matplotlib.lines import Line2D
        from matplotlib.ticker import MaxNLocator
        mylocator = MaxNLocator(6) # set max locator

        out = '{}-{}.multiLineRegress.pdf'.format(xlabel, ylabel) \
                    if not out else out
        columns = ('chr1', 'start1', 'end1', 'gene1', 'value1',
                    'chr2', 'start2', 'end2', 'gene2', 'value2',
                    'type')
        hue_order = ('AA','BB','AB','BA')
        color_palette = ('#BB4853','#209093','#032F49','#797596')
        df = pd.read_csv(annotated_file, sep='\t', header=None,
                        names=columns)
        
        ## calculate rvalue by scipy.stats.lineregress
        def calc(stype):
            a = df.loc[df.type == stype]['value1']
            b = df.loc[df.type == stype]['value2']

            slope, intercept, rvalue, pvalue, stderr = scipy.stats.linregress(a, b)
            text = '{} (r={:.2f})'.format(stype, rvalue)
            return text
        
        rvalues = list(map(calc, hue_order))
        fig, ax = plt.subplots(figsize=(7, 6))
        #plt.subplots_adjust(left=0, bottom=0, right=1, top=1, wspace=0, hspace=0)
        with sns.plotting_context(rc={'legend.fontsize': 7}):
            g = sns.lmplot(x='value1', y='value2', hue='type', truncate=True,
                            height=5, data=df, hue_order=hue_order, 
                            palette=color_palette, scatter_kws={'s': 2})
            g.set_axis_labels(xlabel, ylabel)
            ax = plt.gca() 
            ax.xaxis.set_major_locator(mylocator) # set xaxis major locator to 6
            #legend_elements = ([Line2D([0], [0], color=i, lw=2) for i in color_palette])
            #legend_labels = list(map(calc, hue_order))
            #ax.legend(legend_elements, legend_labels, bbox_to_anchor=(1.05, 0.5), loc=2)
            g._legend.set_title('Switch type') 
            g._legend.set_bbox_to_anchor([1.05, 0.5])
            
            #g._legend.fontsize = 4
            #g._legend._loc = 2
            for t, l in zip(g._legend.texts, rvalues):
                t.set_text(l)
            
            plt.savefig(out, dpi=300, bbox_inches='tight')


def import_ab_data(path):
    """
    Import data from bed4 as `pandas.DataFrame`, 
    and name columns with `chrom    start    end    pca1`.
    """
    df = pd.read_csv(path, 
                    sep='\t',
                    header=None, 
                    names=[ 'chrom', 'start', 'end', 'pca1'], 
                    index_col=0)
    return df


def ab_dist(df):
    """
    Caclculate the distribution of A/B compartments in chromosome and whole chromosomes.
    """
    dist_db = defaultdict(lambda :[0, 0])
    for chrom in sorted(set(df.index)):
        df_chrom = df.loc[chrom]
        chrom_pca1 = df_chrom.pca1
        distances = df_chrom.end - df_chrom.start
        distances = distances.tolist()
        for n, pca1 in enumerate(chrom_pca1):
            if pca1 > 0:
                dist_db[chrom][0] += distances[n]
            elif pca1 < 0:
                dist_db[chrom][1] += distances[n]
        result_df = pd.DataFrame(dist_db).T
        result_df.loc['Total'] = result_df.sum()
    return result_df


def is_conserved(v1, v2):
    """
    judge the conservative of ab compartments in two species
    >>> is_conserved(-0.1, -0.2)
    True
    >>> is_conserved(-0.1, 0.1)
    False
    """
    if v1 * v2 <= 0:
        return False
    if v1 *v2 > 0:
        return True

    
def switch_type(v1, v2):
    """
    judge the switch type of ab compartments in two species: `AA`, `BB`, `AB`, `BA`
    >>> switch_type(-0.1, 0.2)
    'BA'
    >>> switch_type(0.1, 0.1)
    'AA'
    """
    if v1 > 0:
        if v2 > 0:
            return 'AA'
        else:
            return 'AB'
    else:
        if v2 > 0:
            return 'BA'
        else:
            return "BB"

def switch_type_in_df(df):
    """
    judge the switch type of ab compartments in two dataframe of pca1 bedGraph
    >>> df1-df2.head(1)
    chrom1 start1 end1 value1 chrom2 start2 end2 value2

    >>> df1-df2.apply(switch_type_in_df, axis=1)
    """
    if df.value1 >= 0 and df.value2 >= 0:
        return 'AA'
    elif df.value1 >= 0 and df.value2 < 0:
        return 'AB'
    elif df.value1 < 0 and df.value2 < 0:
        return 'BB'
    elif df.value1 < 0 and df.value2 >= 0:
        return 'BA'

def two_species_conserved_compare(tgy_df, jgy_df):
    """
    compare two file ab compartments  switch type and return a directory.
    """
    db = {"AA": 0, "BB": 0, "AB": 0, "BA": 0}
    length_df = tgy_df.end - tgy_df.start
    for i in range(len(tgy_df)):
        v1, v2 = tgy_df.pca1.iloc[i], jgy_df.pca1.iloc[i]
        db[switch_type(v1, v2)] += length_df.iloc[i]
    
    return db


def cut_chrom_range(data, chrom, start=None, end=None):
    """
    Get chrom range from data.
    """
    data_chrom = data.loc[chrom]
    if start or end:
        tmp_data = data_chrom[ (data_chrom.end <= end) & (data_chrom.end >= start)]
    else:
        tmp_data =  data_chrom
    return tmp_data

    
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




def plot_ab(ax, ab_data, xaxis_pos='bottom'):
    """
    Plot A/B compartments
    """
    xdata = ab_data.end
    ax.fill_between(xdata, ab_data.pca1, where=ab_data.pca1>0, facecolor='#BB4853')
    ax.fill_between(xdata, ab_data.pca1, where=ab_data.pca1<0, facecolor='#209093')
    
    ax.set_xticks(np.linspace(xdata.iloc[0], xdata.iloc[-1], 8))
    ax.set_xticklabels(chrom_ticks_convert(np.linspace(xdata.iloc[0], xdata.iloc[-1], 8) ))
    
    tb_spines = ['bottom', 'top']
    tb_spines.remove(xaxis_pos)
    tb_spine, = tb_spines
    for pos in [tb_spine, 'right']:
        ax.spines[pos].set_visible(False)

    if xaxis_pos == 'top':
        ax.xaxis.tick_top()
        ax.xaxis.set_label_position('top')
    ax.spines['left'].set_bounds(-0.05, 0.05)
    ax.yaxis.set_major_locator(plt.MaxNLocator(3))
    ax.spines[xaxis_pos].set_bounds(xdata.iloc[0], xdata.iloc[-1])
    ax.set_xlabel('{}'.format(ab_data.index[0]))
    
    return ax



def two_species_ab_compare_plot(tgy_df, jgy_df, chrom,
                           start=None, end=None):
    """
    Plot
    """

    fig, axs = plt.subplots(2,1, figsize=(16, 4))
    (ax1, ax2) = axs

    tgy_ab_data = cut_chrom_range(tgy_df,  chrom, start, end)
    jgy_ab_data = cut_chrom_range(jgy_df, chrom, start, end)
    plot_ab(ax1, tgy_ab_data, 'top')
    plot_ab(ax2, jgy_ab_data)



def plot_two_species_conserved_pie(tgy_df, jgy_df):
    """
    Pie plot of two species conserved.
    """
    db = two_species_conserved_compare(tgy_df, jgy_df)
    conserved = db['A Stable'] + db['B Stable']
    non_conserved = db['A2B'] + db['B2A']
    conserved_label = chrom_size_convert(conserved)
    non_conserved_label = chrom_size_convert(non_conserved)
    labels = ['Conserved ({})'.format(conserved_label), 'Non-Conserved ({})'.format(non_conserved_label)]
    colors = [ '#209093', '#BB4853']
    with mpl.rc('font', size=16.0):
        plt.figure(figsize=(10,10))
        _, texts, autotexts = plt.pie((conserved, non_conserved), labels=labels, colors=colors, autopct='%1.2f%%')
        for autotext in autotexts:
            autotext.set_color('white')
            autotext.set_fontsize(20)


def plot_two_species_ab_pie(tgy_df, jgy_df, out):
    """
    Pie plot of two species A/B switch type.
    """
    db = two_species_conserved_compare(tgy_df, jgy_df)
    func = lambda x: "{} ({})".format(x[0], chrom_size_convert(x[1]))
    labels = list(map(func, db.items()))
    colors = ['#265e8a', '#032F49','#BB4853', '#209093','#a83836',]
    mpl.rc('font', size=16.0)
    plt.figure(figsize=(10,10))
    _, texts, autotexts = plt.pie(db.values(), labels=labels, colors=colors, autopct='%1.2f%%')
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontsize(20)
    plt.savefig(out, dpi=300)
    plt.savefig(out.rsplit('.', 1)[0] + '.png', dpi=300)

## outside command ##
def plotLineRegress(args):
    """
    %(prog)s a.pca1 b.pca1 [Options]\n
        To plot lineregress of two species' pca1 eigenvector\n
    
    Examples:\n
        %(prog)s a.pca1 b.pca1 --xlabel a --ylabel b\n
        - plotting per chromosome\n
        %(prog)s a.bg b.bg --xlabel a --ylabel b --perchrom\n
        - plotting marked color by chromosomes\n
        %(prog)s a.bg b.bg --xlabel a --ylabel b --withColor\n
    """
    p = argparse.ArgumentParser(prog=plotLineRegress.__name__,
                        description=plotLineRegress.__doc__, 
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('a')
    pReq.add_argument('b')
    pOpt.add_argument('--xlabel', default='species1', metavar='species1',
            help='the xlabel of picture [default: %(default)s]' )
    pOpt.add_argument('--ylabel', default='species2', metavar='species2',
            help='the ylabel of picture [default: %(default)s]' )
    pOpt.add_argument('-o', '--out', default=None,
            help='output file [default: xlabel-ylabel_ab_lineregress.pdf]')
    pOpt.add_argument('--perchrom', default=False, action='store_true',
            help='if plot per chromosome [default: %(default)s]')
    pOpt.add_argument('--withColor', default=False, action='store_true',
            help='if plot color for scatter by chromosome [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help', 
            help='show help message and exit.')

    args = p.parse_args(args)
    from TDGP.graphics.ploty import plotLineRegress as plr
    def checkFileType(infile):
        with open(infile) as fp:
            for line in fp:
                if len(line.strip().split()) == 1:
                    return 'list'
                elif len(line.strip().split()) == 4:
                    return 'bedGraph'
                else:
                    return 
                break
    def import_pca1(bgfile):
        df = pd.read_csv(bgfile, sep='\t', 
                        header=None, index_col=0, 
                        names=['chrom', 'start', 'end', 'pca1'])
        return df

    def trim_axs(axs, N):
        """
        little helper to message the axs list to have correct length...
        """
        axs = axs.flat
        for ax in axs[N:]:
            ax.remove()

        return axs[:N]

    def plotPerChrom(xdf, ydf, out=None, xlabel='species1', 
                        ylabel='species2', ncols=4):
        if not out:
            out = "{}-{}_ab_lineregress_perchrom".format(xlabel, ylabel)
        chrom_list = sorted(set(xdf.index))
        nrows = math.ceil(len(chrom_list) / ncols)
        fig, axes = plt.subplots(int(math.ceil(len(chrom_list)*1.0/ncols)),
                        ncols, figsize=(ncols*4 + 3,
                        int(math.ceil(len(chrom_list)*1.0/ncols))*4 + 1), 
                        constrained_layout=True)
        axes = trim_axs(axes, len(chrom_list))
        for ax, chrom in zip(axes, chrom_list):
            plr(ax, xdf.loc[chrom].pca1, ydf.loc[chrom].pca1, chrom, chrom)
        fig.text(0.5, -0.04, xlabel, ha='center', va='center', 
                    fontsize=16)
        fig.text(-0.02, 0.5, ylabel, ha='center', va='center', 
                    rotation='vertical', fontsize=16)
        plt.savefig(out.replace('.pdf', '') + '.png', 
                        dpi=300, bbox_inches='tight')
        plt.savefig(out.replace('.png', '') + '.pdf', 
                        dpi=300, bbox_inches='tight')

    def plotAll(xdf, ydf, out=None, xlabel='species1', ylabel='species2'):
        if not out:
            out = "{}-{}_ab_lineregress_all".format(xlabel, ylabel)
        fig, ax = plt.subplots(figsize=(5, 5))
        ax = plr(ax, xdf.pca1, ydf.pca1, xlabel, ylabel)
        plt.savefig(out.replace('.pdf', '') + '.png', 
                        dpi=300, bbox_inches='tight')
        plt.savefig(out.replace('.png', '') + '.pdf', 
                        dpi=300, bbox_inches='tight')

    def plotWithColor(xdf, ydf, 
                    out=None, 
                    onlySwitch=True, 
                    xlabel='species1', 
                    ylabel='species2'):
        if not out:
            out = "{}-{}_ab_lineregress_withColor".format(xlabel, ylabel)
        fig, ax = plt.subplots(figsize=(5, 5))
        chrom_list = sorted(set(xdf.index))
        df = pd.concat([xdf, ydf], axis=1)
        df.columns = ['start1', 'end1', 'pca11', 'start2', 'end2', 'pca12']
        for chrom in chrom_list:
            if not onlySwitch:
                ax.scatter(xdf.loc[chrom].pca1, ydf.loc[chrom].pca1, s=2, label=chrom)
            else:
                ax.scatter(xdf.loc[chrom].pca1, ydf.loc[chrom].pca1, s=2, color='#209093')
                tmp_df = df.loc[chrom]
                a2b = tmp_df.loc[(tmp_df.pca11 > 0 ) & (tmp_df.pca12 < 0)]
                b2a = tmp_df.loc[(tmp_df.pca11 < 0 ) & (tmp_df.pca12 > 0)]
                x = pd.concat([a2b.pca11, b2a.pca11])
                y = pd.concat([a2b.pca12, b2a.pca12])
                ax.scatter(x, y, s=2, label="{} ({})".format(chrom, len(x)))
            leg = ax.legend(bbox_to_anchor=(1, 0.5))

        ax = plr(ax, xdf.pca1, ydf.pca1, 
                            xlabel=xlabel, ylabel=ylabel, scatter=False)
        ax.add_artist(leg)
        plt.savefig(out.replace('.pdf', '') + '.png', 
                        dpi=300, bbox_inches='tight')
        plt.savefig(out.replace('.png', '') + '.pdf', 
                        dpi=300, bbox_inches='tight')

    if checkFileType(args.a) == 'list' and checkFileType(args.b) == 'list':
               
        ABC = ABComparisionSpecies()
        check_file_exists(args.a)
        check_file_exists(args.b)
        a = [float(i.strip()) for i in open(args.a)]
        b = [float(i.strip()) for i in open(args.b)]
        ABC.plotLineRegress(a, b, args.out, xlabel=args.xlabel, 
                ylabel=args.ylabel)
    elif checkFileType(args.a) == 'bedGraph' and checkFileType(args.b) == 'bedGraph':
        xdf = import_pca1(args.a)
        ydf = import_pca1(args.b)
        if args.perchrom:
            plotPerChrom(xdf, ydf, args.out, 
                        xlabel=args.xlabel, 
                        ylabel=args.ylabel)
        else:
            if args.withColor:
                plotWithColor(xdf, ydf, args.out, 
                        onlySwitch=True, 
                        xlabel=args.xlabel, 
                        ylabel=args.ylabel)
            else:
                plotAll(xdf, ydf, args.out, 
                        xlabel=args.xlabel, 
                        ylabel=args.ylabel)
            
    else:
        logging.error('incorrect of input file, only support for `list` and `bedGraph`')
        sys.exit()

    

def getSyntenyGenePca(args):
    """
    %(prog)s species1.eigen1.bg species2.eigen1.bg 
                species1.synteny.bed species2.synteny.bed [Options]

        the pipeline of get synteny gene ab pc1 value.
    """

    p = p=argparse.ArgumentParser(prog=getSyntenyGenePca.__name__,
                        description=getSyntenyGenePca.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')

    pReq.add_argument('bg1', help="species1 ab eigen1 bedGraph file")
    pReq.add_argument('bg2', help="species2 ab eigen1 bedGraph file")
    pReq.add_argument('bed1', help="species1 synteny bedfile")
    pReq.add_argument('bed2', help="species2 synteny bedfile")

    pOpt.add_argument('--plot', action='store_true', default=False,
            help='is plot the lineregress [default: %default]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    ABC = ABComparisionSpecies
    ABC.getSyntenyGenePca(args.bg1, args.bg2, args.bed1, args.bed2,
                            args.plot)


def annotateType(args):
    """
    %(prog)s <eigen1.bg> [options]
        To annotate AB type of a compartment pca results
        ** infile format:
            chrom start end value
    """

    p = p=argparse.ArgumentParser(prog=annotateType.__name__,
                        description=annotateType.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bg', help='Compartment PCA analysis result, '
        'bedgraph formats')
    pOpt.add_argument('-e', '--exclude', action='store_true', 
            default=False, help='exclude the pca value, [default: %(default)s]')
    pOpt.add_argument('-o', '--out', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    with open(args.bg) as fp:
        for line in fp:
            line_list = line.strip().split()
            _type = 'A' if float(line_list[3]) >= 0 else 'B'
            if args.exclude:
                line_list[3] = _type
            else:
                line_list.append(_type)
            print('\t'.join(line_list), file=args.out)

def annotateSwitchType(args):
    """
    %(prog)s <species1-species2.synteny.eigen1.tsv> [Options]
        
        To annotate the A/B switch type to a synteny gene pairs eigen file.
        ** infile format:
            chrom1 start1 end1 value1 chrom2 start2 end2 value2
    """
    p = p=argparse.ArgumentParser(prog=annotateSwitchType.__name__,
                        description=annotateSwitchType.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('infile', help='The synteny gene pairs file with eigen')
    pOpt.add_argument('--stype', nargs="*", 
            default=('AA', 'BB', 'BA', 'AB'),
            choices=('AA', 'BB', 'BA', 'AB'), 
            help='the switch type of output [default: %(default)s]'
            )
    pOpt.add_argument('-o', '--out', type=argparse.FileType('w'), 
            default=sys.stdout, help='output file [default: %(default)s]')
    pOpt.add_argument('--same', action='store_true', default=False,
            help='if same species [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)

    ABC = ABComparisionSpecies if not args.same else ABComparision
    ABC.annotateSwitchType(args.infile, args.out, args.stype)


def plotMultiLineRegress(args):
    """
    %(prog)s <in.annotated.tsv> [Options]
        To plot four switch type line regression plot in a picture

    """
    p = p=argparse.ArgumentParser(prog=plotMultiLineRegress.__name__,
                        description=plotMultiLineRegress.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('tsv', help='The annotated syntenygene pairs pca1 tsv')
    pOpt.add_argument('-x', '--xlabel', default='species1', 
            help='the xlabel of picture [default: %(default)s]')
    pOpt.add_argument('-y', '--ylabel', default='species2', 
            help='the ylabel of picture [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)
    ABC = ABComparisionSpecies
    ABC.plotMultiLineRegress(args.tsv, out=None, xlabel=args.xlabel, 
            ylabel=args.ylabel)


def plotSwitchPie(args):
    """
    %(prog)s <infile1> <infile2> [Options]

        Plot pie of compartments switch.
        
        **infile:
            chrom start end value
    """
    p = p=argparse.ArgumentParser(prog=plotSwitchPie.__name__,
                        description=plotSwitchPie.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('infile', nargs=2, 
            help='infile of pca1')
    pReq.add_argument('-o', '--out', required=True,
            help='output file')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.') 

    args = p.parse_args(args)
    infile1, infile2 = args.infile
    df1 = import_ab_data(infile1)
    df2 = import_ab_data(infile2)
    plot_two_species_ab_pie(df1, df2, args.out)


def plotEnrichment(args):
    """
    %(prog)s [Options]
    """
    p = p=argparse.ArgumentParser(prog=plotEnrichment.__name__,
                        description=plotEnrichment.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-m', '--matrix', nargs="+", required=True,  
            help='all chromosome matrix path list')
    pReq.add_argument('-e', '--eig', nargs='+', required=True,
            help='eig files')
    pReq.add_argument('-g', '--genome', nargs='+', required=True,
            help='genome file')
    pReq.add_argument('-w', '--window', type=int, required=True,
            help='the resolution of matrix')
    pReq.add_argument('-o', '--out', required=True, 
            help='out picture file.')
    #pOpt.add_argument('-t', '--thread', type=int, default=24,
    #        help='the thread of programe [default: %(default)s]')
    pOpt.add_argument('-c', '--chrom', nargs='*', default=[],
            help='only plot these chromosomes [default: %(default)s]')
    pOpt.add_argument('--exclude', nargs="*", default=[],
            help='exclude these chromosome [default: %(default)s]')
    pOpt.add_argument('--exclude_contig', nargs='*', 
            default=['tig', 'scafflod', 'Un', 'Sy', 'Mt', 'Pt'], 
            help='exclude these chromosome if it contain these string'
                ' [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    map(check_file_exists, args.matrix)
    data_list = []
    for species in args.matrix:
        tmp = OrderedDict()

        with open(species) as fp:
            for line in fp:
                line_list = line.strip().split()
                if len(line_list) < 2:
                    continue
                chrn, m = line_list
                tmp[chrn] = np.loadtxt(m, dtype=float)
        data_list.append(tmp)
    
    map(check_file_exists, args.genome)
    genome_list = []
    if len(args.genome) == 1 and len(genome_list) < len(data_list):
        logging.debug('Using same genome file for all matrixes')
        fasta, = args.genome
        genome = Genome(fasta, exclude=args.exclude, 
                exclude_contig=args.exclude_contig)
        genome_list = [genome] * len(data_list)
    elif len(genome_list) != len(data_list):
        logging.error('genome files must be equal to '
            'matrix files or set only one')
        sys.exit()
    else:
        for fasta in args.genome:
            genome = Genome(fasta, exclude=args.exclude, 
                    exclude_contig=args.exclude_contig)
            genome_list.append(genome)
    
    for genome in genome_list:
        genome.makeWindows(args.window)

    if 'GC' in args.eig:
        eig = ['GC'] * len(data_list)
    else:
        if len(args.eig) != len(data_list):
            logging.error('eig file must be equal to matrix files')
            sys.exit()
        eig = [BedGraph(i).values for i in args.eig]
    
    def mysum(x):
        import functools
        return functools.reduce(lambda x,y:x+y, x)

    fig, axes = plt.subplots(1, len(data_list))
    for i, matrix in enumerate(data_list):
        strengthes, permutted = Compartment().getStrengthError(matrix, 
            eig[i], genome_list[i], args.window, correct=True, 
            iterCorrect=False)
        ax = axes[i] if len(data_list) > 1 else axes
        ax.imshow(np.log(mysum(strengthes)), interpolation = 'none', 
                cmap='coolwarm')
        ax.set_xticks([])
        ax.set_yticks([])
        #plt.colorbar()
    
    plt.savefig(args.out, dpi=300, bbox_inches='tight')


def plotStrength(args):
    """
    %(prog)s [Options]
    """
    p = p=argparse.ArgumentParser(prog=plotEnrichment.__name__,
                        description=plotEnrichment.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-m', '--matrix', nargs="+", required=True,  
            help='all chromosome matrix path list')
    pReq.add_argument('-e', '--eig', nargs='+', required=True,
            help='eig files')
    pReq.add_argument('-g', '--genome', nargs='+', required=True,
            help='genome file')
    pReq.add_argument('-w', '--window', type=int, required=True,
            help='the resolution of matrix')
    pReq.add_argument('-o', '--out', required=True, 
            help='out picture file.')
    #pOpt.add_argument('-t', '--thread', type=int, default=24,
    #        help='the thread of programe [default: %(default)s]')
    pOpt.add_argument('-c', '--chrom', nargs='*', default=[],
            help='only plot these chromosomes [default: %(default)s]')
    pOpt.add_argument('--exclude', nargs="*", default=[],
            help='exclude these chromosome [default: %(default)s]')
    pOpt.add_argument('--exclude_contig', nargs='*', 
            default=['tig', 'scafflod', 'Un', 'Sy', 'Mt', 'Pt'], 
            help='exclude these chromosome if it contain these string'
                ' [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    map(check_file_exists, args.matrix)
    data_list = []
    for species in args.matrix:
        tmp = OrderedDict()

        with open(species) as fp:
            for line in fp:
                line_list = line.strip().split()
                if len(line_list) < 2:
                    continue
                chrn, m = line_list
                tmp[chrn] = np.loadtxt(m, dtype=float)
        data_list.append(tmp)
    
    map(check_file_exists, args.genome)
    genome_list = []
    if len(args.genome) == 1 and len(genome_list) < len(data_list):
        logging.debug('Using same genome file for all matrixes')
        fasta, = args.genome
        genome = Genome(fasta, exclude=args.exclude, 
                exclude_contig=args.exclude_contig)
        genome_list = [genome] * len(data_list)
   
    elif len(args.genome) != len(data_list):
        logging.error('genome files must be equal to '
            'matrix files or set only one')
        sys.exit()
    else:
        for fasta in args.genome:
            genome = Genome(fasta, exclude=args.exclude, 
                    exclude_contig=args.exclude_contig)
            genome_list.append(genome)
    
    for genome in genome_list:
        genome.makeWindows(args.window)

    if 'GC' in args.eig:
        eig = ['GC'] * len(data_list)
    else:
        if len(args.eig) != len(data_list):
            logging.error('eig file must be equal to matrix files')
            sys.exit()
        eig = [BedGraph(i).values for i in args.eig]
    
    def mysum(x):
        import functools
        return functools.reduce(lambda x,y:x+y, x)
    
    
    values = []
    errors = []
    for i, matrix in enumerate(data_list):
        strengthes, permutted = Compartment().getStrengthError(matrix, 
            eig[i], genome_list[i], args.window, correct=True, 
             iterCorrect=True)
        value, error = Compartment().calStrength(strengthes, permutted)
        values.append(value)
        errors.append(error)
    
    print("values:", end="\t")
    print("\t".join(map(str, values)))
    print("errors:", end="\t")
    print("\t".join(map(str, errors)))
    names = list(map(lambda x: x.split('.')[0], args.matrix))
    index = np.arange(len(values))
    ax = sns.barplot(index, values, yerr=errors, error_kw=dict(ecolor='0.3'), 
            palette=sns.color_palette('Set2'))
    if 'GC' in args.eig :
        ymax = max(values) + 0.1
        plt.ylim(0, ymax)
    ax.set_ylabel('Compartments Strength', fontsize=14)
    ax.set_xticklabels(names, fontsize=14)
    ax.tick_params(labelsize=14)
    sns.despine(trim=True)
    change_width(ax, .5)
    plt.savefig(args.out, dpi=300, bbox_inches='tight')


def plotLRPerChrom(args):
    """
    %(prog)s sample1.bg sample2.bg [Options]

        To plot barplot of two sample pca1 linregression per chromosome.
    """
    
    p = p=argparse.ArgumentParser(prog=plotLRPerChrom.__name__,
                        description=plotLRPerChrom.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bg', nargs=2, help='The pca1 bg file')
    pReq.add_argument('-o', '--out', 
            help='output file name')

    pOpt.add_argument('--plot', action='store_true', default=False, 
            help='plot the all lineregression picture of per chromosome [default: %(default)s')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    args = p.parse_args(args)
    infile1, infile2 = args.bg
    list(map(check_file_exists, args.bg))
    df1 = pd.read_csv(infile1, header=None, sep='\t', index_col=0, 
            names=[ 'chrom', 'start', 'end', 'pca1'])
    df2 = pd.read_csv(infile2, header=None, sep='\t', index_col=0, 
            names=[ 'chrom', 'start', 'end', 'pca1'])
    df1['type'] = df1['pca1'].map(lambda x: 'A' if x>=0 else 'B')
    df2['type'] = df2['pca1'].map(lambda x: 'A' if x>=0 else 'B')
    
    chroms = sorted(set(df1.index), key=lambda x: x[3:])

    def calc_linregress(value1, value2):
        return  scipy.stats.linregress(value1, value2).rvalue
    def calc_perchrom(chrom, df1, df2):
        return calc_linregress(df1.loc[chrom].pca1, df2.loc[chrom].pca1)

    rvalues = [calc_perchrom(chrom, df1, df2) for chrom in chroms]
    r_df = pd.DataFrame({0: dict(zip(chroms, rvalues))})
    r_df.to_csv(args.out.rsplit(".", 1)[0] + ".tsv", sep='\t', header=None)

    rvalues = np.array(rvalues)
    fig, ax = plt.subplots(figsize=(7, 4))
    colors = ['#D33B00' if _y >=0 else '#00426D' for _y in rvalues]
    sns.barplot(list(range(15)), rvalues, palette=colors, ax=ax)
    ax.set_xticklabels(chroms, rotation=45, ha='right')
    ax.set_ylabel('r',fontsize=18)
    plt.tick_params(labelsize=12, width=1.5)

    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    sns.despine(trim=True)
    plt.savefig(args.out, dpi=300, bbox_inches='tight')  
    plt.savefig(args.out.rsplit('.', 1)[0] + '.png', dpi=300, 
            bbox_inches='tight')
    logging.debug('Successful, picture is in `{}`'.format(args.out))

    if args.plot:
        prefix1 = infile1.split("_", 1)[0]
        prefix2 = infile2.split("_", 1)[0]
        outdir = prefix1 + "-" + prefix2 + "out"
        try:
            os.makedirs(outdir)
        except OSError:
            logging.warning('The such file of `{}` is already exists'.format(outdir))
        
        for chrom in chroms:
            outfile = "{}/{}-{}_{}.pdf".format(outdir, 
                            prefix1, prefix2, chrom)
            xlabel = "{} {} (PC1)".format(prefix1, chrom)
            ylabel = "{} {} (PC1)".format(prefix2, chrom)
            v1 = df1.loc[chrom].pca1
            v2 = df2.loc[chrom].pca1
            abc = ABComparisionSpecies
            abc.plotLineRegress(v1, v2, outfile, 
                xlabel=xlabel, ylabel=xlabel)

        logging.debug('Successful, result is in `{}`'.format(outdir))


def get_chromsize(df):
    """
    Get total chromosize from a bedgraph file with pandas dataframe format
    
    Params:
    --------
    df: pandas.DataFrame 

    Return:
    --------
    out: `int` genome sizes

    Examples:
    ---------
    >>> get_chromsize(df)
     304930202
    """
    chroms = set(df.index.to_list())
    size = 0
    for chrom in chroms:
        size += df.loc[chrom].end.max()
    return size

def statData(bg):
    """
    Stat A/B compartments bedgraph file.

    Params:
    --------
    bg: `str` bedgraph file of compartments pca1
    
    Returns:
    --------
    out: `pandas.DataFrame` of compartments informations

    Examples:
    --------
    >>> df = statData('sample.bg')
    >>> df.head()
    chrom genome_size a b a_rate b_rate
    chr1 3000000000 1600000000 1400000000 0.53 0.47
    ...
    """
    df = pd.read_csv(bg, header=None, sep='\t', index_col=0, 
        names=[ 'chrom', 'start', 'end', 'pca1'])
    df['type'] = df['pca1'].map(lambda x: 'A' if x>=0 else 'B')
    chroms = sorted(set(df.index.to_list()), key=lambda x: x[3:])

    header = ('chrom', 'genome_size', 'a', 'b', 'a_rate', 'b_rate')
    db = {i:[] for i in header}
    
    for chrom in chroms:
        tmp_df = df.loc[chrom]
        a_df = tmp_df[tmp_df.type == 'A']
        b_df = tmp_df[tmp_df.type == 'B']
        alen = sum(a_df.end - a_df.start)
        blen = sum(b_df.end - b_df.start)
        for i, j in zip(header, [chrom, alen+blen, alen, blen, 
                alen*1.0/(alen+blen), blen*1.0/(alen+blen)]):
            db[i].append(j)
        res_df = pd.DataFrame(db)
        
    return res_df

def statAB_old(args):
    """
    %(prog)s eigen1.bg

        Stat A/B compartments per chromosome and total.

    """
    p = p=argparse.ArgumentParser(prog=statAB_old.__name__,
                        description=statAB_old.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bg', help='input A/B compartments bedgraph file')
    pOpt.add_argument('-o', '--out', type=argparse.FileType('w'), 
            default=sys.stdout, help='output file [default: stdin]')
    pOpt.add_argument('--plot', action='store_true', default=False, 
            help='plot the result [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    bg = BedGraph(args.bg)
    a_db = {chrom:0 for chrom in bg.chromList}
    b_db = {chrom:0 for chrom in bg.chromList}

    for chrom in bg.chromList:
        for item in bg.bedGraphDict[chrom]:
            size = item.end - item.begin
            if item.data > 0:
                a_db[chrom] += size
            if item.data < 0:
                b_db[chrom] += size
        chromsize = bg.chromSizes[chrom]
        a_size = a_db[chrom]
        a_perc = a_size / chromsize
        b_size = b_db[chrom]
        b_perc = b_size / chromsize
        print("\t".join(map(str, [chrom, a_size, a_perc, 
                b_size, b_perc, chromsize])), file=args.out)
        
    
    logging.debug('Successful ... output file is in `{}`'.format(args.out.name))

def statAB(args):
    """
    %(prog)s eigen1.bg [eigen1.bg] [Options]

        Stat A/B compartments per chromosome and total.
        Also can visualization the use through `--plot`.

    """
    p = p=argparse.ArgumentParser(prog=statAB.__name__,
                        description=statAB.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bg', nargs="+", help='input A/B compartments bedgraph file')
    pOpt.add_argument('-o', '--out', default='statAB.pdf', 
            help='output file [default: %(default)s]')
    pOpt.add_argument('--plot_total', action='store_true', default=False,
            help='plot the total result [default: %(default)s]')
    pOpt.add_argument('--labels', nargs="*", 
            help='labels of total barplot ylabel [default: bg_prefix]')
    pOpt.add_argument('--plot', action='store_true', default=False, 
            help='plot the result [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    colors = ("#bb4853", "#209093")
    
    def plot(ax, idx, df):
        data = np.array([df.a_rate, df.b_rate]).T
        data_cum = data.cumsum(axis=1)
        
        labels = np.array(df.chrom)[::-1]
        for i, (colname, color) in enumerate(zip(_types, colors)):
            widths = data[:, i][::-1]
            starts = data_cum[:, i][::-1] - widths
            ax.barh(labels, widths, left=starts, height=0.5, 
                        label=colname, color=color)
        if idx == 0:
            ax.tick_params(left='off', labelsize=12, width=0)
        else:
            ax.yaxis.set_visible(False)
        
        ax.xaxis.set_visible(False)
        return ax
    
    def plot_total(ax, idx, df, label):
        genome_size = df.genome_size.sum()
        a_len = df.a.sum()
        b_len = df.b.sum()
        a_rate = a_len / genome_size
        b_rate = b_len / genome_size
        data = np.array([a_rate, b_rate])
        data_cum = data.cumsum()

        for i, (colname, color) in enumerate(zip(_types, colors)):
            widths = data[i]
            starts = data_cum[i] - widths
            ax.barh(label, widths, left=starts, height=0.5, label=colname)
            ax.tick_params(left='off', labelsize=16, width=0)
            xcenters = starts + widths/2
            #for y, (x, c) in enumerate(zip(xcenters, widths)):
            x = xcenters
            y = 0
            c = widths
            ax.text(x, y, "{:.2%}".format(c), ha='center', va='center', fontsize=16)
            ax.xaxis.set_visible(False)
            ax.tick_params(labelsize=12, width=0)
            #ax.yaxis.set_visible(False)

    _types = ('A', 'B')
    df_list = []
    for bg in args.bg:
        df = statData(bg)
        df.to_csv(bg.rsplit('.', 1)[0] + '_stat.tsv', sep='\t', header=None, index=None)
        df_list.append(df)
    
    if args.plot:
        fig, axes = plt.subplots(1, len(df_list), figsize=(3, 3*len(df_list)))
        axes = np.array(axes)
        print(list(map(type, axes)))
        for i, (df, ax) in enumerate(zip(df_list, axes)):
            plot(ax, i, df)
        
        axes[1].legend(ncol=len(_types), bbox_to_anchor=(-1.3, -0.1), 
                loc='lower left')
        sns.despine(bottom=True,left=True)
        plt.savefig(args.out, dpi=300, bbox_inches='tight')
        plt.savefig(args.out.rsplit('.', 1)[0] + '.png', 
                dpi=300, bbox_inches='tight')
    
    if args.plot_total:
        if not args.labels:
            labels = list(map(lambda x: x.split("_")[0], args.bg))
        else:
            labels = args.labels
        
        fig, axes = plt.subplots(len(df_list), 1, figsize=(10, len(df_list)))
        axes = np.array(axes)
        for i, (df, ax) in enumerate(zip(df_list, axes)):
            plot_total(ax, i, df, labels[i])
        plt.legend(bbox_to_anchor=(1, 0.5))
        sns.despine(bottom=True, left=True)
        plt.savefig(args.out, dpi=300, bbox_inches='tight')
        plt.savefig(args.out.rsplit('.', 1)[0] + '.png', 
                dpi=300, bbox_inches='tight')
        

def buscoGeneDist(args):
    """
    %(prog)s sample_eigen1.bg all_gene.bed busco_gene.bed [Options]

        To analysis busco gene distribution between A and B.
    """

    p = p=argparse.ArgumentParser(prog=buscoGeneDist.__name__,
                        description=buscoGeneDist.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bg', help='bedGraph file of compartment pca.')
    pReq.add_argument('genebed', help='bed file of all genes')
    pReq.add_argument('buscobed', help='bed file of busco genes.')
    pOpt.add_argument('-l', '--label', type=str, default='BUSCO',
            help='the label of legend [default: %(default)s')
    pOpt.add_argument('--include', action='store_true', default=False, 
            help='include the buscogene to random selece [default: %(default)s')
    pOpt.add_argument('-f', '--fraction', default=0.5, type=float, 
            help='Minimum overlap requred as a fraction of B [default: %(default)s]')
    pOpt.add_argument('-o', '--out', 
            help='output picture [default: sample_dist.pdf]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    prefix = op.basename(args.bg).rsplit("_")[0]
    out = prefix + "_busco_gene_dist.pdf" if not args.out else args.out
    
    ## annotate compartment
    tsv = args.bg.replace('.bg', '.annotated.tsv')
    annotateType([args.bg, '-o', tsv])
    
    ## intersect gene bed and busco gene bed to tsv
    cmd_formatter = 'bedtools intersect -a {} -b {} -wo -F {} > {}'
    cmd_gene = cmd_formatter.format(tsv, args.genebed, args.fraction, 
            prefix + '_gene_compartment.tsv')
    cmd_busco = cmd_formatter.format(tsv, args.buscobed, args.fraction,
            prefix + '_{}_compartment.tsv'.format("_".join(args.label.split())))
    os.system(cmd_gene)
    os.system(cmd_busco)

    def importData(infile):
        df = pd.read_csv(infile, sep='\t', header=None, index_col=0, 
                    names=('chrom', 'start', 'end', 'pca1', 'type', 
                            'gene_chrom', 'gene_start', 'gene_end',
                            'gene', 'length'))
        return df
    
    def random_gene(df, nums):
        return  df.sample(nums)

    def iter_random(df, nums, _type='A', times=1000):
        res = []
        for i in range(times):
            random_df = random_gene(df, nums)
            res.append(len(random_df.loc[random_df['type'] == _type]))
        return res
    ## import data
    gene_compartment = importData(prefix + '_gene_compartment.tsv')
    busco_compartment = importData(prefix + 
                        '_{}_compartment.tsv'.format("_".join(args.label.split())))

    ## remove busco gene from all genes
    if not args.include:
        gene_compartment = gene_compartment[~gene_compartment['gene'].isin(
                                                    busco_compartment['gene'])]
    ## iter random select gene 
    randomA = iter_random(gene_compartment, len(busco_compartment), _type= 'A')
    randomB = iter_random(gene_compartment, len(busco_compartment),_type = 'B')
    buscoA = len(busco_compartment[busco_compartment['type'] == 'A'])
    buscoB = len(busco_compartment[busco_compartment['type'] == 'B'])
    sns.set()
    ax = sns.distplot(randomA, color='#a83836', label='Random genes in A')
    sns.distplot(randomB, ax=ax, color='#275e8c', label='Random genes in B')
    ax.axvline(buscoA, 1, 0,linestyle='-.', color='#a83836',
            label='{} genes in A'.format(args.label))
    ax.axvline(buscoB, 1, 0, linestyle='-.', color='#275e8c', 
            label='{} genes in B'.format(args.label))
    ax.tick_params(labelsize=12)
    ax.set_xlabel('Gene Numbers', fontsize=14)
    ax.set_ylabel('Frequency', fontsize=14)
    ax.legend()
    plt.savefig(out, dpi=300, bbox_inches='tight')
    logging.debug('Done, picture is in `{}`'.format(out))


def getSwitchLink(args):
    """
    %(prog)s <sample1.bg> <sample2.bg> [Options]

        To get two sample AB swith links for plotting circos,
            also can return 
        **sample.bg
            chrom start end value
        
    """
    p = p=argparse.ArgumentParser(prog=getSwitchLink.__name__,
                        description=getSwitchLink.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bg1', help='bedGraph of pca1 values')
    pReq.add_argument('bg2', help='bedGraph of pca1 values')
    pReq.add_argument('-o', '--out', type=str, required=True,
            help='output file of all annotated table')
    pOpt.add_argument('-1', '--suffix1', type=str, default='',
            help='suffix of bg1 chromosome id [default: %(default)s]')
    pOpt.add_argument('-2', '--suffix2', type=str, default='',
            help='suffix of bg2 chromosome id [default: %(default)s]')
    pOpt.add_argument('--shrink', default=0, type=int,
            help='shrink the link ranges [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    df1 = pd.read_csv(args.bg1, 
                        sep='\t', 
                        header=None, 
                        names=('chrom1', 'start1', 
                                'end1', 'value1'))
    df2 = pd.read_csv(args.bg2, 
                        sep='\t', 
                        header=None, 
                        names=('chrom2', 'start2',
                                'end2', 'value2'))
    df_concatted = pd.concat([df1, df2], axis=1)
    df_concatted['type'] = df_concatted.apply(switch_type_in_df, 
                                                                axis=1)
    df_concatted.to_csv(args.out, sep='\t', header=None, index=None)
    logging.debug('Successfule output all annotated file to `{}`'.format(args.out))
    for _type in set(df_concatted.type):
        tmp_df = df_concatted[df_concatted.type == _type].copy()
        
        tmp_df['chrom1'] = tmp_df['chrom1'].map(lambda x: x + args.suffix1)
        tmp_df['chrom2'] = tmp_df['chrom2'].map(lambda x: x + args.suffix2)
        link_df = tmp_df[['chrom1', 'start1', 'end1',
                            'chrom2', 'start2', 'end2']].copy()
        if args.out is not sys.stdout:
            link_out = args.out.rsplit('.', 1)[0] + '_{}.links'.format(_type)
        else:
            link_out = args.out
        link_df['start1'] = link_df['start1'] + args.shrink
        link_df['end1'] = link_df['end1'] - args.shrink
        link_df['start2'] = link_df['start2'] + args.shrink
        link_df['end2'] = link_df['end2'] + args.shrink
        link_df.to_csv(link_out, sep='\t', header=None, index=None)
        logging.debug('Successful output `{}` to `{}`'.format(_type, link_out))



def import_bg5(infile):
    """
    To import ab data of bedgraph five columns.
        **chrom start end pca1 score
    """
    df = pd.read_csv(infile, 
                     sep='\t', 
                     header=None, 
                     index_col=None, 
                     names=('chrom', 'start', 
                            'end', 'pca1', 
                            'score'))
    df['type'] = df.pca1.map(lambda x: 'A' if x >= 0 else 'B')
    return df


def plotBoxPerChrom(args):
    """
    %(prog)s <infile.bg> [Options]

    To plot boxplot of some data of per A/B compartment per chromosomes.
    bg file: chrom start end pca1 score


    ** if using xlabel or ylabel options, 
        you can also use latex to set special format (italic, up ...)
        italic: `"$\mathic{{trans}}$"`
        scientific count: `"$10 \times 6$"`
    """

    p = p=argparse.ArgumentParser(prog=plotBoxPerChrom.__name__,
                        description=plotBoxPerChrom.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bg', help='bedgraph file of ab_data (five columns)')
    pOpt.add_argument('--chrom', nargs="+", 
            help='chromosomes order of picture [default: all]')
    pOpt.add_argument('--scale', type=int, default=1,
            help='scale of yticks [default: %(default)s]')
    pOpt.add_argument('--xlabel', default='Chromosomes',
            help='xlabel of picture [default: %(default)s]')
    pOpt.add_argument('--ylabel', default='', 
            help='ylabel of picture [default: %(default)s]')
    pOpt.add_argument('-o', '--out', default='',
            help='output file of picture [default: inputfile.pdf]' )
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    from matplotlib.cbook import boxplot_stats
    from ..apps.utilities import wi_test

    whiskerprops = dict(linestyle='--')
    ab_df = import_bg5(args.bg)
    ab_df['score'] = ab_df['score'] / float(args.scale)
    if not args.chrom:
        chrom_list = sorted(set(ab_df.chrom))
    
    fig, ax = plt.subplots(figsize=(10, 5))
    sns.boxplot(x='chrom', y='score', 
            data=ab_df, hue='type', 
            ax=ax, hue_order=['A', 'B'],
            palette=['#a83836', '#265e8a'],
            showfliers=False,
           saturation=1, whiskerprops=whiskerprops)
    ## beautiful
    ax.tick_params(width=1.5, labelsize=13)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.set_xticklabels(chrom_list, rotation=45)
    ax.set_xlabel(args.xlabel, fontsize=16)
    ax.set_ylabel(args.ylabel, fontsize=16)
    sns.despine(trim=True)
    plt.legend(fontsize=13)

    for i, chrom in enumerate(chrom_list):
        a_data = ab_df[(ab_df.chrom == chrom) & (ab_df.type == 'A')].score 
        b_data = ab_df[(ab_df.chrom == chrom) & (ab_df.type == 'B')].score 
        
        pvalue = wi_test(a_data, b_data)
        if 0.01 <= pvalue < 0.05:
            star = '*'
        elif pvalue < 0.01:
            star = "**"
        else:
            star = ''
        
        a_whishi = boxplot_stats(a_data)[0]['whishi']
        b_whishi = boxplot_stats(b_data)[0]['whishi']
        h = max(a_whishi, b_whishi)
        plt.text(i, h , star, fontsize=12, ha='center')
    if not  args.out:
        out = op.basename(args.bg).rsplit('.', 1)[0] + '_perchrom.pdf'
    
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.savefig(out.rsplit('.', 1)[0] + '.png', 
            dpi=300, bbox_inches='tight')


def plotBoxMultiSamples(args):
    """
    %(prog)s <infile.bg> [Options]run

    To plot boxplot of some data of per A/B compartment per samples 
            in one chromosome.
    bg file: chrom start end pca1 score


    ** if using xlabel or ylabel options, 
        you can also use latex to set special format (italic, up ...)
        italic: `"$\mathic{{trans}}$"`
        scientific count: `"$10 \times 6$"`
    """

    p = p=argparse.ArgumentParser(prog=plotBoxMultiSamples.__name__,
                        description=plotBoxMultiSamples.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bg', nargs="+", 
            help='bedgraph file of ab_data (five columns)')
    pReq.add_argument('--chrom', nargs=1, 
            help='chromosomes of picture', required=True)
    pOpt.add_argument('--scale', type=int, default=1,
            help='scale of yticks [default: %(default)s]')
    pOpt.add_argument('--xticklabels', nargs='*', 
            help='xticklabels of picture, length must equal to `input file`'
                    '[default: bg file prefix')
    pOpt.add_argument('--xlabel', default='Chromosomes',
            help='xlabel of picture [default: %(default)s]')
    pOpt.add_argument('--ylabel', default='', 
            help='ylabel of picture [default: %(default)s]')
    pOpt.add_argument('--plotWidth', type=float, default=4, 
            help='Plot width in cm. [default: %(default)s]')
    pOpt.add_argument('--plotHeight', type=float, default=5,
            help='Plot height in cm. [default: %(default)s]')
    pOpt.add_argument('-o', '--out', default='',
            help='output file of picture [default: inputfile.pdf]' )
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    from matplotlib.cbook import boxplot_stats
    from ..apps.utilities import wi_test

    whiskerprops = dict(linestyle='--')


    if not args.xticklabels:
        prefixes = list(map(lambda x: x.split('_', 1)[0], args.bg))
    else:
        prefixes = args.xticklabels
    
    df_list = []
    for i, bg in enumerate(args.bg):
        
        df = import_bg5(bg)
        tmp_df = df[df['chrom'] == args.chrom[0]]
        tmp_df['chrom'] = tmp_df.chrom.map(lambda x: prefixes[i])
        tmp_df['score'] = tmp_df['score'] / float(args.scale)
        df_list.append(tmp_df)

    ab_df = pd.concat(df_list, axis=0)
    fig, ax = plt.subplots(figsize=(args.plotWidth, args.plotHeight))
    sns.boxplot(x='chrom', y='score', 
            data=ab_df, hue='type', 
            ax=ax, hue_order=['A', 'B'],
            palette=['#a83836', '#265e8a'],
            showfliers=False,
           saturation=1, whiskerprops=whiskerprops)
    ## beautiful
    ax.tick_params(width=1.5, labelsize=13)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.set_xticklabels(prefixes, rotation=45)
    ax.set_xlabel(args.xlabel, fontsize=16)
    ax.set_ylabel(args.ylabel, fontsize=16)
    sns.despine(trim=True)
    plt.legend(fontsize=13)

    for i, prefix in enumerate(prefixes):
        a_data = ab_df[(ab_df.chrom == prefix) & (ab_df.type == 'A')].score 
        b_data = ab_df[(ab_df.chrom == prefix) & (ab_df.type == 'B')].score 
        
        pvalue = wi_test(a_data, b_data)
        if 0.01 <= pvalue < 0.05:
            star = '*'
        elif pvalue < 0.01:
            star = "**"
        else:
            star = ''
        
        a_whishi = boxplot_stats(a_data)[0]['whishi']
        b_whishi = boxplot_stats(b_data)[0]['whishi']
        h = max(a_whishi, b_whishi)
        plt.text(i, h , star, fontsize=12, ha='center')
    if not  args.out:
        out = '{}_ab_boxplot_{}.pdf'.format( "-".join(map(op.basename, 
                                        prefixes)), args.chrom[0])
    
    plt.savefig(out, dpi=300, bbox_inches='tight')
    plt.savefig(out.rsplit('.', 1)[0] + '.png', 
            dpi=300, bbox_inches='tight')


def quickPlot(args):
    """
    %(prog)s <eigen1.bg> [gene.density.bg] [Options]
        to quick plot picture to visualized A/B compartments

    """

    p = argparse.ArgumentParser(prog=quickPlot.__name__,
                            description=quickPlot.__doc__,
                            conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('eigen1', help='bg file of eigen1')
    pReq.add_argument('chromsizes', help='chromosome sizes file')

    pOpt.add_argument('-g', '--gene', 
            help='gene density bg file')
    pOpt.add_argument('--TE', help='bg file of TE density')
    pOpt.add_argument('-o', '--outdir', default='quickPlot',
            help='outdir [default: %(default)s]')
    pOpt.add_argument('--pdf', default=False, action='store_true',
            help='if output pdf format [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    outdir=op.abspath(args.outdir)
    import configparser
    cf = configparser.ConfigParser()
    
    cf.add_section("ab")
    cf.set('ab', 'file', args.eigen1)
    cf.set('ab', 'color', '#BB4853')
    cf.set('ab', 'negative color', '#209093')
    cf.set('ab', 'title', 'A/B Compartments')
    cf.set('ab', 'fontsize', '14')
    cf.set('ab', 'height', '3')
    cf.set('ab', 'max_value', '0.15')
    cf.set('ab', 'min_value', '-0.15')

    if args.gene:
        cf.add_section('spacer')
        cf.add_section('gene')
        cf.set('gene', 'file', args.gene)
        cf.set('gene', 'color', '#D33B00')
        cf.set('gene', 'title', 'Gene Density')
        cf.set('gene', 'height', '3')
        cf.set('gene', 'fontsize', '14')
        cf.set('gene', 'max_value', '30')
    
    if args.TE:
        cf.add_section('spacerTE')
        cf.add_section('te')
        cf.set('te', 'file', args.TE)
        cf.set('te', 'color', '#00426D')
        cf.set('te', 'title', 'TE Density')
        cf.set('te', 'height', '3')
        cf.set('te', 'fontsize', '14')
        #cf.set('gene', 'max_value', '1500')
    
    
    cf.add_section('x-axis')
    if not op.exists(outdir):
        os.makedirs(outdir)
    with open('{}/quickplot.ini'.format(outdir), 'w+') as f:
        cf.write(f)
    
    os.system("sed -i 's/spacerTE/spacer/g' {}/quickplot.ini".format(outdir))

    chromsizes = dict(i.strip().split() for 
                    i in open(args.chromsizes) 
                    if i.strip())
    plot_cmd_formatter = 'pyGenomeTracks --tracks {}/quickplot.ini '.format(outdir)
    plot_cmd_formatter += '-o {1}/{0}.{2} --region {0}'

   

    fmt = 'pdf' if args.pdf else 'png'
    for chrom in chromsizes:
        start = 1
        end = chromsizes[chrom]
        region = '{}:{}-{}'.format(chrom, start, end)
        print(plot_cmd_formatter.format(region, args.outdir, fmt))

    
if __name__ == "__main__":
    main()
