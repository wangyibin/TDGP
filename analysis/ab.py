#!/usr/bin/env python
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

debug()
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def main():

    actions = (
            ("test", "test"),
            ("plotLineRegress", "Plot two species pca1 lineregress"),
            ('plotMultiLineRegress', "plot two species pca1 lineregress per switch type"),
            ('plotEnrichment', 'plot compartment strength enrichment'),
            ('plotStrength', 'plot compartment strength in multi samples'),
            ('getSyntenyGenePca', "get the synteny gene pairs pca value"),
            ('annotateSwitchType', "annotate swithch type for synteny gene pairs"),
            
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
            return np.log(i[0,0] * i[-1,-1] / i[0,-1]**2)
            t1 = i[:3,:3].sum() + i[:2, :2].sum()
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

            
            gc = np.array(eig[genome.label2idx[chrom]]) if eig == 'GC' \
                else np.array(eig[chrom])
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
        map(check_file_exists, (bg1, bg2, bed1, bed2))
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
        map(os.system, (bedtools_cmd1, bedtools_cmd2, cut_cmd1, cut_cmd2, paste_cmd))
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
        label = [r"R$^2$ = {:.2f}  $\mathit{{P}}$ = {:.2e}".format(rvalue, pvalue)]
        fig, ax = plt.subplots(figsize=(5, 5))
        regplot(a, b, ax=ax, truncate=True, 
                scatter_kws=scatter_params, line_kws=line_params)
        legend_elements = [Line2D([0], [0], **line_params)]
        #ax.set_title('The regression of PC1')
        ax.set_xlabel("{}".format(xlabel))
        ax.set_ylabel("{}".format(ylabel))
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
            text = '{} ($R^2$={:.2f})'.format(stype, rvalue)
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
    df = pd.read_csv(path, sep='\t', header=None, names=[ 'chrom', 'start', 'end', 'pca1'], index_col=0)
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

## outside command ##
def plotLineRegress(args):
    """
    %(prog)s a.pca1 b.pca1 [Options]
        To plot lineregress of two species' pca1 eigenvector
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
    pOpt.add_argument('-h', '--help', action='help', 
            help='show help message and exit.')

    args = p.parse_args(args)

    ABC = ABComparisionSpecies()
    check_file_exists(args.a)
    check_file_exists(args.b)
    a = [float(i.strip()) for i in open(args.a)]
    b = [float(i.strip()) for i in open(args.b)]
    ABC.plotLineRegress(a, b, args.out, xlabel=args.xlabel, 
            ylabel=args.ylabel)


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


def annotateSwitchType(args):
    """
    %(prog)s species1-species2.synteny.eigen1.bg [Options]
        
        To annotate the A/B switch type to a synteny gene pairs eigen file.
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
    pOpt.add_argument('--same', action='store_true', default=True,
            help='if same species [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')

    args = p.parse_args(args)

    ABC = ABComparisionSpecies if not args.same else ABComparision
    ABC.annotateSwitchType(args.infile, args.out, args.stype)


def plotMultiLineRegress(args):
    """
    %(prog)s in.annotated.tsv [Options]
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
    %(prog)s infile1 infile2 [Options]

        Plot pie of compartments switch.
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
    pOpt.add_argument('-t', '--thread', type=int, default=24,
            help='the thread of programe [default: %(default)s]')
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
            thread=args.thread, iterCorrect=False)
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
    pOpt.add_argument('-t', '--thread', type=int, default=24,
            help='the thread of programe [default: %(default)s]')
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
    values = []
    errors = []
    for i, matrix in enumerate(data_list):
        strengthes, permutted = Compartment().getStrengthError(matrix, 
            eig[i], genome_list[i], args.window, correct=True, 
            thread=args.thread, iterCorrect=True)
        value, errorCompartment().calStrength(strengthes, permutted)
        values.append(value)
        errors.append(error)
    index = np.arrange(len(values))
    plt.bar(index, values, 0.8, yerr=errors, error_kw=dict(ecolor='0.3'))
        #plt.colorbar()
    
    plt.savefig(args.out, dpi=300, bbox_inches='tight')


  
    

    

    


if __name__ == "__main__":
    main()
