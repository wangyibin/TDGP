#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
plot allele genes and annotate variants supported by isoseq
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from argparse import Namespace
from Bio import SeqIO
from collections import OrderedDict
from matplotlib.lines import Line2D
from matplotlib import patches
from itertools import accumulate


class GffRecord:
    def __init__(self, gff):
        self.gff3 = gff
        self.record = self.import_gff3()
        
    def import_gff3(self):
        names = ['chrom', 'source', 'feature', 
                'start', 'end', 'score', 'strand', 
                'frame', 'attribute']
        df = pd.read_csv(self.gff3, sep='\t', header=None, 
                         index_col=None, names=names)
        df.chrom = df.chrom.astype('category')
        df.source = df.source.astype('category')
        df.feature = df.feature.astype('category')
        df.strand = df.strand.astype('category')
        df.frame = df.frame.astype('category')
        df.attribute =  df.attribute.str.rstrip(";").str.split(";")
        df.attribute = df.attribute.apply(lambda x: Namespace(**dict(list(map(lambda y: y.split("="), x)))))

        return df
    
    @property
    def gene(self):
        return self.record[self.record.feature == 'gene']
    @property
    def gene_start(self):
        return self.gene.start.loc[0]
    @property
    def gene_name(self):
        return self.gene.attribute.loc[0].ID
    @property
    def sequence_length(self):
        length = self.gene['end'] - self.gene['start'] + 1
        return length.loc[0]
    
    @property
    def chromosome(self):
        return self.gene.loc[0, 'chrom']
    
    @property
    def exons(self):
        return self.record[self.record.feature == 'exon'].sort_values(by='start')
    @property
    def exon_positions(self):
        return self.exons[['start', 'end']].values
    @property
    def exon_length(self):
        return self.exon_positions[:, 1] - self.exon_positions[:, 0] + 1
    @property
    def cds(self):
        return self.record[self.record.feature == 'CDS'].sort_values(by='start')
    @property
    def cds_positions(self):
        return self.cds[['start', 'end']].values
    @property
    def cds_start(self):
        return self.cds_positions[0, 0]
    @property
    def cds_length(self):
        return self.cds_positions[:, 1] - self.cds_positions[:, 0] + 1

def broken_line(ax, pos1, pos2, text='', text_prop={},
                hshift=0.1, 
                length_fraction=0.2, color='#000000',
               *args, **kwargs):
    line_length = (pos2[1] - pos1[1]) * length_fraction
    first_pos1 = (pos1[0], pos1[0])
    first_pos2 = (pos1[1] , pos1[1] + line_length)
    third_pos1 = (pos2[0] + hshift, pos2[0] + hshift)
    third_pos2 = (pos2[1] - line_length, pos2[1])
    second_pos1 = (first_pos1[1], third_pos1[0])
    second_pos2 = (first_pos2[1], third_pos2[0])
    ax.plot(first_pos1, first_pos2, color=color, *args, **kwargs)
    ax.plot(third_pos1, third_pos2, color=color, *args, **kwargs)
    ax.plot(second_pos1, second_pos2, color=color, *args, **kwargs)
    
    if text:
        ax.text(third_pos1[1], third_pos2[1],# + line_length, 
                text, ha='left', va='top', **text_prop)
    return ax

def addGeneTicks(ax, xpos, ypos, length, tick_length=0.05,
                 scale=500, *args, **kwargs):
    
    ticks_num = length // scale
    ticks_redundance = length % scale
    ticks = np.linspace(0, ticks_num * scale, ticks_num + 1, dtype='int64')
    #ticks = np.hstack([ticks, [ticks_num * scale + ticks_redundance]])
    
    for tick in ticks:
        tick_label = "{:,}".format(tick) 
        xy = (tick / length * 100, ypos)
        xytext = (xy[0] , xy[1]  +  tick_length)
        ax.annotate(tick_label, xy=xy, xytext=xytext, ha='center', 
                    va='bottom', color='#636363',
                   arrowprops=dict(arrowstyle='-', color='#636363'), 
                    *args, **kwargs)
    
    return ax

def getSites(site_file, aln_file, index):
    sites = [i.strip().split() for i in open(site_file) if i.strip()]
    sites = np.array(sites, dtype='int32')
    site_df = pd.DataFrame(columns=sites[:, 0])
    with open(aln, 'r') as fp:
        fa = SeqIO.parse(fp, 'fasta')
        for record in fa:
            seq = []
            for site in sites:
                seq.append(str(record.seq[site[0] - 1: site[0] - 1 + site[1]]).upper())
            site_df.loc[record.name] = seq

    
    site_df = site_df.loc[index]

    return site_df

def getCDSTotalPositions(gffRecord):
    end = list((accumulate(gffRecord.cds_length)))
    start = [0] + end[:-1]
    start = np.array(start) + 1
    end = np.array(end)
    cds_total_positions = pd.DataFrame(np.array(list(zip(start, end))), 
                                       columns=['start', 'end'])
    return cds_total_positions

def convertAlnPos2GenePos(value, cds_total_positions, gffRecord):
    idx, = cds_total_positions[(cds_total_positions.start <= value) & 
                               (cds_total_positions.end >= value)].index.values
    length = value - cds_total_positions.loc[idx].start + 1
    position = length + gffRecord.cds_positions[idx][0] - 1
    
    return ((position - gffRecord.gene_start) / gffRecord.sequence_length) * 100

def plot(gene_name, gene_length, site_df, if_hshift=True):
    strand = "+"
    gene_name = df.gene_name.strip('ab')
    gene_length = df.sequence_length
    index = site_df.index
    annotate_fontsize = 3
    gene_fontsize = 9
    text_prop = dict(fontsize=annotate_fontsize)
    text_pos = 2
    mpl.rcParams['font.style'] = 'arial'
    mpl.rcParams['font.style'] = 'normal'
    # from matplotlib.backends.backend_pgf import FigureCanvasPgf
    # matplotlib.backend_bases.register_backend('pdf', FigureCanvasPgf)
    pgf_with_latex = {
        "text.usetex": False,            # use LaTeX to write all text
        "pgf.rcfonts": False,           # Ignore Matplotlibrc
        "pgf.preamble": [
            r'\usepackage{color}'     # xcolor for colours
        ]
    }
    mpl.rcParams.update(pgf_with_latex)
    
    x, y, dx, dy = 0, 10, 100.3, 0
    if strand == "-":
        x, y, dx, dy = 100, 10, -100.3, 0
    
    fig, ax = plt.subplots()

    ax.broken_barh(exon_data, (9, 2), facecolor="#ffd383", edgecolor='#000000')
    arrow = ax.arrow(x, y, dx, dy, head_width=1.4, 
                    head_length=1, linewidth=0.5,
                    fc="#bcbcbc", ec='#000000'
                    )
    arrow.set_zorder(0)

    if len(site_df.columns) >= 20:
        iso_index_pos = -40
        hshift = -30
    else:
        iso_index_pos = -5
        hshift = 0
    init_hshift = hshift
    ax.text(-5, 10, gene_name, va='center', ha='right', fontsize=gene_fontsize)
    ax.text(iso_index_pos + 5,  text_pos, "\n".join(index),va='top', ha='right', fontsize=annotate_fontsize )
    #ax.hlines(5, 0, 100)
    previous_site = 0
    hshift_num = 0
   
    for i, site in enumerate(site_df.columns):
        seq = "\n".join(site_df[site])
        single_seq_length = len(site_df[site].iloc[1])
        if if_hshift:
            if (previous_site != 0):
                if (site - previous_site) < 0.3:
                    factor = 1.4 + single_seq_length * 0.2
                elif (site - previous_site) < 0.3:
                    factor = 1.2 + single_seq_length * 0.2
                elif (site - previous_site) < 0.5:
                    factor = 1 + single_seq_length * 0.2
                elif (site - previous_site) < 1:
                    factor = 0.6 + single_seq_length * 0.2
                elif (site - previous_site) < 2:
                    factor = 0.4 + single_seq_length * 0.2
                elif (site - previous_site) < 2.5:
                    factor = 0.2 + single_seq_length * 0.2
                elif (site - previous_site) < 3:
                    factor = 0.1 + single_seq_length * 0.2
                elif (site - previous_site) < 4:
                    factor = -0.2 + single_seq_length * 0.2
                else:
                    factor = -1.4 * (site - previous_site)


                hshift += site - previous_site  + factor

        broken_line(ax, (site, 9), (site, text_pos), text=seq, 
                    text_prop=text_prop, hshift=hshift, linewidth=0.3)
        previous_site = site
    addGeneTicks(ax, (0, 100), 11.5, gene_length, 
                scale=500, fontsize=annotate_fontsize + 2, 
                tick_length=1)
    ax.set_ylim(-20, 40)
    ax.set_xlim(-2 + init_hshift, 102 + hshift)
    ax.axis('off')
    plt.savefig('{}.svg'.format(gene_name), bbox_inches='tight', dpi='300')


def plotAlleleGeneModel(args):
    p = argparse.ArgumentParser(prog=plotAlleleGeneModel.__name__,
                        description=plotAlleleGeneModel.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('gff',  help='gff file only one gene record')
    pReq.add_argument('sites', help='variants site information '
                    'two columns `position\tlength`')
    pReq.add_argument('aln', help='iso seq aln file')
    pReq.add_argument('geneName', help='gene name without hap suffix')
    pOpt.add_argument('--hshift', default=False, action='store_true',
            help='if hshift [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    gene_name = args.geneName
    df = GffRecord(args.gff)
    exon_data = np.vstack([(df.exon_positions[:, 0] - df.gene_start + 1) , 
                            df.exon_length]).T
    exon_data = (exon_data / df.sequence_length) * 100

    index = [f'{gene_name}a', 'IsoReads07', 'IsoReads08', 
             'IsoReads09', 'IsoReads10', 'IsoReads01', 
             'IsoReads02', 'IsoReads03', 'IsoReads04', 
             'IsoReads05', 'IsoReads06', f'{gene_name}b']
    site_df = getSites(args.sites, args.aln, index)
    cds_total_positions = getCDSTotalPositions(df)
    site_df.columns = site_df.columns.to_series().apply(convertAlnPos2GenePos,  
                                            args=(cds_total_positions, df))
    
    gene_length = df.sequence_length
    strand = "+"
    
    plot(gene_name, gene_length, site_df, if_hshift=args.hshift)

if __name__ == "__main__":
    plotAlleleGeneModel(sys.argv[1:])