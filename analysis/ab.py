#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
A/B compartments analysis libraries
"""

import logging
import sys
import os
import os.path as op
import pandas as pd

from collections import defaultdict
from TDGP.apps.base import debug

debug()
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)




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
    judge the switch type of ab compartments in two species: `AA`, `BB`, `A2B`, `B2A`
    >>> switch_type(-0.1, 0.2)
    'B2A'
    >>> switch_type(0.1, 0.1)
    'AA'
    """
    if v1 > 0:
        if v2 > 0:
            return 'A Stable'
        else:
            return 'A2B'
    else:
        if v2 > 0:
            return 'B2A'
        else:
            return "B Stable"


def two_species_conserved_compare(tgy_df, jgy_df):
    """
    
    """
    db = {"A Stable": 0, "B Stable": 0, "A2B": 0, "B2A": 0}
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
    mpl.rcParams['font.size'] = 16.0
    plt.figure(figsize=(10,10))
    _, texts, autotexts = plt.pie((conserved, non_conserved), labels=labels, colors=colors, autopct='%1.2f%%')
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontsize(20)


def plot_two_species_ab_pie(tgy_df, jgy_df):
    """
    Pie plot of two species A/B switch type.
    """
    db = two_species_conserved_compare(tgy_df, jgy_df)
    func = lambda x: "{} ({})".format(x[0], chrom_size_convert(x[1]))
    labels = list(map(func, db.items()))
    colors = ['#265e8a', '#032F49','#BB4853', '#209093','#a83836',]
    mpl.rcParams['font.size'] = 16.0
    plt.figure(figsize=(10,10))
    _, texts, autotexts = plt.pie(db.values(), labels=labels, colors=colors, autopct='%1.2f%%')
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontsize(20)




if __name__ == "__main__":
    main()
