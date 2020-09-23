#!/usr/bin/env python
# -*- coding:utf-8 -*-

# My data valization libraries
"""
ploty: my data valization libraries. 
"""
from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns
from scipy.stats import linregress
from matplotlib.lines import Line2D

from TDGP.apps.base import listify

def main():

    actions = (
            ("test", "test"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

class Venn(object):
    COLOR2 = ()
        
    def __init__(self):
        pass

    @staticmethod
    def plotVenn2(F):
        pass

# https://stackoverflow.com/questions/34888058/
# changing-width-of-bars-in-bar-chart-created-using-seaborn-factorplot
def change_width(ax, new_value) :
    """
    Change barplot width in seaborn.barplot

    Examples:
    --------
    >>> ax = sns.barplot(x=[0, 1], y=[0.39380692778678045, 0.37079504302939925], 
    yerr=[0.009796196230459669, 0.009569196380718775])
    >>> change_width(ax, .5)

    """
    for patch in ax.patches :
        current_width = patch.get_width()
        diff = current_width - new_value

        # we change the bar width
        patch.set_width(new_value)

        # we recenter the bar
        patch.set_x(patch.get_x() + diff * .5)


def trim_axes(axes, N):
    """
    little helper to message the axs list to have correct length...

    Params:
    --------
    axes: `list`
            list of axes
    N: `int`
            number of axes to return
    
    Returns:
    --------
    axes: `list`
            list of trimed axes
    
    Examples:
    --------
    >>> fig, axes = plt.subplots(5, 2)
    >>> axes = trim_axes(axes, 7)
    array([<matplotlib.axes._subplots.AxesSubplot object at 0x7f5f48365198>,
       <matplotlib.axes._subplots.AxesSubplot object at 0x7f5f49275438>,
       <matplotlib.axes._subplots.AxesSubplot object at 0x7f5f4712e320>, ...]
    """

    axes = axes.flat
    for ax in axes[N:]:
        ax.remove()
    
    return axes[:N]


def savefig(figname, dpi=300, bbox_inches=None,
            formats=['pdf', 'png'], cleanup=True):
    """
    function for savefig, can save multi format

    Params:
    -------
    figname: `str`
            figname
    dpi: `int`
            dpi for figure [default: 300]
    formats: `list` or `str`
            picture formats [default: ["pdf", "png"]]
    bbox_inches: `str`
            bbox_inches params for plt.savefig [defult: None]
    cleanup: `bool`
            if cleanup rcParms after savefig [default: True]
    
    Returns:
    --------
    out: output figure

    Examples:
    --------
    >>> savefig('out')
    
    """
    formats = listify(formats)

    try:
        if not op.splitext(figname)[-1]:
            raise
        else:
            fmt = op.splitext(figname)[-1].lower()
    except:
        fmt = "pdf"
    if fmt not in formats:
        formats.append(fmt)
    
    for fmt in formats:
        figprefix = op.splitext(figname)[0]
        outname = "{}.{}".format(figprefix, fmt)
        plt.savefig(outname, dpi=dpi, format=fmt, 
                bbox_inches=bbox_inches)
        msg = "Figure saved to `{}`".format(outname)
        logging.debug(msg)
    
    ## reset rcParams after savefig
    if cleanup:
        plt.rcdefaults()


    

def plotLineRegress(ax, xdata, ydata, 
                    xlabel='species1', 
                    ylabel='species2', 
                    scatter=True):
    """
    Function of lineregress ploting

    Params
    --------
    ax: `ax`
    xdata, ydata: `array-like`
    xlabel, ylabel: `str` label of axis
    scatter: `boolean` if to plot scatter

    Returns:
    --------
    out: `ax`

    Examples:
    --------
    >>> plotLineRegress(ax, xdata, ydata)
    """

    assert len(xdata) == len(ydata), "`xdata` and ``ydata length must euqual"
    scatter_params = dict(color='#209093', s=2)
    line_params = dict(color='#032F49', lw=2)
    slope, intercept, rvalue, pvalue, stderr = linregress(xdata, ydata)
    #r2 = rvalue ** 2 if rvalue > 0 else -1 * rvalue ** 2
    label = [r"r = {:.2f}  $\mathit{{P}}$ = {:.2e}".format(rvalue, pvalue)]
    sns.regplot(xdata, ydata, ax=ax, truncate=True, scatter=scatter,
            scatter_kws=scatter_params, line_kws=line_params)
    legend_elements = [Line2D([0], [0], **line_params)]
    #ax.set_title('The regression of PC1')
    ax.set_xlabel("{}".format(xlabel), fontsize=13)
    ax.set_ylabel("{}".format(ylabel), fontsize=13)
    ax.legend(legend_elements, label, loc='best')
            #bbox_to_anchor=(1, 0.5))
    #ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    #ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    return ax


if __name__ == "__main__":
    main()
