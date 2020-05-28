#!/usr/bin/env python
# -*- coding:utf-8 -*-

# My data valization libraries
"""
ploty: my data valization libraries. 
"""

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

import os
import os.path as op
import sys

import pandas as pd
import seaborn as sns
from scipy.stats import linregress
from matplotlib.lines import Line2D

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
    ax.set_xlabel("{}".format(xlabel), fontsize=12)
    ax.set_ylabel("{}".format(ylabel), fontsize=12)
    ax.legend(legend_elements, label, loc='best')
            #bbox_to_anchor=(1, 0.5))
    
    return ax


if __name__ == "__main__":
    main()
