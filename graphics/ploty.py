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

