#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Some colors or colormaps for plot.
"""

from __future__ import print_function

import os
import os.path as op
import sys

import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

from matplotlib.colors import LinearSegmentedColormap

BlueBlackRed_dict = {
        'red':   ((0.0, 0.0, 0.0),
                   (0.5, 0.0, 0.1),
                   (1.0, 1.0, 1.0)),

        'green': ((0.0, 0.0, 0.0),
                   (1.0, 0.0, 0.0)),

        'blue':  ((0.0, 0.0, 1.0),
                   (0.5, 0.1, 0.0),
                   (1.0, 0.0, 0.0))
        }

blue_red = LinearSegmentedColormap('BlueBlackRed', BlueBlackRed_dict)
plt.register_cmap(cmap=blue_red) 


class MyColor(object):
    posneg = ('#797596', '#0B1D51')
    posneg_r = tuple(reversed(list(posneg)))
    blrd = ("#89023E", "#b4b4b4")
    blrd_r = tuple(reversed(list(blrd)))
