#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Some properties for plot, such as color, font, colormap.
"""

from __future__ import print_function

import argparse 
import logging
import os
import os.path as op
import sys

import matplotlib as mpl 
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import yaml
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from cycler import cycler
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import rcParams




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


# https://stackoverflow.com/questions/3380726/
# converting-a-rgb-color-tuple-to-a-six-digit-code-in-python
def clamp(x):
    return max(0, min(x, 255))
def rgb2hex(x):
    r, g, b = x
    return "#{0:02x}{1:02x}{2:02x}".format(clamp(r), clamp(g), clamp(b))


class wesColors(object):
    """
    colors from wes
    """
    def __init__(self):
        _colors_path = op.dirname(op.realpath(__file__))
        self._colors_path = op.join(_colors_path, 'colors.yaml')
        self._colors = yaml.load(open(self._colors_path), Loader=Loader)
        
    @property
    def colors_path(self):
        return self._colors_path

    @property
    def colors(self):
        return self._colors

    def _set_colors(self, ccycle):
        rcParams['axes.prop_cycle'] = cycler(color=ccycle)
    
    @classmethod
    def set_palette(self, cname="GrandBudapest1"):
        """
        set matplotlib default colormap to specify colors

        Params:
        --------
        cname: `str` 
                Name of palette shown in wesColors.available()
        """
        try:
            if cname.endswith('_r'):
                colorList = wesColors().colors[cname[:-2]][::-1]
            else:
                colorList = wesColors().colors[cname]
            wesColors()._set_colors(colorList)
        except KeyError:
            raise KeyError("{cname} is not in available. Check wesColors.available()")
    
    @classmethod
    def available(self, plot=False):
        """
        show all availabel color palettes

        """
        colors = wesColors().colors
        if not plot:
            return colors
        else:
            fig, axes = plt.subplots(3, 6, figsize=(12, 12), sharex=True)
            for i, cname in enumerate()
    
    @classmethod
    def plot_palettes(*args):
        if len(args) > 1:
            fig, axes = plt.subplots(1, len(args))
            for i, cname in enumerate(args):
                color = 0

    def get_color_cycle(self):
        sns.palettes.get_color_cycle()