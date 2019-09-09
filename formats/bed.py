#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Time: 2019/7/22 8:38

"""
bed format
"""
from __future__ import print_function

import os.path as op
import os
import sys


class BedpeLine(object):
    def __init__(self):
        self.chrom1 = ''
        self.start1 = 0
        self.end2 = 0
        self.chrom2 = ''
        self.start2 = 0
        self.end2 = 0
        self.name = '.'
        self.score = '.'
        self.strand1 = '.'
        self.strand2 = '.'
        self.other = []
