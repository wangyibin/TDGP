#!/usr/bin/env python
# -*- coding:utf-8 -*-


"""
the libraries of blast out formats
"""

from __future__ import print_function

import logging
import os
import os.path as op
import sys

from TDGP.apps.base import debug
from TDGP.formats.base import BaseFile

debug()

class BlastLine(object):
    def __init__(self, line):
        line_list = line.strip().split()
        self.qseqid= line_list[0]
        self.sseqid = line_list[1]
        self.identity = line_list[2]
        self.length = line_list[3]
        self.mismatch = line_list[4]
        self.gap = line_list[5]
        self.qstart = line_list[6]
        self.qend = line_list[7]
        self.sstart = line_list[8]
        self.send = line_list[9]
        self.evalue = line_list[10]
        self.bitscore = line_list[11]




