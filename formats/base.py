#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
the base libraries of formats module
"""


from __future__ import print_function

import logging
import os
import os.path as op
import sys

from TDGP.apps.base import debug


debug()


class BaseFile(object):
    def __init__(self, filename):
        if filename:
            logging.debug("Loading file `{}`".format(filename))


class LineFile(BaseFile, list):
    def __init__(self,filename, comment=None, load=False):
        super(LineFile, self).__init__(filename)
        if load:
            pass


