#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Basic support for running library as script
"""

from __future__ import print_function

import logging
import os
import os.path as op
import sys


def debug(level=logging.DEBUG):
    """
    Basic config logging format
    """
    from 3
    formats = "%(asctime)s [%(module)s]"
    formats += " %(message)s"
    logging.basicConfig(level=level, format=formats, datefmt="%H:%M:%S")

debug()

