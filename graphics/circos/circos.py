#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
tools for circos plot.
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from TDGP.apps.base import ActionDispatcher
from TDGP.apps.base import check_file_exists, debug
from TDGP.apps.base import listify



def main():

    actions = (
            ("test", "test"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())  


## out command
def karyotype(args):
    """
    %(prog)s <chrom.sizes> [Options]

        generate karyotype file for circos
    """

    p = p=argparse.ArgumentParser(prog=karyotype.__name__,
                        description=karyotype.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('sizes', help='chromosome sizes file `chrom size`')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    pass 