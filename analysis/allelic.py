#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Allelic genes analysis.
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from TDGP.apps.base import debug, check_file_exists, listify
from TDGP.apps.base import ActionDispatcher

debug()

def main():

    actions = (
            ("test", "test"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())




## out command ##

def allelicTableToList(args):
    """
    %(prog)s allelic.table [Options]

    """
    p = p=argparse.ArgumentParser(prog=allelicTableToList.__name__,
                        description=allelicTableToList.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('table', 
            help='the table of allelic genes'
            '#num geneA geneB geneC geneD')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    
    