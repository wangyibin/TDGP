#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
generate allele table from synteny results
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 

from util import (pairs, 
                    create_empty_allele_table,
                    import_anchor,
)


def generate_allele_table(args):
    
    p = p=argparse.ArgumentParser(prog=generate_allele_ta
    ble.__name__,
                        description=generate_allele_table.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('', 
            help='')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)