#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
RepeatMasker out format
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 

from collections import OrderedDict
from TDGP.apps.base import ActionDispatcher
from TDGP.apps.base import check_file_exists, debug
from TDGP.apps.base import listify



def main():

    actions = (
            ("rename", "rename chromosome name in repeat file"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())

REPEATMASKER_OUT_HEADER = ['score', 'div', 'del', 'ins', 'chrom',
                            'start', 'end', 'left', 'strand',
                            'name', 'class_type', 'repeat_start', 
                            'repeat_end', 'repeat_left', 'repeat_ID' ]

def import_repeat_out(inputfile):
    """
    import repeatMasker out as DataFrame
    """
    df = pd.read_csv(inputfile, sep='\s*', skiprows=5, skipinitialspace=True,
                    index_col=4, usecols=[i for i in range(16)],
                    header=None, engine='python')
    
    return df 


## out command ##

def rename(args):
    """
    %(prog)s <repeat.out> [Options]

    """
    p = argparse.ArgumentParser(prog=rename.__name__,
                        description=rename.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('repeat', 
            help='repeat file from RepeatMasker')
    
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    df = import_repeat_out(args.repeat)
    print(df.head())


if __name__ == "__main__":
    main()

