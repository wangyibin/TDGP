#!/usr/bin/env python

import argparse
import logging
import os.path as op
import os
import sys


def ALLHiC_Pipeline(args):
    """
    %(prog)s [Options]
    """
    p = p=argparse.ArgumentParser(prog=ALLHiC_Pipeline.__name__,
                        description=ALLHiC_Pipeline.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('',  help='')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


if __name__ == "__mian__":
    ALLHiC_Pipeline(sys.argv[1:])