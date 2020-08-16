#!/usr/bin/env python
#! -*- coding -*-


from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from validpairs2links import validpairs2links
from agp2assembly import agp2assembly

def parse_args(name=None, usage=None):

    parser = argparse.ArgumentParser(name,
                        description=usage,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    subparsers = parser.add_subparsers(title='Commands',
                                dest='command')
    subparsers.add_parser('valid2links', 
            help='convert validpairs to 3d-dna links')
    subparsers.add_parser('agp2assembly',
            help='convert agp to 3d-dna assembly file')
    
    
    return parser


def ALLHiC_adjust(args):
    
    
    # valid2links_parser.set_defaults(func=validpairs2links, )
    # pReq = p.add_argument_group('Required arguments')
    # pOpt = p.add_argument_group('Optional arguments')

    # pOpt.add_argument('-h', '--help', action='help',
    #         help='show help message and exit.')
    p = parse_args(name=ALLHiC_adjust.__name__,
                    usage=ALLHiC_adjust.__doc__)
    arg = p.parse_args(args)
    if arg.command == 'valid2links':
        validpairs2links(args[1:])
    elif arg.command == 'agp2assembly':
        agp2assembly(args[1:])
    else:
        p.print_help()

if __name__ == "__main__":
    if len(sys.argv) == 1:
        sys.argv.append("--help")

    ALLHiC_adjust(sys.argv[1:])