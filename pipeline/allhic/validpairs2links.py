#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
convert hicpro validpairs to 3d-dna links
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from TDGP.analysis.qc import ValidPairs

class LinkLine(object):
    """
    LinkLine object for 3d-dna links line.
    Examples:
    -------
    >>> link = LinkLine()
    >>> link.read1_name = 'Chr1
    >>> link.read1_pos = 123123
    >>> print(link)
    0 Chr1 123123 0 ...
    """
    def __init__(self):
        self.strandA = 0
        self.read1_name = ""
        self.read1_pos = None
        self.order1 = 0
        self.strandB = 0
        self.read2_name = ""
        self.read2_pos = None
        self.order2 = 1
        self.suffix = " 1 - - 1 - - -"
    def __str__(self):
        res = " ".join(map(str, [self.strandA, self.read1_name,
                        self.read1_pos, self.order1,
                        self.strandB, self.read2_name,
                        self.read2_pos, self.order2,
                        self.suffix]))
        return res
    __repr__ = __str__ 
    def vpl2links(self, vpl):
        """
        Convert validpairs line to Link line.

        Params:
        -------
        vpl: `ValidPairsLine`

        Examples:
        ---------
        >>> links = LinkLine()
        >>> links.vpl2links(vpl)
        >>> links
        0 Chr1 123123 0 ...
        """
        self.strandA = 0 if vpl.strand1 == '+' else 16
        self.read1_name = vpl.chr1
        self.read1_pos = vpl.pos1  
        self.order1 = 0
        self.strandB = 0 if vpl.strand2 == '+' else 16
        self.read2_name = vpl.chr2
        self.read2_pos = vpl.pos2
        self.order2 = 1


def validpairs2links(args):
    """
    %(prog)s <sample.ValidPairs> [Options]

        convert hicpro validpairs to 3d-dna pairs
    """
    p = argparse.ArgumentParser(prog=validpairs2links.__name__,
                        description=validpairs2links.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('validpairs', 
            help='validpairs from hicpro')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    vp = ValidPairs(args.validpairs)
    for i, vpl in enumerate(vp.getLine()):
        links = LinkLine()
        links.vpl2links(vpl)
        print(links, file=args.output)
        if i % 100000 == 0:
            logging.debug('Parse {} pairs'.format(i))

if __name__ == "__main__":
    validpairs2links(sys.argv[1:])
    