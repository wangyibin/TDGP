#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
make a paierd bam which from allhic pipeline

"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys
import pysam


class LinkLine(object):
    """
    LinkLine object for 3d-dna links line.
    Example
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


def bam2links(args):
    """
    %(prog)s <sample.sorted.bam> [Options]
        
    """
    p = p = argparse.ArgumentParser(prog=bam2links.__name__,
                                    description=bam2links.__doc__,
                                    formatter_class=argparse.RawTextHelpFormatter,
                                    conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('bamfile',
                      help='bam file from allhic')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
                      default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
                      help='show help message and exit.')

    args = p.parse_args(args)
    bamfile = args.bamfile
    handle = pysam.Samfile(bamfile, 'rb')
    pre_name = ""
    count = 0
    links_num = 0
    for read in handle.fetch(until_eof=True):
        links = LinkLine()
        count += 1
        if read.query_name == pre_name:
            continue
        pre_name = read.query_name
        if not read.is_unmapped and not read.is_unmapped \
                and read.is_paired:
            links.strandA = 16 if read.is_reverse else 0
            links.read1_name = read.reference_name
            links.read1_pos = read.reference_start
            links.strandB = 16 if read.mate_is_reverse else 0
            links.read2_name = read.next_reference_name
            links.read2_pos = read.next_reference_start

            if read.reference_id <= read.next_reference_id:
                links.order1 = 0
                links.order2 = 1
            else:
                links.order1 = 1
                links.order2 = 0

            links_num += 1
            print(links, file=args.output)

            if count % 100000 == 0:
                logging.info('Parse {} reads'.format(count))

    else:
        logging.info('Done...Output {} pairs in {}'.format(
            links_num, args.output.name))


if __name__ == "__main__":
    bam2links(sys.argv[1:])
