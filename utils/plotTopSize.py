#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
%(prog)s <tad.txt> [numbers] [options]
    plotting the top size of tads 
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from TDGP.apps.grid import CMD

def topSize(infile, minimum=False, num=10):

    tailorhead = 'head' if minimum else 'tail'
    cmd = ("awk '{{print $0,$3-$2}}' OFS='\t' {infile} | "
            "sort -rk5,5n | {tailorhead} -{num}".format(infile=infile, 
            tailorhead=tailorhead, num=num))
    
    return os.popen(cmd)

def main(args):
    
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('infile', 
            help='tad file from hitad')
    pReq.add_argument('track', 
            help='track file of pyGenomeTracks')
    pOpt.add_argument('-o', '--outdir', default='./', 
            help='output directory [default: %(default)s]')
    pOpt.add_argument('--slide_window', default=100000, type=int,
            help='the slide window of tad down and up [default: %(default)s]')
    pOpt.add_argument('--minimum', action='store_true', default=False,
            help='plotting the minimum [default: %(default)s]')
    pOpt.add_argument('-n', '--numbers', default=10, type=int,
            help='numbers of tads to plot [default: %(default)s]')
    pOpt.add_argument('-t', '--threads', default=4, type=int,
            help='number of program threads to launch [default:s %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    
    if args.outdir != './':
        try:
            os.makedirs(args.outdir)
        except FileExistsError:
            logging.warning('the directory of `{}` is already exists'.format(args.outdir))

    slide_window = args.slide_window
    cmds = []
    with topSize(args.infile, args.minimum, args.numbers) as fp:
        for line in fp:
            chrom, start, end = line.strip().split()[:3]
            slide_start = int(start) - slide_window
            slide_end = int(end) + slide_window
            regions = "{}:{}-{}".format(chrom, slide_start, slide_end)
            outfile = "{}/{}:{}-{}.png".format(args.outdir, chrom, start, end)

            cmd = 'pyGenomeTracks --tracks {} -o {} --region {}'.format(args.track, outfile, regions)
            cmds.append(cmd)
    
    CMD(cmds, threads=args.threads)

if __name__ == "__main__":
    main(sys.argv[1:])