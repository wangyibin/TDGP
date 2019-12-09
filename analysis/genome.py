#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Genome analysis libraries.
"""

from __future__ import print_function

import gzip
import logging
import os.path as op
import os
import sys

from optparse import OptionParser
from TDGP.apps.base import ActionDispatcher
from TDGP.apps.base import debug, listify, check_file_exists
from TDGP.apps.base import BaseFile, Line

debug()

def main():
    actions =  (
        ('getTSS', 'get a bed file of TSS'),
        ('getTSSbw', 'get bigwig of TSS density per window'), 
        
        
    )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def guess_filetype(infile):
    """
    To guess infile is gff or bed.
    """
    gff_end = ['gff3', 'gff3.gz', 'gff', 'gff.gz']
    bed_end = ['bed', 'bed.gz']
    for type_ in gff_end:
        if infile.endswith(type_):
            return 'gff'
    
    for type_ in bed_end:
        if infile.endswith(type_):
            return 'bed'
    
    return None

def must_open(infile, mode='r'):
    """
    To parse gz or not gz file, and return a file handle.
    """
    if infile[-3:] == ".gz":
        handle = gzip.open(infile, mode)
    else:
        handle = open(infile, mode)
    return handle


def getTSS(args):
    """
    %prog infile [Options]
        To get TSS from gff or bed.
    """

    p = OptionParser(getTSS.__doc__)
    p.add_option('--type', default='gene',
            help='the type of sequence [default: %default]')
    p.add_option('-o', '--out', default=sys.stdout,
            help='output file. [default: stdout]')
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(p.print_help())
    
    infile, = args
    check_file_exists(infile)
    out = opts.out
    filetype = guess_filetype(infile)
    if not filetype:
        logging.error('Input filetype must a gff or bed')
        sys.exit()
    
    if filetype == "gff":
        with must_open(infile, 'r') as fp:
            for line in fp:
                if line.startswith("#"):
                    continue
                line_list = line.strip().split('\t')
                (chrom, _, type_, start, end, _, strand, 
                    _, info) = line_list[:9]
                
                if type_ == opts.type:
                    print('\t'.join(map(str, (chrom, start, 
                            int(start) + 1))), file=out)
    
    elif filetype == 'bed':
        with must_open(infile, 'r') as fp:
            for line in fp:
                if line.startswith("#"):
                    continue
                line_list = line.strip().split()
                chrom, start, end = line_list[:3]
                print('\t'.join(map(str, (chrom, start, 
                            int(start) + 1))), file=out)
    
    out = 'stdout' if isinstance(out, file) else out
    logging.debug('Done, output is in `{}`'.format(out))


def getTSSbw(args):
    """
    %prog <tss.gff/bed> <chrom.sizes> <out_prefix> [options]
        To obtain a bedgraph file of tss sites 
            density of per windows
    """

    p = OptionParser(getTSSbw.__doc__)
    p.add_option('-w', '--window', type=int, default=1000,
            help='the window of tss density calculation')
    p.add_option('-o', '--out', default=sys.stdout,
            help='output [default: stdout]')
    p.add_option('--qsub', default=False, action='store_true',
            help='if qsub to sge [default: %default]')
    opts, args = p.parse_args(args)
    if len(args) != 3:
        sys.exit(p.print_help())
    
    gff, chrom_sizes, sample = args
    window = opts.window
    check_file_exists(chrom_sizes)
    command = 'qsub -pe mpi 1 -cwd -j y -S /bin/bash' if opts.qsub \
            else 'sh'
        
    cmd = 'python -m TDGP.analysis.genome getTSS {} > \
        {}.gene.tss.bed\n'.format(gff, sample)
    cmd += 'bedtools makewindows -g {} -w {} > {}.{}.window\n'.format(
                chrom_sizes, window, sample, window)
    cmd += 'bedtools intersect -a {sample}.{window}.window -b \
            {sample}.gene.tss.bed -c | sort -k1,1 -k2,2n > \
                {sample}.gene.tss.{window}.bg\n'.format(sample=sample, 
                window=window)
    cmd += 'bedGraphToBigWig {sample}.gene.tss.{window}.bg {sizes} \
            {sample}.gene.tss.{window}.bw\n'.format(sample=sample, 
                window=window, sizes=chrom_sizes)
    with open('run_{}_tss.sh'.format(sample), 'w') as out:
        out.write(cmd)
    
    os.system('{} run_{}_tss.sh'.format(command, sample))
    logging.debug('Successful')





if __name__ == "__main__":
    main()






