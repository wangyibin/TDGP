#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Genome analysis libraries.
"""

from __future__ import print_function

import argparse
import gzip
import logging
import numpy as np 
import os.path as op
import os
import sys

from Bio import SeqIO, SeqUtils
from collections import OrderedDict
from joblib import Parallel, delayed, Memory
from optparse import OptionParser
from TDGP.apps.base import ActionDispatcher
from TDGP.apps.base import debug, listify, check_file_exists
from TDGP.apps.base import BaseFile, Line

debug()

def main():
    actions =  (
        ('calGC', 'calculate gc content of a genome'),
        ('getTSS', 'get a bed file of TSS'),
        ('getTSSbw', 'get bigwig of TSS density per window'), 
        
        
    )

    p = ActionDispatcher(actions)
    p.dispatch(globals())



class Genome(object):
    """
    Object for genome file with some uitl function to analysis genome.
    
    Params:
    --------
    filename: `str`

    Returns:
    --------

    Examples:
    --------

    """

    def __init__(self, filename, exclude=None, 
        exclude_contig=['tig', 'Un', 'Sy', 'scaffold', 'ctg', 'Pt', 'Mt'], 
        mem_cache='.'):
        check_file_exists(filename)
        self.filename = filename
        self.exclude = listify(exclude)
        self.exclude_contig = listify(exclude_contig)
        self.getChrSizes()
        self.idx2label = dict((i, chrom) 
            for i, chrom in enumerate(self.chromLabels))
        self.label2idx = dict((chrom, i) 
            for i, chrom in enumerate(self.chromLabels))

        self.mem_cache = mem_cache
        self.memory = Memory(mem_cache, verbose=0)
        self.getGCBin = self.memory.cache(self._getGCBin)
    @property
    def handle(self):
    
        if self.filename[-3:] == ".gz":
            self._handle = gzip.open(self.filename, 'rt')
        else:
            self._handle = open(self.filename, 'r')

        return self._handle
        

    @property
    def seqs(self):
        """
        A OrderedDict of sequences.
        """
        if not hasattr(self, '_seqs'):
            self._seqs = []
            fa = SeqIO.parse(self.handle, 'fasta')
            for record in fa:
                if self.exclude:
                    if record.id in self.exclude:
                        continue 
                if self.exclude_contig:
                    for contig in self.exclude_contig:
                        if contig in record.id:
                            break
                    else:
                        self._seqs.append(record.seq)
                else:
                    self._seqs.append(record.seq)
        return self._seqs


    @property
    def chromLabels(self):
        if not hasattr(self, '_chromLabels'):
            self._chromLabels = []
            import pyfaidx
            fa = pyfaidx.Fasta(self.filename)
            for record in fa:
                if self.exclude:
                    if record.name in self.exclude:
                        continue 
                if self.exclude_contig:
                    for contig in self.exclude_contig:
                        if contig in record.name:
                            break
                    else:
                        self._chromLabels.append(record.name)
                else:
                    self._chromLabels.append(record.name)
        return self._chromLabels

    @property
    def chroms(self):
        return list(range(len(self.chromLabels)))
   
    @property
    def chromCount(self):
        return len(self.chroms)

    def getChrSizes(self):
        """
        Calculate the length of chromosome.
        """
        self.chromSizes = np.array([len(self.seqs[i]) 
                for i in range(self.chromCount)])
        return self.chromSizes


    def makeWindows(self, window):
        """
        make chromosome window

        Params:
        --------
        window: `int` window of chromosome

        Returns:
        --------
        out: `list` a list of  windows:

        Examples:
        ---------
        >>> makeWindows(10000)
        [('Chr1', 0, 100000) ...]
        """
        self.window = window
        if not hasattr(self, 'windows'):
            self.windows = OrderedDict()
            for idx, size in enumerate(self.chromSizes):
                temp = []
                chrom = self.idx2label[idx]
                for i in range(0, size + 1, window):
                    temp.append((i , i + window))
                else:
                    if temp[-1][1] > size:
                        temp[-1] = (temp[-1][0], size)
                self.windows[chrom] = temp
            self.chromBins = list(map(len, self.windows.values()))
            self.chromStartBins = np.r_[0, np.cumsum(self.chromBins[:-1])]
            self.chromEndBins = np.cumsum(self.chromBins)
            self.numBins = self.chromEndBins[-1]
            self.chromBinsDict = OrderedDict(zip(self.windows.keys(), 
                tuple(zip(self.chromStartBins, self.chromEndBins))))
        
            logging.debug('Successful makewindow')
        return self.windows


    def getGapBase(self, chrom, start, end):
        """
        Calculate the percentage of gap base number in a region
        """
        seq = self.seqs[chrom][start: end]
        if len(seq) == 0:
            return 0.0
        else:
            gap = seq.count('N') + seq.count('n')
            percent = 100.0 * gap / float(len(seq))
        
            return percent
        

    def getGC(self, chrom, start, end, correct=True):
        """
        Calculate the GC content of a sequence.
        """
        seq = self.seqs[chrom][start: end]
        
        gc = SeqUtils.GC(seq)
        
        gap = self.getGapBase(chrom, start, end) if correct \
            else 0.0

        if gap == 100.0:
            return -1.0
        else:
            corrected_gc = gc * 100.0 / (100.0 - gap)
            #logging.debug('Calculated GC content in {}:{}-{}'.format(
            #        chrom, start, end))
            return corrected_gc
    

    def _getGCBin(self, window, chr=[], correct=True, thread=24):
        """
        Calculate GC content of a series of windows, and return a OrderedDict

        Params:
        --------
        window: `int` window of bin
        chr: `list` default: `[]`
        thread: `int` thread of parallel running default: `24`
        Returns:
        --------
        out: `list` and gc store in array-like

        Examples:
        --------
        >>> getGCbin(1000000)
        [[0.5, 0.2, 0.5 ...], ...]

        """
        self.gcBin = []
        chroms = listify(chr) if chr else self.chromLabels
        _chromsidx = [self.label2idx[i] for i in chroms]
        """
        def subgc(chrom):
            chromWindow = int(self.chromSizes[chrom] // self.window) + 1
            _gc = np.ones(chromWindow, dtype=np.float)
            for i in range(chromWindow):
                _gc[i] = self.getGC(chrom,  i*self.window, 
                                    (i+1)*self.window, correct=correct)
                
            return _gc
        res = Parallel(thread)(delayed(subgc)(args)
                for args in _chromsidx)
        """
        for chrom in _chromsidx:
            chromWindow = int(self.chromSizes[chrom] // self.window) + 1
            self.gcBin.append(np.ones(chromWindow, dtype=np.float))
            for i in range(chromWindow - 1):
                 
                self.gcBin[chrom][i] = self.getGC(chrom,  i*self.window, 
                                        (i+1)*self.window, correct=correct)
            else:
                self.gcBin[chrom][chromWindow-1] = self.getGC(chrom, 
            (chromWindow-1)*self.window, chromWindow*self.window, correct=correct)
        
        return self.gcBin
        
    def clearCache(self):
        """
        clear Memory cache data in the `{}`.
        """.format(self.mem_cache)

        if hasattr(self, 'memory'):
            self.memory.clear()

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



## out command ##
def calGC(args):
    """
    %(prog)s in.fasta [Options]
        calculate GC content.
    """

    p = p=argparse.ArgumentParser(prog=calGC.__name__,
                        description=calGC.__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fasta', help='input fasta file')
    pOpt.add_argument('-w', '--window', type=int, default=100000,
            help='size of calculation window [default: %default]')
    pOpt.add_argument('-o', '--out', type=argparse.FileType('w'), 
            default=sys.stdout, help='output file [default: %(default)s]')
    pOpt.add_argument('--exclude', nargs="*", default=[],
            help='exclude these chromosome [default: %(default)s]')
    pOpt.add_argument('--exclude_contig', nargs='*', 
            default=['tig', 'scafflod', 'Un', 'Sy', 'Mt', 'Pt'], 
            help='exclude these chromosome if it contain these string'
                ' [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    args = p.parse_args(args)
    genome = Genome(args.fasta, exclude=args.exclude, 
            exclude_contig=args.exclude_contig)

    genome.makeWindows(args.window)
    gcBins = genome.getGCBin(args.window, correct=True)
    
    for i, chrom in enumerate(genome.chromLabels):
        for n, windows in enumerate(genome.windows[chrom]):
            
            print("\t".join((chrom, "\t".join(map(str, windows)), 
                str(gcBins[i][n]))), file=args.out)
    


if __name__ == "__main__":
    main()

