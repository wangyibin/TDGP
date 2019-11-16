#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
%prog hicpro_output [Options]
    
    To stat Hi-C libraries quality from HiCPro results.
"""
from __future__ import print_function

import glob
import gzip
import os
import os.path as op
import pandas as pd
import sys

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(__name__)

def fqSize(fastq):
    if fastq[-3:] == ".gz":
        open = gzip.open
    
    reads_num = 0
    base_num = 0
    with open(fastq) as fp:
        for i, line in enumerate(fp, 1):
            if i % 4 == 2:
                base_num += len(line.strip())
                reads_num += 1

    return reads_num, base_num
header = ("Lib", "Data(bp)", "Total Read Pairs", "Mapped Reads Pairs", 
            "Unique Mapped Read Pairs", "Valid Interaction Pairs")


def stat(hic_output, species, sample, read_length=150, lib_func='lambda x: x'):
    lib_func = eval(lib_func)
    lib = lib_func(sample)

    pairstat_file = "{}/bowtie_results/bwt2/{}/{}_reference.fasta.bwt2pairs.pairstat".format(
            hic_output, species, sample)
    df = pd.read_csv(pairstat_file, sep="\t", header=None, index_col=0)
    total_read_pairs = df.loc["Total_pairs_processed"][1]
    data = total_read_pairs * read_length * 2
    mapped_reads_pairs = df.loc["Total_pairs_processed"][1] - df.loc["Unmapped_pairs"][1]
    unique_mapped_read_pairs = df.loc["Unique_paired_alignments"][1]
     
    valid_stat_file = "{}/hic_results/data/{}/{}_reference.fasta.bwt2pairs.RSstat".format(
            hic_output, species, sample)
    df = pd.read_csv(valid_stat_file, sep="\t", header=None, index_col=0, comment="#")
    valid_interaction_pairs = df.loc['Valid_interaction_pairs'][1]
    result_list = (lib, data, total_read_pairs, mapped_reads_pairs, 
            unique_mapped_read_pairs, valid_interaction_pairs)
    result = dict(zip(header, result_list))

    return result


def filter_inputfile(inputfiles, species):
    for inputfile in inputfiles:
        if op.dirname(inputfile) == species:
            return True
        else:
            return False


def valid_rate(df):
    return "{:,.0f}({:.2%})".format(df['Valid Interaction Pairs'],df['Valid Interaction Pairs']/df['Unique Mapped Read Pairs'])

def unique_rate(df):
    return "{:,.0f}({:.2%})".format(df['Unique Mapped Read Pairs'], df['Unique Mapped Read Pairs']/df['Mapped Reads Pairs'])

def mapped_rate(df):
    return "{:,.0f}({:.2%})".format(df['Mapped Reads Pairs'], df['Mapped Reads Pairs']/df['Total Read Pairs'])





def main(hic_output, suffix="_R1.fastq.gz", read_length=150, lib_func='lambda x: x'):
    inputfiles_path = glob.glob('{}/inputfiles_*.txt'.format(hic_output))[0]
    inputfiles = [i.strip() for i in open(inputfiles_path) if i.strip()]
    species_list = list(set(map(op.dirname, inputfiles)))

    
    for species in species_list:
        stat_result = []
        for sample in inputfiles:
            if op.dirname(sample) != species:
                continue
            
            sample = op.basename(sample).replace(suffix, "")
            result = stat(hic_output, species, sample, read_length, lib_func)
            stat_result.append(result)

        df = pd.DataFrame(stat_result)
        df.set_index(df['Lib'], inplace=True)
        df.drop('Lib', axis=1, inplace=True)
        
        df1 = df.copy()
        df1 ['Valid Interaction Pairs'] = df.apply(valid_rate, axis=1)
        df1 ['Unique Mapped Read Pairs'] = df.apply(unique_rate, axis=1)
        df1['Mapped Reads Pairs'] = df.apply(mapped_rate, axis=1)
        
        df1['Total Read Pairs'] = df1['Total Read Pairs'].astype(int).map(lambda x: format(x, ','))
        df1['Data(bp)'] = df1['Data(bp)'].astype(int).map(lambda x: format(x, ','))
        df1.loc['Total'] = df.sum().astype(int).map(lambda x: format(x, ','))
        
        
        df1 = df1[list(header[1:])]   
        df1.to_csv('{}_hic_lib_stat.csv'.format(species), sep='\t')
        df1.to_excel('{}_hic_lib_stat.xls'.format(species), sheet_name=species)




if __name__ == "__main__":
    from optparse import OptionParser

    p = OptionParser(__doc__)
    p.add_option('--suffix', default="_R1.fastq.gz",
            help='The suffix of read1 fastq [default: %default]')
    p.add_option('--read_length', default=150, type=int,
            help='The length of reads ')
    p.add_option('--lib_func', default='lambda x:x',
            help='The function of lib name')

    opts, args = p.parse_args()
    if len(args) == 0:
        sys.exit(p.print_help())

    hic_output, = args
    suffix = opts.suffix
    log.info('Start stat ...')
    main(hic_output, opts.suffix, opts.read_length, opts.lib_func)
    log.info('Done')
