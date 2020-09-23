#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
split genome to homologous groups
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys
import gzip

from Bio import SeqIO
from collections import OrderedDict
from multiprocessing import Pool

def extract_fasta(infasta, gene_set, output_handle):
    """
    extract sequences by list
    
    Params:
    --------
    infasta: `str` 
            fasta file
    gene_set: `set` or `list-like` 
            gene lists
    output_handle: `handle` of output

    Returns:
    --------
    write extracted fasta to a fasta file

    Examples:
    --------
    >>> extract_fasta("sample.fasta", gene_set, output_handle)
    """
    if infasta[-2:] == "gz":
        fp = gzip.open(infasta)
    else:
        fp = open(infasta)
    
    fa = SeqIO.parse(fp, 'fasta')
    for record in fa:
        if record.id in gene_set:
            SeqIO.write(record, output_handle, 'fasta')


def read_groups(groups):
    groups_db = OrderedDict()
    with open(groups) as fp:
        for line in fp:
            group, chrom = line.strip().split()
            if group not in groups_db:
                groups_db[group] = []
            groups_db[group].append(chrom)

    return groups_db

def extract_fasta_by_chroms(args):
    fasta, chroms, output = args
    with open(output, 'w') as out:
        extract_fasta(fasta, chroms, out)
   

def splitGenomebyGroups(args):
    """
    %(prog)s <genome.fasta> <groups.db> [Options]
        split genome by groups
    """
    p = p=argparse.ArgumentParser(prog=splitGenomebyGroups.__name__,
                        description=splitGenomebyGroups.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('fasta', 
            help='genome fasta')
    pReq.add_argument('groups_db', help='groups db file')
    pOpt.add_argument('-t', '--thread', type=int, default=4,
            help='number of thread [default: %(default)s]')
    #pOpt.add_argument('--removeContig', action='store_true', default=False,
     #       help='remove contig [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)

    groups_db = read_groups(args.groups_db)
    task_list = []
    for group in groups_db:
        try:
            os.makedirs(group)
        except FileExistsError:
            pass
        chroms = groups_db[group]
        output = "{0}/{0}.fasta".format(group)
        task_list.append((args.fasta, chroms, output))
    
    pool = Pool(args.thread)
    pool.map(extract_fasta_by_chroms, task_list)
    pool.close()

        
    
if __name__ == "__main__":
    splitGenomebyGroups(sys.argv[1:])



