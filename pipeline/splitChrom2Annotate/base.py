#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
base packages for pipeline of splitHomo2Annotate
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import gzip
#import pandas as pd

from Bio import SeqIO

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

# def import_blast_table(tsv):
    
#     df = pd.read_csv(tsv, sep='\t', header=None,
#                 index_col=1, usecols=[0, 1],
#                 names=['gene', 'chrom'])

#     return df


def blast2fasta(gene_set, infasta, output):
   
    with open(output, 'w') as out: 
        extract_fasta(infasta, gene_set, out)