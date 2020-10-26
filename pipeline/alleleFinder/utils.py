#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
util function for alleleFinder
"""
from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys
import string
import pandas as pd

from Bio import SeqIO
from itertools import combinations


PBS_HEADER = """#!/bin/bash
#PBS -j oe {}
#PBS -q {}
#PBS -V 
#PBS -l nodes=1:ppn={} {}
if [[ ! -z $PBS_O_WORKDIR ]]; then
    cd $PBS_O_WORKDIR
fi
"""

SGE_HEADER = """#!/bin/bash
#$ -j y
#$ -S /bin/bash
#$ -cwd {}
#$ -q {}
#$ -pe mpi {} {}
"""
class Cluster(object):
    """
    class of cluster operation
        in order to execute successful should set 
                the `CLUSTER` variants into ENV
    Params:
    --------
    
    Returns:
    ---------
    out: `str`: CLUSTER

    Functions:
    ---------
    get_header: get the header of cluster system with parameters
    get_raw_header: get the raw of cluster system with parameters

    """
    
    def __init__(self, cluster=None, 
                    name=None, queue=None, 
                    threads=1, array=None):
        self.CLUSTER = cluster if cluster else None
        if not self.CLUSTER:
            self.get()
        self.get_header(name, queue, threads, array)
        self.get_raw_header()

    def get(self):
        """
        To obtain the environment of `CLUSTER`,
            if not found will be set default `SGE`.
        """
        try:
            self.CLUSTER = os.environ['CLUSTER']
        except KeyError:
            self.CLUSTER = 'SGE'
            logging.warning('There is not environment `CLUSTER` in PATH')

        return self.CLUSTER


    def get_header(self, name=None, queue=None, 
                        threads=1, array=None):
        """
        According to the environment of `CLUSTER` to 
            return a header of cluster system
        """
        if self.CLUSTER.upper() == "SGE":
            name = "\n#$ -N " + name  if name else ""
            queue = queue if queue else "all.q"
            array = "\n#$ -t " + array if array else ""
            self.header = SGE_HEADER.format(name, queue, threads, array)   
        elif self.CLUSTER.upper() == "PBS":
            name = "\n#PBS -N " + name if name else ""
            queue = queue if queue else "workq"
            array = "\n#PBS -J " + array if array else ""
            self.header = PBS_HEADER.format(name, queue, threads, array)
        elif self.CLUSTER.upper() == "TORQUE":
            name = "\nPBS -N " + name if name else ""
            queue = queue if queue else "share"
            array = "\n#PBS -J " + array if array else ""
            self.header = PBS_HEADER.format(name, queue, threads, array)
        else:
            logging.warning("there is not of header "
                            "of cluster:`{}`".format(self.CLUSTER))
            sys.exit()
        return self.header

    def get_raw_header(self):
        """
        According to the environment of `CLUSTER` to 
            return a header of cluster system
        """
        if self.CLUSTER.upper() == "SGE":
            self.raw_header = SGE_HEADER
        elif self.CLUSTER.upper() == "PBS":
            self.raw_header = PBS_HEADER
        elif self.CLUSTER.upper() == "TORQUE":
            self.raw_header = PBS_HEADER
        else:
            logging.warning("there is not of header "
                            "of cluster:`{}`".format(self.CLUSTER))
            sys.exit()
        return self.raw_header


    def __str__(self):
        return self.CLUSTER

    __retr__ = __str__

def pairs(samples):
    samples = list(combinations(samples, 2))

    return list(map(lambda x:"{}-{}".format(x[0], x[1]), samples))

def run_synteny(pairs):
    sample1, sample2 = pairs.split("-")
    os.mkdir(pairs)
    os.chdir(pairs)
    header = Cluster().header
    os.system('ln -s ../data/{}* .'.format(sample1))
    os.system('ln -s ../data/{}* .'.format(sample2))
    script = 'run_{}.sh'.format(pairs)
    with open(script, 'w') as out:
        cmd = header + "\n"
        cmd += 'python -m jcvi.compara.catalog ortholog {} {} --no_strip_names --nochpf \n'.format(sample1, sample2)
        out.write(cmd)
    os.system('qsub {}'.format(script))
    os.chdir('../')


def split_bed_by_chr(bedfile, outdir='./', prefix=None):
    """
    split bed by chromosome name
    """
    # if prefix is None:
    #     prefix = op.splitext(op.basename(bedfile))[0]
           
    df = pd.read_csv(bedfile, sep='\t', header=None, index_col=None)
    for chrom, data in df.groupby(0):
        out = '{}/{}.bed'.format(outdir, chrom)
        data.to_csv(out, sep='\t',
                    header=None, index=None)
        logging.debug('output new bed in `{}`'.format(out))


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


def create_empty_allele_table(homoTags):
    """
    create a empty dataframe for allele table

    Params:
    --------
    homoTags: `list` or `array-like` list for homologs  tags
                e.g. ['g1', 'g2', 'g3', 'g4']
    """

    i = len(homoTags)
    columns = []
    for s in string.uppercase[:i]:
        columns.append('gene', + i)
    
    df = pd.DataFrame(columns=columns)

    return df

def import_anchor(anchor):
    """
    import anchor file
    """

    df = pd.read_csv(anchor, sep='\t', header=None, 
            index_col=None, comment="#", usecols=[0, 1])
    
    return df