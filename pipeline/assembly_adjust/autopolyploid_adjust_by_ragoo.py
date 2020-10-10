#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
adjust autopolyploid chromomsome by ragoo.
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

from collections import OrderedDict
from shutil import copyfile
from TDGP.apps.grid import Cluster

header = Cluster('PBS', threads=1).header

def import_db(infile):
    
    db = OrderedDict()
    with open(infile) as fp:
        for line in fp:
            if line.strip():
                chrom, contig = line.strip().split()
                if chrom not in db:
                    db[chrom] = []
                db[chrom].append(contig)
    return db

def write_scripts(adjust_db, relatives_dir, tig_dir, 
            chr_dir, tour_path, keep=False):
    scripts = []
    for chrom in adjust_db:
        for contig in adjust_db[chrom]:
            
            if contig == "NA":
                contig = chrom[:5]
                contig_path = op.abspath(op.join(relatives_dir, contig + ".fasta"))
            else:
                contig_path = op.abspath(op.join(chr_dir, contig + ".fasta"))
            contig_fasta = op.basename(contig_path)
            
            directory = chrom[:5]
            chrom_path = op.abspath(op.join(tig_dir, chrom + ".tig.fasta"))
            chrom_fasta = chrom + ".tig.fasta"
            
            if not keep:
                os.system("rm -rf ragoo/{}/{}".format(directory, chrom))
            try:
                os.makedirs("ragoo/{}/{}".format(directory, chrom))
            except FileExistsError:
                pass
            if contig != chrom:
                cmd = header + "\n"
                cmd += "ln -s {}\n".format(chrom_path)
                cmd += "ln -s {}\n".format(contig_path)
                #cmd += "seqkit replace -p {0} -r {1} {2} > {3}\n".format(contig, chrom, contig_path, chrom + ".fasta")
                cmd += "ragoo.py {} {}\n".format(chrom_fasta, contig_fasta)
                cmd += "ragooOrder2tour.py ragoo_output/orderings/{0}_orderings.txt > {1}.tour\n".format(contig, chrom)
                script = "ragoo/{0}/{1}/run_{1}-{2}.sh".format(directory, chrom, contig)

                with open(script, 'w') as fp:
                    fp.write(cmd) 
                scripts.append(script)
                pass
            else:
                tour_file = op.join(tour_path, contig + ".tour")
                order_dir = "ragoo/{}/{}".format(directory, chrom)
                order_dir_tour = op.join(order_dir, op.basename(tour_file))
                copyfile(tour_file, order_dir_tour)
    return scripts


def run_ragoo(scripts):
    scripts = list(map(lambda x: op.split(x), scripts))
    wrkdir = os.getcwd()
    for path, script in scripts:
        os.chdir(path)
        os.system("qsub {}".format(script))
        os.chdir(wrkdir)


def autopolyploid_adjust_by_ragoo(args):
    """
    %(prog)s <adjust.table> [Options]
        run autopolyploid adjust by ragoo
    """
    p = p=argparse.ArgumentParser(prog=autopolyploid_adjust_by_ragoo.__name__,
                        description=autopolyploid_adjust_by_ragoo.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('adjustTable', 
            help='table of adjust, two columns')
    pOpt.add_argument('--relatives_dir', default='relatives_dir', 
            help='relatives chromosome fasta directory [default: %(default)s]')
    pOpt.add_argument('--chr_dir', default='chr_ref',
            help='chromosome fasta directory [default: %(default)s]')
    pOpt.add_argument('--tig_dir', default='tig_ref',
            help='contig fasta directory [default: %(default)s]')
    pOpt.add_argument('--tour_dir', default='tour', 
            help='chromosome tour file path [default: (default)s]')
    pOpt.add_argument('--qsub', action='store_true', default=False,
            help='qsub all jobs to cluster [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    adjust_db = import_db(args.adjustTable)
    scripts = write_scripts(adjust_db, args.relatives_dir, 
                                args.tig_dir, args.chr_dir, args.tour_dir)
    if args.qsub:
        run_ragoo(scripts)

if __name__ == "__main__":
    autopolyploid_adjust_by_ragoo(sys.argv[1:])