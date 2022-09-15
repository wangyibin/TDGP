#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
convert agp file to juicerbox assemble tools assembly file.
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys
import pandas as pd

from collections import OrderedDict

AGP_NAMES_tig = ['chrom', 'start', 'end', 'number',
                 'type', 'id', 'tig_start', 'tig_end', 'orientation']
AGP_NAMES_gap = ['chrom', 'start', 'end', 'number',
                 'type', 'length', 'type', 'linkage', 'evidence']


def import_agp(agpfile, split=True):
    """
    import agp file and return a dataframe
    """
    df = pd.read_csv(agpfile, sep='\t', comment='#',
                     header=None, index_col=None)
    if split:
        tig_df = df[df[4] == 'W']
        gap_df = df[df[4] == 'U']
        tig_df.columns = AGP_NAMES_tig
        gap_df.columns = AGP_NAMES_gap
        tig_df.set_index('chrom', inplace=True)
        gap_df.set_index('chrom', inplace=True)

        return tig_df, gap_df
    else:
        return df


def agp2assembly(args):
    """
    %(prog)s <groups.agp> [Options]

        convert agp to assembly, for juicerbox assemble tool adjust.
    """
    p = p = argparse.ArgumentParser(prog=agp2assembly.__name__,
                                    description=agp2assembly.__doc__,
                                    conflict_handler='resolve',
                                    formatter_class=argparse.RawTextHelpFormatter)
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('agpfile', help='agpfile from allhic or other ways')
    pOpt.add_argument('--add_gap', default=False, action="store_true",
                      help='add gap into assembly [default: %(default)s]')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
                      default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
                      help='show help message and exit.')

    args = p.parse_args(args)
    agpfile = args.agpfile
    tig_df, gap_df = import_agp(agpfile)
    gap_length = gap_df.iloc[0]['length']
    chrom_matrix_db = OrderedDict()
    tig_list = []
    for i, item in enumerate(tig_df.iterrows(), 1):
        chrom = item[0]
        tig = item[1].id
        tig_length = item[1].tig_end
        orientation = item[1].orientation

        tig_list.append([tig, i, tig_length])
        if chrom not in chrom_matrix_db:
            chrom_matrix_db[chrom] = []
        orientation = orientation if orientation == '-' else ""
        chrom_matrix_db[chrom].append(orientation + str(i))
    else:
        hic_gap_number = i + 1

    for item in tig_list:
        print(">{}".format(" ".join(map(str, item))), 
                            file=args.output)
    else:
        if args.add_gap:
            print(">{}".format(" ".join(['hic_gap_{}'.format(
                hic_gap_number), str(hic_gap_number), str(gap_length)])), 
                file=args.output)
            _gap = " {} ".format(hic_gap_number)
        else: 
             _gap = " "
        for chrom in chrom_matrix_db:
            print(_gap.join(
                map(str, chrom_matrix_db[chrom])), file=args.output)


if __name__ == "__main__":
    agp2assembly(sys.argv[1:])
    
