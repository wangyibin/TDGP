#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
convert agpfile to clusters file
"""

import argparse
import logging
import os
import os.path as op
import sys

import pandas as pd 

from pyfaidx import Fasta

AGP_NAMES_tig = ['chrom', 'start', 'end', 'number',
                 'type', 'id', 'tig_start', 'tig_end', 'orientation']
AGP_NAMES_gap = ['chrom', 'start', 'end', 'number',
                 'type', 'length', 'type', 'linkage', 'evidence']
def import_agp(agpfile, split=True):
    """
    import agp file and return a dataframe
    """
    df = pd.read_csv(agpfile, sep='\t', comment='#',
                     header=None, index_col=None,)
    logging.info('load agp file: `{}`'.format(agpfile))
    
    if split:
        tig_df = df[df[4] == 'W']
        gap_df = df[df[4] == 'U']
        tig_df.columns = AGP_NAMES_tig
        gap_df.columns = AGP_NAMES_gap
        tig_df = tig_df.astype(
            {'chrom': 'category', 'orientation': 'category'})
        gap_df = gap_df.astype({'chrom': 'category'})

        tig_df.set_index('chrom', inplace=True)
        gap_df.set_index('chrom', inplace=True)
        
        return tig_df, gap_df
    else:
        return df


def agp2contig(agp, contig, outdir):
    agp_df, _ = import_agp(agp)
    
    # remove contig
    agp_df.reset_index(inplace=True)
    agp_df = agp_df[agp_df['chrom'] != agp_df['id']]
    agp_df['chrom'] = agp_df['chrom'].astype('category')
    cluster_df = agp_df.groupby('chrom')['id'].apply(lambda x: list(x))
    
    fasta = Fasta(contig)

    if not op.exists(outdir):
        os.makedirs(outdir)


    for i, cluster in cluster_df.iteritems():
        if not cluster:
            continue
        with open(f'{outdir}/{i}.contig.fasta', 'w') as out:
            for ID in cluster:
                seq = fasta[ID]
                out.write(f">{ID}\n{seq}\n")
       

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('agp', help='agp file')
    pReq.add_argument('contig', help='contig-level fasta')
    pOpt.add_argument('-o', '--outdir', default='agp2contig_dir',
            help='output dir [default: %(default)s')
    
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)


    agp = args.agp
    contig = args.contig
    agp2contig(agp, contig, args.outdir)

if __name__ == "__main__":
    main(sys.argv[1:])
