#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
%prog allPairs bam
convert the hicpro map result to binless.
"""


from __future__ import print_function
from collections import defaultdict

import os
import os.path as op
import pandas
import pysam
import sys


def read_all_pairs(all_pairs):
    header = ('chr1', 'pos1', 'strand1', 'chr2', 'pos2', 'strand2',
                'distance', 'enz1', 'enz2', 'qual1', 'qual2','nan')
    df = pandas.read_csv(all_pairs,sep='\t',header=None,index_col=0,
                    names=header)
    return df.drop('nan', axis=1)

def read_enzyme_bed(enzyme_bed):
    header = ('chrn','start', 'end', 'name', 'zero','strand')
    df = pandas.read_csv(enzyme_bed, sep='\t',header=None,index_col=3,
                    names=header)
    
    return df


def calculate_site_distance(all_pairs, enzyme_bed):
    enzyme_df = read_enzyme_bed(enzyme_bed)
    all_pairs_df = read_all_pairs(all_pairs)
    
    def calc(ids):

        pos1, chr1, pos2,chr2, enz1, enz2 = \
                all_pairs_df.loc[ids,['pos1','chr1','pos2','chr2','enz1','enz2']]
            
        data1 = enzyme_df.loc[(enzyme_df['chrn'] == chr1) & (enzyme_df['start'] <= pos1) & (enzyme_df['end'] >= pos1)]
        data2 = enzyme_df.loc[(enzyme_df['chrn'] == chr2) & (enzyme_df['start'] <= pos2) & (enzyme_df['end'] >= pos2)]
        
        start1, end1 = data1['start'][0],data1['end'][0]
        start2, end2 = data2['start'][0],data2['end'][0]
        rd1 = pos1 - start1 
        ld1 = end1 - pos1
        rd2 = pos2 - start2
        ld2 = end2 - pos2

        return rd1, ld1, rd2, ld2


    distance_dict = {}
    for ids in all_pairs_df.index:
        distance_dict[ids] = calc(ids)
    header = ['rup1','rdn1','rup2','rdn2']
    df = pandas.DataFrame(distance_dict)
    df = df.T
    df.columns = header
    
    return (pandas.concat([all_pairs_df,df],axis=1))


def read_bam(bamfile):
    db = defaultdict(lambda:[])
    sam = os.popen('samtools view {}'.format(bamfile))
    for record in sam:
        data = record.split()
        db[data[0]].append(len(data[9]))

    return db


def filter_read(db, read_list):
    read_set = set(read_list)
    def is_in(x):
        if x[0] in read_set and len(x[1]) == 2:
            return True

    return dict(filter(is_in, db.items()))

def out_tsv(allpairs, enzyme_bed, bamfile,out_csv):
    df = calculate_site_distance(allpairs, enzyme_bed)
    bam_db = read_bam(bamfile)
    read_list = list(df.index)
    db = filter_read(bam_db, read_list)
    read_len_df = pandas.DataFrame(data=db)
    read_len_df = read_len_df.T
    read_len_df.columns = ['length1','length2']
    out_df = pandas.concat([df,read_len_df],axis=1)
    
    out_header = ['chr1','pos1','strand1','length1','rup1',
                'rdn1','chr2','pos2','strand2','length2',
                'rup2','rdn2']
    out_df = out_df[out_header]
    out_df.to_csv(out_csv,sep='\t')



        




if __name__ == "__main__":
    from optparse import OptionParser
    p = OptionParser(__doc__)
    
    opts, args = p.parse_args()

    if len(args) != 4:
        sys.exit(p.print_help())

    all_pairs, enzyme_bed, bamfile, out_csv = args
    out_tsv(all_pairs, enzyme_bed, bamfile, out_csv)


