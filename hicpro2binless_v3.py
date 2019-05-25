#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
%prog allPairs enzyme.bed bam out.csv
convert the hicpro map result to binless.
"""


from __future__ import print_function
from collections import defaultdict

import multiprocessing
import os
import os.path as op
import pandas
import pysam
import sys


def split_pairs(all_pairs):

    df = pandas.read_csv(all_pairs, chunksize=1000000)
    
    chunk_df_list = []
    for n,chunk in enumerate(df):
        #chunk.to_csv('split_pairs/chunk{}.csv'.format(n), header=False, index=False)
        chunk_df_list.append(chunk)
    print(chunk_df_list)
    return chunk_df_list
    


def read_pairs(chunk_pairs):
    header = ('ids','chr1', 'pos1', 'strand1', 'chr2', 'pos2', 'strand2',
                'distance', 'enz1', 'enz2')
    df = pandas.read_csv(chunk_pairs,sep='\t',header=None,index_col=[1,4,0], 
                names=header, usecols=[0,1,2,3,4,5,6,7,8,9], chunksize=100000)
    
    return df

def read_enzyme_bed(enzyme_bed):
    header = ('chrn','start', 'end', 'name')
    df = pandas.read_csv(enzyme_bed, sep='\t',header=None,index_col=[0,3],
                    names=header, usecols=[0,1,2,3])
    
    return df


def enzyme_name_enlarge(name, num, max_enz):
    name = name.rsplit("_",1)
    max_num = int(max_enz.rsplit("_", 1)[-1])
    up_num = int(name[-1]) - num if int(name[-1]) - num > 0 else 0
    dn_num = int(name[-1]) + num if int(name[-1]) - num < max_num else max_num
    name_list = list(map(lambda x: name[0] + "_" + str(x), 
                    (i for i in range(up_num, dn_num + 1)) ))
    return name_list
    


def calculate_site_distance(all_pairs_df, enzyme_df):

    
    def calc(ids):

        pairs_df = all_pairs_df.loc[(slice(None),slice(None),ids),['pos1','pos2','enz1','enz2']]
        for k,v in pairs_df.iterrows():
            chr1,chr2,ids = k
            pos1, pos2, enz1, enz2 = v
            break
        enzyme_df1 = enzyme_df.loc[chr1]
        enzyme_df2 = enzyme_df.loc[chr2]
        enz1_list = enzyme_name_enlarge(enz1, 10, enzyme_df1.index[-1])
        enz2_list = enzyme_name_enlarge(enz2, 10, enzyme_df2.index[-1])
        enzyme_df1 = enzyme_df1.loc[enz1_list]
        enzyme_df2 = enzyme_df2.loc[enz2_list]
        data1 = enzyme_df1.loc[(enzyme_df1['start'] <= pos1) & (enzyme_df1['end'] >= pos1)]
        data2 = enzyme_df2.loc[(enzyme_df2['start'] <= pos2) & (enzyme_df2['end'] >= pos2)]
        
        start1, end1 = data1['start'][0],data1['end'][0]
        start2, end2 = data2['start'][0],data2['end'][0]
        rd1 = pos1 - start1 
        ld1 = end1 - pos1
        rd2 = pos2 - start2
        ld2 = end2 - pos2
        return rd1, ld1, rd2, ld2

    distance_dict = {}
    for ids in all_pairs_df.index.levels[2]:
        distance_dict[ids] = calc(ids)
    header = ['rup1','rdn1','rup2','rdn2']
    df = pandas.DataFrame(distance_dict)
    df = df.T
    df.columns = header
    all_pairs_df.reset_index(inplace=True)
    all_pairs_df.set_index(['ids'],inplace=True)
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

def out_tsv(allpairs, enzyme_bed, bamfile,out_csv, thread=12):
    split_pairs_df = read_pairs(allpairs)
    enzyme_df = read_enzyme_bed(enzyme_bed)
    pool = multiprocessing.Pool(thread)
    task_list = [ (df, enzyme_df) for df in split_pairs_df ]
    result_df_list = []
    for task in task_list:
        result_df_list.append(pool.apply_async(calculate_site_distance, task))
    pool.close()
    pool.join()
    result_df_list = [ i.get() for i in result_df_list ]
    df = pandas.concat(result_df_list, axis=1)
    bam_db = read_bam(bamfile)
    read_list = list(df.index)
    db = filter_read(bam_db, read_list)
    read_len_df = pandas.DataFrame(data=db)
    read_len_df = read_len_df.T
    read_len_df.columns = ['length1','length2']
    #df.set_index('ids')
    out_df = pandas.concat([df,read_len_df],axis=1,join='inner')

    
    out_header = ['chr1','pos1','strand1','length1','rup1',
                'rdn1','chr2','pos2','strand2','length2',
                'rup2','rdn2']
    out_df = out_df[out_header]
    out_df.to_csv(out_csv,sep='\t')



if __name__ == "__main__":
    from optparse import OptionParser
    p = OptionParser(__doc__)
    
    p.add_option('-t','--thread',dest='thread',type=int,default=4,
                    help='the thread numbe [default: %default]')
    opts, args = p.parse_args()
    if len(args) != 4:
        sys.exit(p.print_help())

    all_pairs, enzyme_bed, bamfile, out_csv = args
    out_tsv(all_pairs, enzyme_bed, bamfile, out_csv,opts.thread)
    #print(read_enzyme_bed(enzyme_bed))
    #print(read_all_pairs(all_pairs).loc[('ChrUn','ChrUn')])
    #print(enzyme_name_enlarge('HIC_Chr1_10',10,'HIC_Chr1_15'))
