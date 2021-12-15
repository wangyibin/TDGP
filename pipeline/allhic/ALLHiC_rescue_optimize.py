#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
rescue and optimize in a scaffolded group, and rescue unanchor contigs.
   
    -contig1-contig2-contig3-
    contig4
    
    contig4 -> contig2 and linked contig3
    insert contig4 into contig2 and contig3

"""

import argparse
import logging
import os
import os.path as op
import sys

import cooler 
import numpy as np
import pandas as pd 

from collections import OrderedDict


def read_tour(tourfile):
    
    with open(tourfile, 'r') as fp:
        for line in fp:
            pass
        last_line = line

    tour_db = OrderedDict()
    line_list = last_line.strip().split()
    for contig in line_list:
        contig, orientation = contig[:-1], contig[-1]
        tour_db[contig] = orientation

    return tour_db



def add_unanchor_contig(cool, tour, unanchors):

    matrix = cool.matrix(balance=False, sparse=True)
    df = pd.Series.sparse.from_coo(matrix[:])
    df = df.reset_index()

    chrom_sizes = cool.chromsizes
    chrom_bins = cool.bins()[:]
    chrom2idx = dict(zip(chrom_bins.chrom, chrom_bins.index))
    tour_list = list(tour.keys())
    tour_idx = [chrom2idx[i] for i in tour_list]
    anchored_df = df[df['level_1'].isin(tour_idx)]
    
    anchored_df.set_index(['level_0', 'level_1'], inplace=True)
    #df.set_index(['level_0', 'level_1'], inplace=True)
    unanchors_idx = [chrom2idx[i] for i in unanchors]
    unanchor_signals = []
    insert_pos_db = {}
    for idx1 in unanchors_idx:
        try:
            tmp = anchored_df.loc[idx1].sort_values(by=0, ascending=False).head(1)
        except KeyError:
            continue
        idx2 = tmp.index[0]
        count = tmp.values[0][0]
        max_pairs = (idx1, idx2, count)
        
        tour_map = tour_idx.index(idx2)

        if tour_map == 0:
            insert_pos_db[idx1] = (-1, 0)
        elif tour_map == len(tour_idx):
            insert_pos_db[idx1] = (len(tour_idx)- 1, len(tour_idx))
        else:
            up_flag, down_flag = True, True
            tour_up, tour_down = tour_map - 1, tour_map + 1
            up_idx, down_idx = tour_idx[tour_up], tour_idx[tour_down]
            
            try:
                up_signal = anchored_df.loc[(idx1, up_idx)].values[0]
            except KeyError:
                up_flag = False
            
            try:
                down_signal = anchored_df.loc[(idx1, down_idx)].values[0]
            except KeyError:
                down_flag= False
            

            if up_flag and down_flag:
                if up_signal >= down_signal:
                    down_flag = False
                else:
                    up_flag = False
           
            elif up_flag and not down_flag:
                insert_pos_db[idx1] = (tour_up, tour_map) 
            elif not up_flag and down_flag:
                insert_pos_db[idx1] = (tour_map, tour_down)
            else:
                insert_pos_db[idx1] = (tour_up, tour_map)
        unanchor_signals.append(max_pairs)

    insert_pos_db = dict(sorted(insert_pos_db.items(), key=lambda x: x[1][0]))
    for idx in insert_pos_db:
        pos = insert_pos_db[idx]
        print(pos)

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('cool', 
            help='whole contig cool file')
    pReq.add_argument('tour', 
            help='tour file')
    pReq.add_argument('unanchors',
            help='unanchor contigs list')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    coolfile = args.cool

    cool = cooler.Cooler(coolfile)

    tour_db = read_tour(args.tour)
    unanchors = [i.strip() for i in open(args.unanchors) if i.strip()]
    
    add_unanchor_contig(cool, tour_db, unanchors)

if __name__ == "__main__":
    main(sys.argv[1:])