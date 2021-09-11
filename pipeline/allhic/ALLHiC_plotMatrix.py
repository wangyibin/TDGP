#!/usr/bin/env python

"""
plot the hic matrix after assembly from ALLHiC

Examples:
    %(prog)s -m sample.cool -a groups.agp
"""


import argparse
import logging

import os
import os.path as op
from re import T, split
import sys
import time
import warnings

from pandas.core.common import SettingWithCopyWarning
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=DeprecationWarning)
warnings.simplefilter(action="ignore", category=SettingWithCopyWarning)
import cooler 
import gc
import numpy as np
import pandas as pd
import tempfile

from collections import OrderedDict
from hicexplorer.reduceMatrix import reduce_matrix
from intervaltree import IntervalTree, Interval 
from itertools import combinations
from scipy.sparse import coo_matrix

logging.getLogger('matplotlib').setLevel(logging.ERROR)
logging.getLogger('cooler').setLevel(logging.ERROR)
logging.getLogger('numexpr').setLevel(logging.ERROR)
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(op.basename(__file__))

def parse_arguments(args=None):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-m', '--matrix', 
            help='cool matrix from ALLHiC_buildMatrix',
            metavar='cool',
            required=True)
    pReq.add_argument('-a', '--agp', 
            help='agp file from ALLHiC.',
            metavar='agp',
            required=True)
    pOpt.add_argument('--binSize', '-bs',
            default=500000, type=int,
            help='Size in bp for the plotting heatmap bins. [default: %(default)s]')
    pOpt.add_argument('--chromSize', default=None,
            help='chromosome size and chromosome order of heatmap, [default: from agp]')
    pOpt.add_argument('--cmap', default='YlOrRd',
            help='colormap of heatmap [default: %(default)s]')
    pOpt.add_argument('-o', '--outprefix', 
            default='prefix_of_cool',
            help='prefix of output heatmap [default: %(default)s]')
    pOpt.add_argument('-om', '--outmatrix', default=None,
            help='output the chromosomal-level matrix [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    return p  

AGP_NAMES_tig = ['chrom', 'start', 'end', 'number',
                 'type', 'id', 'tig_start', 'tig_end', 'orientation']
AGP_NAMES_gap = ['chrom', 'start', 'end', 'number',
                 'type', 'length', 'type', 'linkage', 'evidence']

def import_agp(agpfile, split=True):
    """
    import agp file and return a dataframe
    """
    df = pd.read_csv(agpfile, sep='\t',
                     header=None, index_col=None)
    log.info('load agp file: `{}`'.format(agpfile))
    
    if split:
        tig_df = df[df[4] == 'W']
        gap_df = df[df[4] == 'U']
        tig_df.columns = AGP_NAMES_tig
        gap_df.columns = AGP_NAMES_gap
        tig_df.loc[:, 'chrom'] = tig_df['chrom'].astype('category')
        gap_df.loc[:, 'chrom'] = gap_df['chrom'].astype('category')
        tig_df.set_index('chrom', inplace=True)
        gap_df.set_index('chrom', inplace=True)
        
        return tig_df, gap_df
    else:
        return df

def get_bins(bin_size, chrom_size, start=0, orientation="+", 
                reverse=False, region=None):
    """
    Split the chromosomes into bins of length by specify bin_size

    Params:
    -------
    bin_size: `int` size of bins
    chrom_sizes: `list-like` list of chromosome size
    
    Returns:
    -------
    return a list of tuples 

    Examples:
    -------
    >>> chrom_sizes
    [('Chr1', 29999), ('Chr2', 23333)]
    >>> get_bins(50000, chrom_size, region='contig2’)
    [('Chr2', 0, 23333)]
    """
    bin_intervals = []
    if region:
        tmp = []
        for chroms in chrom_size:
            chrom, size = chroms
            if chrom == region:
                tmp.append(chroms)
        chrom_size = tmp 
    
    for chrom, size in chrom_size:
        length = size - start
        if orientation == "-":
            #print(length, length % bin_size, start, start + bin_size)
            if length % bin_size > 0:
                old_start = start 
                start = start + (length % bin_size)
                bin_intervals.append((chrom, old_start, start))
        for interval in range(start, size, bin_size):
                bin_intervals.append((chrom, interval, 
                                    min(size, interval + bin_size)))
       
    return bin_intervals if not reverse else bin_intervals[::-1]


def get_chrom_index_data(bin_interval):

    chrom_index_data = OrderedDict()
    previous_chrom = ''
    for i, (chrom, _, _) in enumerate(bin_interval):
        if chrom != previous_chrom:
            if previous_chrom:
                chrom_index_data[previous_chrom].append(i - 1)
            chrom_index_data[chrom] = [i]

        previous_chrom = chrom 
    
    return chrom_index_data   


def bedListToIntervalTree(bed_list):
    r"""
    convert a bed list to interval tree

    >>> bed_list = [('Chr1', 0, 10000), ('Chr1', 10000, 20000)]
    >>> res = bedListToIntervalTree(bed_list)
    >>> res['Chr1']
    [Interval(0, 10000), Interval(10000, 20000)]
    """
    bed_interval_tree = {}
    for i, _interval in enumerate(bed_list):
        chrom, start, end = _interval
        if chrom not in bed_interval_tree:
            bed_interval_tree[chrom] = IntervalTree()
        bed_interval_tree[chrom].add(Interval(start, end, i))

    return bed_interval_tree

def getContigIntervalTreeFromAGP(agp_df):
    r"""
    get contig intervaltree from agp
    
    Return:
        {'Chr1': IntervalTree((0, 100, (contig1, "+")))

    """
    contig_bin_interval_tree = {}
    for i, row in agp_df.iterrows():
        
        chrom = row.name
        contig = row.id
        start = row.start
        end = row.end
        orientation = row.orientation
        if chrom not in contig_bin_interval_tree:
            contig_bin_interval_tree[chrom] = IntervalTree()
        
        contig_bin_interval_tree[chrom].add(Interval(start, end, (contig, orientation)))
    
    return contig_bin_interval_tree

def splitAGPByBinSize(agp_df, bin_size=10000):
    r"""
    split chromosome regions into specifial intervals for agp file

    """
    db = []
    for i, row in agp_df.iterrows():
        row.start = row.start - 1
        row.tig_start = int(row.tig_start) - 1
        row.tig_end = int(row.tig_end) + 1
        if int(row.end - row.start + 1) <= bin_size:
            tmp_row = row.copy()
            db.append(tmp_row)
        else:
            tmp_chrom_bins = get_bins(bin_size,  [(row.name, int(row.end))], 
                                start=row.start, orientation=row.orientation)
            tmp_contig_bins = get_bins(
                bin_size, [(row.id, int(row.tig_end))], start=0, 
                reverse=False if row.orientation == "+" else True)
    
            for (_, start, end), (_, tig_start, tig_end) in \
                     zip(tmp_chrom_bins, tmp_contig_bins):
                
                tmp_row = row.copy()
                tmp_row['start', 'end', 'tig_start', 'tig_end'] = start, \
                                    end, tig_start, tig_end
                db.append(tmp_row)
    
    res_df = pd.concat(db, axis=1).T

    return res_df
          

def findIntersectRegion(a, b, fraction=0.51):
    r"""
    find the region by minimum fraction
    >>> a
    IntervalTree([0, 100], [995, 1200])
    >>> b
    (0, 1000)
    >>> findIntersectRegion(a, b)
    [(0, 100)]
    """
    overlaps = sorted(a.overlap(*b))
    tree_begin, tree_end = b[0], b[1]
    
    if not overlaps:
        return overlaps
    start, end = overlaps[0], overlaps[-1]
    if (start.begin - tree_begin) < (start.length() * 0.5):
        overlaps.pop(0)
    if (end.end - tree_end) > (end.length() * 0.5):
        overlaps.pop(-1)
    
    return overlaps

def getContigOnChromBins(chrom_bins, contig_on_chrom_bins):

    _file1 = tempfile.NamedTemporaryFile(delete=False)
    _file2 = tempfile.NamedTemporaryFile(delete=False)
    _file3 = tempfile.NamedTemporaryFile(delete=False)
    
    chrom_bins.to_csv(_file1.name, sep='\t', header=None, index=None)
    contig_on_chrom_bins.to_csv(_file2.name, 
                                    sep='\t', index=True, header=None)

    os.system('bedtools intersect -a {} -b {} -F 0.5 -wo > {} 2>/dev/null'.format(
        _file1.name, _file2.name, _file3.name))

    df = pd.read_csv(_file3.name, sep='\t',
                        header=None, index_col=None)
    _file1.close()
    _file2.close()
    _file3.close()

    df = df.drop([3, 4, 5, 10], axis=1)
    df.columns = ['chrom', 'start', 'end', 'contig',
                  'tig_start', 'tig_end', 'orientation']
    
    return df

def chrRangeID(args, axis=0):
    """
    Chrom range transformation.
    Examples:
    --------
    >>> args = ["Chr1", 100, 200]
    >>> chrRangeID(args)
    "Chr1:100-200"
    >>> args = "Chr1:100-200"
    >>> chrRangeID(args, axis=1)
    ("Chr1", "100", "200")
    """
    if axis == 0:
        chrom, start, end = map(str, args)
        return "{}:{}-{}".format(chrom, start, end)
    elif axis == 1:
        chrom, ranges = args.split(':')
        start, end = ranges.split('-')
        return chrom, start, end
    else:
        return 

def ALLHiC_plotMatrix(args=None):

    args = parse_arguments().parse_args(args)

    start_time = time.time()

    cool = cooler.Cooler(args.matrix)
    agp_df, _ = import_agp(args.agp)
    bin_size = int(cool.binsize)
    ## get chromosome size database from arguments or agp file
    if not args.chromSize:
        chrom_sizes = agp_df.groupby(agp_df.index)['end'].max()
        chrom_sizes = pd.DataFrame(chrom_sizes)
        chrom_sizes.reset_index(inplace=True)
        chrom_sizes.columns = ['chrom', 'length']
    else: 
        chrom_sizes = pd.read_csv(args.chromSize, sep='\t', header=None,
                                index_col=0, names=['chrom', 'length'])
    
    chrom_sizes = [i for _, i in chrom_sizes.iterrows()]
    chrom_bin_interval_df = pd.DataFrame(get_bins(bin_size, chrom_sizes), 
                                    columns=['chrom', 'start', 'end'])
    contig_bins = cool.bins()

    new_agp = splitAGPByBinSize(agp_df, bin_size=bin_size)
    
    split_contig_on_chrom_df = getContigOnChromBins(
        chrom_bin_interval_df, new_agp.drop(['number', 'type'], axis=1))
    split_contig_on_chrom_df.drop(['orientation'], inplace=True, axis=1)
    
    # count_type = matrix[0]['count'].dtype
    # use_type = np.int32 if np.iinfo(np.int32).max >= chrom_bin_length else np.int64
    
    # chrom_matrix = coo_matrix((chrom_bin_length, chrom_bin_length), dtype=count_type)
    # row = np.zeros(cool.info['nnz'], dtype=use_type)
    # col = np.zeros(cool.info['nnz'], dtype=use_type)
    # data = np.zeros(cool.info['nnz'], dtype=count_type)
    

    contig_bins_index = contig_bins[:].apply(lambda x: chrRangeID(x.values), axis=1)
    contig_bins_index = contig_bins_index.to_frame().reset_index()
    contig_bins_index.rename(columns={'index': 'contigidx', 0: 'contig_region'}, inplace=True)
    contig_bins_index.set_index('contig_region', inplace=True)

    chrom_bins_index = chrom_bin_interval_df.apply(lambda x: chrRangeID(x.values), axis=1)
    chrom_bins_index = chrom_bins_index.to_frame().reset_index()
    chrom_bins_index.rename(columns={'index': 'chromidx', 0: 'chrom_region'}, inplace=True)
    chrom_bins_index.set_index('chrom_region', inplace=True)

    func = lambda row: chrRangeID(row.values)
    split_contig_on_chrom_df['chrom_region'] = split_contig_on_chrom_df[[
                                                'chrom', 'start', 'end']].apply(func, axis=1)
    split_contig_on_chrom_df['contig_region'] = split_contig_on_chrom_df[['contig',
                                                'tig_start', 'tig_end']].apply(func, axis=1)

    split_contig_on_chrom_df['chromidx'] = chrom_bins_index.loc[
                                               split_contig_on_chrom_df['chrom_region']]['chromidx'].values
    split_contig_on_chrom_df['contigidx'] = contig_bins_index.loc[
                                                split_contig_on_chrom_df['contig_region']]['contigidx'].values


    split_contig_on_chrom_df.drop(['chrom', 'start', 'end',
                                    'contig', 'tig_start', 'tig_end',
                                    'chrom_region', 'contig_region'], 
                                    inplace=True, axis=1)
 
    log.info('starting')

    matrix = cool.matrix(balance=False, sparse=True)
    contig_matrix = matrix[:].todense()

    contig_idx = split_contig_on_chrom_df.contigidx
    # order the contig-level matrix to chromosome-level
    ordered_contig_matrix = contig_matrix[contig_idx, :][:, contig_idx]
    
    # reindex the contigidx to new index
    split_contig_on_chrom_df.contigidx = split_contig_on_chrom_df.index
    grouped_contig_idx = split_contig_on_chrom_df.groupby('chromidx')[
                                        'contigidx'].aggregate(list)
    grouped_contig_border = grouped_contig_idx.apply(
                                        lambda x: x if len(x) == 1 else [x[0]]).values
    grouped_contig_border = [i for sublist in grouped_contig_border 
                                                            for i in sublist]

    chrom_matrix = np.add.reduceat(np.add.reduceat(ordered_contig_matrix, 
                                            grouped_contig_border,  axis=0), 
                                            grouped_contig_border, axis=1)


    ## slowest method
    """
    i = 0
    for contig1, contig2 in combinations(cool.chromnames, 2):
        tmp_contig_pixels = matrix.fetch(contig1, contig2)
        if tmp_contig_pixels.empty:
            continue
        bin1_id = tmp_contig_pixels.bin1_id.values
        bin2_id = tmp_contig_pixels.bin2_id.values
        counts = tmp_contig_pixels['count']
        try:
            chrom_idx1 = split_contig_on_chrom_df.loc[bin1_id].chromidx.values
            chrom_idx2 = split_contig_on_chrom_df.loc[bin2_id].chromidx.values
        except KeyError:
            continue
        
        tmp_chrom_idx1 = chrom_idx1.copy()
        tmp_chrom_idx2 = chrom_idx2.copy()
        expression = tmp_chrom_idx1 > tmp_chrom_idx2
        chrom_idx1[expression] = tmp_chrom_idx2[expression]
        chrom_idx2[expression] = tmp_chrom_idx1[expression]
        size = len(tmp_contig_pixels)
        row[i:i + size] = chrom_idx1
        col[i: i + size] = chrom_idx2
        data[i: i + size] = counts
        i += size

    del tmp_chrom_idx1, tmp_chrom_idx2, bin1_id, bin2_id, counts, tmp_contig_pixels
    """

    chrom_matrix = coo_matrix(np.triu(chrom_matrix))
    chrom_pixels = pd.Series.sparse.from_coo(chrom_matrix)
    chrom_pixels = pd.DataFrame(chrom_pixels).reset_index()
    chrom_pixels.columns = ['bin1_id', 'bin2_id', 'count']
 
    hic_metadata = {}
    hic_metadata['matrix-generated-by'] = np.string_(
        'ALLHiC_plotMatrix'
    )
    hic_metadata['matrix-generated-by-url'] = np.string_(
        'https://github.com/tangerzhang/ALLHiC'
    )
    cooler.create_cooler(args.matrix.rsplit(".", 1)[
                         0] + ".chrom.cool", chrom_bin_interval_df, 
                         chrom_pixels, metadata=hic_metadata)
    
    log.info('Done, elasped time {} s'.format(time.time() - start_time))


if __name__ == "__main__":
    #ALLHiC_plotMatrix(sys.argv[1:])
    os.chdir('/share/home/stu_wangyibin/test/ALLHiC_plot')
    ALLHiC_plotMatrix(['--agp', 'chrn.V2.agp', '--matrix', 'JGY.cool'])