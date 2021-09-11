#!/usr/bin/env python

"""
build hic matrix from bwa mem for ALLHiC plot.

Examples:
    %(prog)s -b sample.bwa_mem.bam -o sample.cool
"""

import argparse
import logging
import os
import os.path as op
import sys
import time

import warnings
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=DeprecationWarning)

import cooler 
import numpy as np
import pysam 
import pandas as pd

from collections import OrderedDict
from ctypes import c_uint
from multiprocessing import Process, Queue
from multiprocessing.sharedctypes import RawArray 
from scipy.sparse import coo_matrix

logging.getLogger('matplotlib').setLevel(logging.ERROR)
logging.getLogger('cooler').setLevel(logging.ERROR)
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(op.basename(__file__))

def parse_arguments(args=None):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-b', '--bamFile', 
            help='The alignment file from reads mapping',
            required=True)
    pReq.add_argument('-o', '--output', 
            help='Output hicmatrix as cool format',
            metavar='cool',
            required=True)
    pOpt.add_argument('--binSize', '-bs',
            default=100000, type=int,
            help='Size in bp for the bins. [default: %(default)s]')
    pOpt.add_argument('--mapq', default=None, type=int,
            help='map quality of alignments [default: %(default)s]')
    pOpt.add_argument('-t', '--threads', type=int, default=8,
            help='number of program threads[default:%(default)s]')
    pOpt.add_argument('--bufferSize', default=400000, 
            type=int,
            help='Size of the buffer of each thread. '
                'Reduce it to decrease memory.'
                ' [default: %(default)s]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    return p  

def getContigSizes(bamFileHandle):
    """
    get contig size from the bam header, and return size list
    [('contig1', 29999), ('contig2', 233333)]
    >>> getContigSizes(pysam.Samfile('sample.bam', 'rb')
    [('contig', 29999), ('contig2', 233333)]
    """
    contig_sizes = OrderedDict(zip(bamFileHandle.references,
                                        bamFileHandle.lengths))
    
    return list(contig_sizes.items())

def get_bins(bin_size, chrom_size, region=None):
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
    [('contig1', 29999), ('contig2', 23333)]
    >>> get_bins(50000, chrom_size, region='contig2’)
    [('contig2', 0, 23333)]
    """
    bin_intervals = []
    start = 0
    if region:
        tmp = []
        for chroms in chrom_size:
            chrom, size = chroms
            if chrom == region:
                tmp.append(chroms)
        chrom_size = tmp 
    
    for chrom, size in chrom_size:
        for interval in range(start, size, bin_size):
            bin_intervals.append((chrom, interval, 
                                min(size, interval + bin_size)))
    
    return bin_intervals


def get_ctg_index_data(bin_interval):

    ctg_index_data = OrderedDict()
    previous_ctg = ''
    for i, (ctg, _, _) in enumerate(bin_interval):
        if ctg != previous_ctg:
            if previous_ctg:
                ctg_index_data[previous_ctg].append(i - 1)
            ctg_index_data[ctg] = [i]

        previous_ctg = ctg 
    
    return ctg_index_data   


def parseBamFiles(bamFileHandle, buffer, mapq):
    """
    parse bam file into severals buffers
    """
    buffer_data = []
    all_read_processed = False
    iter_num = 0
    i = 0
    while i < buffer:
        try:
            read = next(bamFileHandle)
        except StopIteration:
            all_read_processed = True 
            break
        
        buffer_data.append(read)
        i += 1 
    
    if all_read_processed and len(buffer_data) != 0:
        return buffer_data, True, iter_num - len(buffer_data)
    if all_read_processed and len(buffer_data) == 0:
        return None, True

    return buffer_data, False, iter_num - len(buffer_data)

def process_data(buffer_data, queue_out,
                    bin_size, contig_index_db,
                    mapq, row, col, data):
    """
    process read records to a matrix

    """

    pair_added = 0

    if buffer_data is None:
        queue_out.put([None])
        return 

    for i, read in enumerate(buffer_data):

        if read.is_unmapped or read.mate_is_unmapped:
            continue

        if read.is_read2:
            continue

        if mapq:
            if read.mapq < mapq:
                continue
        
        ctg1 = read.reference_name
        read_pos1 = read.reference_start
        ctg2 = read.next_reference_name
        read_pos2 = read.next_reference_start
        
        contig_start_index1 = contig_index_db[ctg1][0]
        contig_start_index2 = contig_index_db[ctg2][0]
        
        contig_index1 = contig_start_index1 + read_pos1 // bin_size 
        contig_index2 = contig_start_index2 + read_pos2 // bin_size

        if contig_index1 > contig_index2:
            contig_index1, contig_index2 = contig_index2, contig_index1
        row[i] = contig_index1
        col[i] = contig_index2
        data[i] = np.uint8(1)

        pair_added += 1
    
    queue_out.put([[pair_added]])

    return
   
def ALLHiC_buildMatrix(args=None):
    
    args = parse_arguments().parse_args(args)
    
    if args.threads < 2:
        args.threads = 2 
        log.warn("At least 2 threads, setting --threads=2 !")
    
    bin_size = args.binSize

    threads = args.threads - 1
    buffer_size = args.bufferSize
    mapq = args.mapq 
    
    start_time = time.time()
    start_times = [0] * threads 
    log.info("import bam file `{}` to build matrix".format(args.bamFile))
    bamHandle = pysam.Samfile(args.bamFile, 'rb', threads=threads)
    contig_sizes = getContigSizes(bamHandle)
    bin_intervals = get_bins(bin_size, contig_sizes)
    matrix_size = len(bin_intervals)
    contig_index_db = get_ctg_index_data(bin_intervals)

    hic_matrix = coo_matrix((matrix_size, matrix_size), dtype='uint32')

    row = []
    col = []
    data = []
    for i in range(threads):
        row.append(RawArray(c_uint, buffer_size))
        col.append(RawArray(c_uint, buffer_size))
        data.append(RawArray(c_uint, buffer_size))
    
    iter_num = 0
    buffer_workers = [None] * threads
    process = [None] * threads
    all_data_processed = False
    
    queue = [None] * threads
    thread_done = [False] * threads 
    all_threads_done = False 

    while not all_data_processed or not all_threads_done:
        for i in range(threads):
            if queue[i] is None and not all_data_processed:
                start_times[i] = time.time()
                buffer_workers[i], all_data_processed, \
                    iter_num_ = parseBamFiles(bamHandle, 
                                              args.bufferSize, 
                                              args.mapq)

                iter_num += iter_num_

                queue[i] = Queue()
                thread_done[i] = False 

                process[i] = Process(target=process_data, 
                                args=(buffer_workers[i],
                                        queue[i],
                                        bin_size,
                                        contig_index_db,
                                        mapq,
                                         row[i], 
                                         col[i], 
                                         data[i]))
                
                process[i].start()
            
            elif queue[i] is not None and not queue[i].empty():
                result = queue[i].get()
                
                if result[0] is not None:
                    pair_added = result[0][0]
                    if hic_matrix is None:
                        hic_matrix = coo_matrix(
                            (data[i], (row[i], col[i])), 
                            shape=(matrix_size, matrix_size))
                    else:
                        hic_matrix += coo_matrix(
                            (data[i], (row[i], col[i])), 
                            shape=(matrix_size, matrix_size))

                ## clear done thread
                buffer_workers[i] = None
                queue[i] = None 
                process[i].join()
                process[i].terminate()
                process[i] = None
                ## raise this thread done signal
                thread_done[i] = True

                elapsed_time = time.time() - start_times[i]
                log.info("{} pairs processed took {:.2f} s".format(buffer_size, elapsed_time))

            elif all_data_processed and queue[i] is None:
                thread_done[i] = True
            else:
                time.sleep(1)
        
        ## whether done all thread works
        if all_data_processed:
            all_threads_done = True
            for thread in thread_done:
                if not thread:
                    all_threads_done = False 
    
    bins = pd.DataFrame(bin_intervals, columns=('chrom', 'start', 'end'))
    pixels = pd.Series.sparse.from_coo(hic_matrix.tocoo())
    pixels = pd.DataFrame(pixels).reset_index()
    pixels.columns = ['bin1_id', 'bin2_id', 'count']

    hic_metadata = {}
    hic_metadata['matrix-generated-by'] = np.string_(
        'ALLHiC_buildMatrix'
    )
    hic_metadata['matrix-generated-by-url'] = np.string_(
        'https://github.com/tangerzhang/ALLHiC'
    )
    cooler.create_cooler(args.output, bins, pixels, metadata=hic_metadata)
    log.info('Successful output hic matrix in `{}`'.format(args.output))
    log.info('Done, elapsed time {:.2f} s'.format(time.time() - start_time))

if __name__ == "__main__":
    ALLHiC_buildMatrix(sys.argv[1:])