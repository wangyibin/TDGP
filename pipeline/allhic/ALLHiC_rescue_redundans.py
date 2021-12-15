#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
ALLHiC rescue intergating with redundans
"""

import argparse
import logging
import os
import os.path as op
from re import sub
import sys
import subprocess 
import pysam
import re
import time
from FastaIndex import FastaIndex

from multiprocessing import Pool, Process
# from pyfaidx import Fasta, Faidx

def time_print(info, type='info'):
    if type != 'info':
        info = "\033[35m%s\033[0m"%info
    print("\033[32m%s\033[0m %s"%(time.strftime('[%H:%M:%S]',time.localtime(time.time())), info))


## redundans from https://github.com/lpryszcz

def run_last(anchor_fasta, unplaced_fasta, identity, threads, verbose=1):
    """Start LAST with multi-threads"""
    if verbose:
        sys.stderr.write(" Running LAST...\n")
    # build db
    
    if not os.path.isfile(anchor_fasta+".suf"):
        os.system("lastdb -P %s -W 11 %s %s" % (threads, anchor_fasta, anchor_fasta))
    # run LAST
    args1 = ["lastal", "-P", str(threads), "-f", "TAB", anchor_fasta, unplaced_fasta]
    proc1 = subprocess.Popen(args1, stdout=subprocess.PIPE, stderr=sys.stderr)
    return proc1

def _qhits_generator(handle, minLength):
    pq, pqsize, hits = '', 0, {}
    for l in handle: 
        
        if not l.strip():
            continue
        if l[0] == '#': 
            continue
    
        # unpack
        (score, t, tstart, talg, tstrand, tsize, q, qstart, qalg, qstrand, qsize, blocks) = l.split()[:12]
        (score, qstart, qalg, qsize, tstart, talg, tsize) = map(int, (score, qstart, qalg, qsize, tstart, talg, tsize))
        # skip reverse matches

        if t==q or tsize<qsize or qsize<minLength or tsize==qsize and t<q: 
            continue
        
        # report previous query
        if pq != q:
            
            yield pq, pqsize, hits
            # reset
            pq, pqsize, hits = q, qsize, {}
        if t not in hits:
            hits[t] = []
        if qstrand=="+":
            s = qstart
            e = s + qalg
        else:
            e = qsize - qstart
            s = qsize - qstart - qalg
        
        hits[t].append((score, t, s, e))
    # make sure to report last bit
    if hits:
        
        yield pq, pqsize, hits

def _overlap(s, e, score, hits, maxfrac=0.1):
    """Return True if at least 10% overlap with existing hit"""
    maxoverlap = maxfrac*(e-s)
    selection = lambda x: x[0]<s<x[1] and s+maxoverlap<x[1] or x[0]<e<x[1] and e-maxoverlap>x[0] or \
                          s<x[0]<e and x[0]+maxoverlap<e or s<x[1]<e and x[1]-maxoverlap>s
    if filter(selection, hits):
        return True
            
def hits2valid(hits, q, qsize, identityTh, overlapTh):
    """Return valid matches for particular query"""
    for t, (score, qalg, se) in hits.items():
        identity = 1.0 * (score+(qalg-score)/2) / qalg
        if qalg>qsize: qalg = qsize
        overlap  = 1.0 * qalg / qsize
        # filter by identity and overlap
        if identity >= identityTh and overlap >= overlapTh:
            yield score, t, q, qalg, identity, overlap

def fasta2hits(anchor_fasta, unplaced_fasta, threads, identityTh, overlapTh, minLength):
    """Return LASTal hits passing identity and overlap thresholds"""
    # execute last
    last = run_last(anchor_fasta, unplaced_fasta, identityTh, threads)
    
    for q, qsize, qhits in _qhits_generator(last.stdout, minLength):
        hits = {}
        
        for t in qhits:
            hits[t] = [0, 0, []]
            for score, t, s, e in sorted(qhits[t], reverse=1):
                if _overlap(s, e, score, hits[t][2]):
               
                    continue
                
                hits[t][0] += score
                hits[t][1] += e-s
                hits[t][2].append((s, e, score))
        for d in hits2valid(hits, q, qsize, identityTh, overlapTh):
            yield d
               
def fasta2skip(anchor_fasta, unplaced_fasta, faidx, threads, identityTh, overlapTh, minLength):
    """Return dictionary with redundant contigs and their best alignments"""
    # get hits generator
    hits = fasta2hits(anchor_fasta, unplaced_fasta, threads, identityTh, overlapTh, minLength)
    # iterate through hits
    identities, sizes = [], []
    contig2skip = {c: 0 for c in faidx} 
    
    for i, (score, t, q, algLen, identity, overlap) in enumerate(hits, 1):
        # store first match or update best match
    
        if q not in contig2skip or not contig2skip[q] or score > contig2skip[q][0]:
            contig2skip[q] = (score, t, algLen, identity, overlap)
        # store identity and alignment for plotting
        identities.append(identity)
        sizes.append(algLen)
    
    return contig2skip


def create_qry_file(source_cds, gff, target_cds, target_bed):
    src_cds_db = read_fasta(source_cds)
    idx = 1
    qry_db = {}
    with open(target_cds, 'w') as fcds:
        with open(target_bed, 'w') as fbed:
            with open(gff, 'r') as fin:
                for line in fin:
                    if line.strip() == '' or line[0] == '#':
                        continue
                    data = line.strip().split()
                    if data[2] != 'gene':
                        continue
                    id = re.findall(r'Name=(.*)', data[8])[0].split(';')[0]
                    new_id = "%s_%d"%(id, idx)
                    idx += 1
                    chrn = data[0]
                    sp = int(data[3])
                    ep = int(data[4])
                    if ep <= sp:
                        continue
                    direct = data[6]

                    if chrn not in qry_db:
                        qry_db[chrn] = set()
                    qry_db[chrn].add(new_id)
                    try:
                        fcds.write(">%s\n%s\n"%(new_id, src_cds_db[id]))
                        fbed.write("%s\t%d\t%d\t%s\t0\t%s\n"%(chrn, sp, ep, new_id, direct))
                    except KeyError:
                        continue
    return qry_db


def read_anchors(anchors_file):
    anchor_db = {}
    with open(anchors_file, 'r') as fin:
        for line in fin:
            if line.strip() == '' or line[0] == '#':
                continue
            data = line.strip().split()
            qry_gn = data[0]
            ref_gn = data[1]
            anchor_db[qry_gn] = ref_gn
    return anchor_db


def convert_query_db(qry_db, anchor_db):
    new_qry_db = {}
    for chrn in qry_db:
        new_qry_db[chrn] = set()
        for gn in qry_db[chrn]:
            if gn in anchor_db:
                new_qry_db[chrn].add(anchor_db[gn])
    return new_qry_db



def get_black_list(args):
    chrn, anchor_fasta, unplaced_fasta, threads, \
                identity, overlaps, minLength = args

    faidx = FastaIndex(unplaced_fasta)
  

    contig2skips = fasta2skip(anchor_fasta, unplaced_fasta,
                        faidx, threads, identity, overlaps,
                        minLength)
    contig2skips =  list(filter(lambda x: x[1] !=0, contig2skips.items()))
    black_list = list(map(lambda x: x[0], contig2skips))
    
    return chrn, black_list

def get_black_db(anchor_fastas, unplaced_fasta, threads=4,
                identity=0.51, overlaps=0.8, minLength=200,
                ):
    raw_black_db = {}
    pool = Pool(10)
    args = []
    for chrn, anchor_fasta in anchor_fastas.items():
        args.append((chrn, anchor_fasta, unplaced_fasta,
                            threads, identity, overlaps, minLength))
    
    res = []
    for arg in args: 
        r = pool.apply_async(get_black_list, (arg,))
        res.append(r)
    results = [ result.get() for result in res]
    
    
    pool.close()
    pool.join()
    raw_black_db = dict(res)
    with open('black.list', 'w') as out:
        for chrn in sorted(raw_black_db):
            out.write('{}\t{}\n'.format(chrn, ' '.join(raw_black_db[chrn])))
    # black_db = {}
    # for chrn in raw_black_db:
    #     if chrn not in black_db:
    #         black_db[chrn] = set()
    #     black_db[chrn].add(chrn)
    
    return raw_black_db    

def get_clusters(clu):
    clu_db = {}
    clu_ctgs = {}
    with open(clu, 'r') as fin:
        for line in fin:
            if line[0] == '#':
                continue
            data = line.strip().split()
            chrn = data[0]
            ctgs = data[2:]
            clu_db[chrn] = ctgs
            for ctg in ctgs:
                clu_ctgs[ctg] = chrn
    return clu_db, clu_ctgs


def get_ovlp(qry_set1, qry_set2):
    return len(qry_set1.intersection(qry_set2))

def get_hic_signal(bam):
    signals = {}
    with pysam.AlignmentFile(bam, 'rb') as fin:
        for line in fin:
            ctg1 = line.reference_name
            pos1 = line.reference_start
            ctg2 = line.next_reference_name
            pos2 = line.next_reference_start
            if pos1==-1 or pos2==-1:
                continue
            if ctg1 not in signals:
                signals[ctg1] = {}
            if ctg2 not in signals[ctg1]:
                signals[ctg1][ctg2] = 0
            signals[ctg1][ctg2] += 1

            if ctg2 not in signals:
                signals[ctg2] = {}
            if ctg1 not in signals[ctg2]:
                signals[ctg2][ctg1] = 0
            signals[ctg2][ctg1] += 1
    return signals

def get_counts(counts):
    header = ""
    counts_db = {}
    with open(counts, 'r') as fin:
        for line in fin:
            if line[0] == '#':
                header = line
            else:
                ctg = line.strip().split()[0]
                counts_db[ctg] = line
    
    return header, counts_db

def read_fasta(in_fa):
    fa_db = {}
    with open(in_fa, 'r') as fin:
        for line in fin:
            if line[0] == '>':
                id = line.strip().split()[0][1:]
                fa_db[id] = []
            else:
                fa_db[id].append(line.strip())
    for id in fa_db:
        fa_db[id] = ''.join(fa_db[id])
    
    return fa_db

def ALLHiC_rescue(ref, bam, black_db, exclude, clu, counts, gff3, jprex, wrk):
    if not os.path.exists(wrk):
        os.mkdir(wrk)
    
    ref = os.path.abspath(ref)
    bam = os.path.abspath(bam)
    clu = os.path.abspath(clu)
    counts = os.path.abspath(counts)
    bed = os.path.abspath(jprex+'.bed')
    cds = os.path.abspath(jprex+'.cds')
    gff3 = os.path.abspath(gff3)

    jprex = jprex.split('/')[-1]
    time_print("Entering: %s"%wrk)
    os.chdir(wrk)

    os.system("ln -sf %s %s.cds"%(cds, jprex))
    os.system("ln -sf %s %s.bed"%(bed, jprex))
    new_cds = "dup.cds"
    new_bed = "dup.bed"

    qry_db = create_qry_file(cds, gff3, new_cds, new_bed)
    
    if not os.path.exists("dup.%s.anchors"%jprex):
        time_print("Running jcvi", type="important")
        cmd = "python -m jcvi.compara.catalog ortholog dup %s &> jcvi.log"%jprex
        os.system(cmd)
    else:
        time_print("Anchors file found, skip", type="important")
    
    time_print("Loading anchors file")
    anchor_db = read_anchors("dup.%s.anchors"%jprex)
    
    time_print("Converting query db")
    qry_db = convert_query_db(qry_db, anchor_db)

    time_print("Loading clusters")
    clu_db, clu_ctgs = get_clusters(clu)
    clu_set = {}
    for chrn in clu_db:
        clu_set[chrn] = set()
        for ctg in clu_db[chrn]:
            if ctg not in qry_db:
                continue
            clu_set[chrn] = clu_set[chrn].union(qry_db[ctg])
    
    
    remain_ctgs = []
    ctg_db = read_fasta(ref)

    for ctg in ctg_db:
        if ctg not in clu_ctgs:
            remain_ctgs.append([ctg, len(ctg_db[ctg])])
    
    time_print("Loading HiC signals")
    signal_db = get_hic_signal(bam)

    time_print("Get best matches")
    for ctg, ctgl in sorted(remain_ctgs, key=lambda x: x[1], reverse=True):
        score_list = []
        if ctg not in signal_db:
            continue
        for ctg2 in signal_db[ctg]:
            if ctg2 not in clu_ctgs:
                continue
            sig = signal_db[ctg][ctg2]
            ovlp = 0
            chrn = clu_ctgs[ctg2]
            if ctg in qry_db:
                ovlp = get_ovlp(qry_db[ctg], clu_set[chrn])
            score_list.append([ovlp, sig, chrn])

        if len(score_list)==0:
            continue
        best_match = None
        i = 0
        while not best_match and i < len(score_list):
            best_match = sorted(score_list, key=lambda x: [x[0], -x[1]])[i]
            if best_match[2] in exclude:
                best_match = ""
                i += 1 
                continue
            if ctg in black_db[best_match[2]]:
                best_match = ""
            i += 1
        if not best_match:
            continue
        if best_match[1] < 10:
            continue
        time_print("\t%s matched %s, sig: %d, ovlp: %d"%(ctg, best_match[2], best_match[1], best_match[0]))
        clu_db[best_match[2]].append(ctg)
        if ctg in qry_db:
            clu_set[best_match[2]] = clu_set[best_match[2]].union(qry_db[ctg])
    
    time_print("Writing new groups")
    header, counts_db = get_counts(counts)
    for chrn in clu_db:
        with open("%s.txt"%chrn, 'w') as fout:
            fout.write(header)
            for ctg in clu_db[chrn]:
                fout.write(counts_db[ctg])     
        
    os.chdir("..")
    time_print("Finished")

def main(args):
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('-r', '--ref', help="Contig level assembly fasta", required=True)
    pReq.add_argument('-b', '--bam', help="Unprunned bam", required=True)
    pReq.add_argument('-c', '--cluster', help="Cluster file of contigs", required=True)
    pReq.add_argument('-n', '--counts', help="count REs file", required=True)
    pReq.add_argument('-g', '--gff3', help="Gff3 file generated by gmap cds to contigs", required=True)
    pReq.add_argument('-j', '--jcvi', help="CDS file for jcvi, bed file with same prefix must exist in the same position", required=True)
    pReq.add_argument('-w', '--workdir', help="Work directory, default=wrkdir", default="wrkdir")
    pReq.add_argument('-e', '--exclude', help="exclude chromosome, comma seperated", default=None)
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    ref = args.ref
    bam = args.bam
    clu = args.cluster
    counts = args.counts
    gff3 = args.gff3
    jprex = args.jcvi
    exclude = args.exclude.split(',') if args.exclude else None
    jprex = '.'.join(jprex.split('.')[:-1])
    wrk = args.workdir

    anchor_fastas = {}
    anchor_ctgs = []
    with open(clu) as fp:
        for line in fp:
            if line.strip() and line[0] == '#':
                continue
            line_list = line.strip().split("\t")
            chrn, num, contigs = line_list 
            contigs = contigs.strip().split()
            anchor_ctgs.extend(contigs)
            if chrn not in anchor_fastas:
                anchor_fastas[chrn] = chrn + ".contig.fasta"
            
            with open(chrn + ".contig.list", 'w') as out:
                out.write('\n'.join(contigs))
            
    command = 'seqkit grep -f {0}.contig.list {1} -o {0}.contig.fasta'
    processes = set()
    max_processes = 10 
    for chrn in anchor_fastas:
        processes.add(subprocess.Popen(command.format(chrn, ref), shell=True))
        if len(processes) >= max_processes:
            processes.difference_update([
                    p for p in processes if p.poll() is not None 
            ])

    with open('all.anchor.contig.list', 'w') as out:
        out.write('\n'.join(anchor_ctgs))
    os.system('seqkit grep -v -f all.anchor.contig.list {} > unplaced.fa'.format(ref))
            
    unplaced_fasta = 'unplaced.fa'
    if not op.exists('black.list'):        
        black_db = get_black_db(anchor_fastas, unplaced_fasta)
    else:
        black_db = {i.strip().split("\t")[0]: i.strip().split("\t")[1].split() 
                        for i in open('black.list')}
    
    ALLHiC_rescue(ref, bam, black_db, exclude, clu, counts, gff3, jprex, wrk)
    
if __name__ == "__main__":
    main(sys.argv[1:])