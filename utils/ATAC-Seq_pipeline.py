#!/usr/bin/env python
# -*- coding:utf-8 -*-


from __future__ import print_function

import glob
import multiprocessing as mp
import os
import os.path as op
import sys


def build_index(ref="reference.fasta", threads=12):
    if not op.exists(ref + ".1.bt2"):
        os.system("bowtie2-build --thread {0} {1} {1}".format(threads, ref))
    else:
        print("[bowtie2-build]: The index is allready exists. skip...", file=sys.stderr)

def is_fastq(fqfile):
    fastq_suffixes = ['fastq', 'fq']
    if fqfile[-2:] == "gz":
        fq_suffix = fqfile[:-3]
        if fq_suffix.split('.')[-1] in fastq_suffixes:
            return True
        else:
            return False
    elif fqfile.split('.')[-1] in fastq_suffixes:
        return True
    else:
        return False


def find_fastq(data_dir):
    if not op.isdir(data_dir):
        sys.exit("There is no directory of {}".format(data_dir))
    
    file_list = os.listdir(data_dir)
    file_list = list(map(lambda x: "{}/{}".format(data_dir, x), file_list))
    file_list = filter(op.exists, file_list)
    file_list = filter(is_fastq, file_list)
    

    return file_list


def glob_fastq(data_dir, fq_suffix="_R1.fastq.gz"):
    fq_1_list = glob.glob('{}/*{}'.format(data_dir, fq_suffix))
    return fq_1_list


def rmdup(name):
    cmd = "java -jar $picard MarkDuplicates INPUT={0}.sorted.bam OUTPUT={0}.rmdup.bam METRICS_FILE={0}.metrics VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000".format(name)
    os.system(cmd)


def bam_coverage(name):
    cmd = "samtools index {0}.rmdup.bam && bamCoverage -b {0}.rmdup.bam -o {0}.bw".format(name)
    os.system(cmd)


def macs(name):
    genome_size = os.popen(r"""samtools view -H {}.rmdup.bam | grep "^@SQ" |awk '{{split($2, SN,":");split($3,LN,":");print SN[2]"\t"LN[2]}}' | awk -F "\t" '{{num+=$2}}END{{print num}}'""".format(name), 'r').next().strip()
    cmd = "macs2 callpeak -t {0}.rmdup.bam -f BAM -n {0} -g {1} --nomodel --extsize 70 -q 0.01".format(name, genome_size)
    os.system(cmd)


def align_by_bowtie2(fq1, ref, fq_suffix="_R1.fastq.gz", is_paired=True, 
        is_coverage=True, is_macs=True, threads=12):
    cmd = "bowtie2 -p {0} -N 1 -X 2000 -q -x {1} {2} 2>{3}.stat | samtools sort -@ {0}  | samtools view -@ 12 -bS > {3}.sorted.bam"
    if is_paired:
        fq2 = fq1.replace('1.', '2.')
        if not op.exists(fq2):
            return None
        cmd_mid = "-1 {} -2 {}".format(fq1, fq2)
    else:
        cmd_mid = "-U {}".format(fq1)
    name = op.basename(fq1).replace(fq_suffix, "")
    cmd = cmd.format(threads, ref, cmd_mid, name )
    if not op.exists(name + ".sorted.bam"):
        os.system(cmd)
    else:
        print("[alignment]: The {}.sorted.bam existed. skip...".format(name), file=sys.stderr)
    
    if not op.exists(name + ".rmdup.bam"):
        rmdup(name)
    else:
        print("[rmdup]: The {}.rmdup.bam existed. skip...".format(name), file=sys.stderr)

    if is_coverage:
        if not op.exists(name + ".bw"):
            bam_coverage(name)
        else:
            print("[bamCoverage]: The {}.bw existed. skip...".format(name), file=sys.stderr)

    if is_macs:
        if not op.exists(name + "_peaks.narrowPeak"):
            macs(name)
        else:
            print("[stringtie]: The {}_peaks.narrowPeak existed. skip...".format(name), file=sys.stderr)


    
def main(data_dir, ref, fq_suffix="_R1.fastq.gz", is_paired=True, 
        is_coverage=True, is_macs=True, threads=12):
    build_index(ref, threads)
    fq_list = glob_fastq(data_dir, fq_suffix)
    
    for fq in fq_list:
        align_by_bowtie2(fq, ref, fq_suffix, is_paired, is_coverage, is_macs, threads)



if __name__ == "__main__":
    from optparse import OptionParser

    p = OptionParser(__doc__)


    p.add_option("--fq_suffix", default='_R1.fastq.gz', 
            help="The suffix of the fastq file, if paired-end data, "
            "should give R1. [default: %default]")
    p.add_option("--is_paired", default=True, action="store_false",
            help="is single data [default: %default]")
    p.add_option("--is_coverage", default=True, action="store_false",
            help="is calculate the coverage of bam file [defaule: %default]")
    p.add_option("--is_macs", default=True, action="store_false",
            help="is to run macs [default: %default]")
    p.add_option("-t", "--threads", default=12, type=int,
            help="the thread number of program [default: %default]")
    
    opts, args = p.parse_args()

    if len(args) != 2:
        sys.exit(p.print_help())

    data_dir, ref = args

    main(data_dir, ref, opts.fq_suffix, opts.is_paired, opts.is_coverage, opts.is_macs, opts.threads)
