#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
%(prog)s <vcf> <gff> [Options]

    Calculate heterozygosity for per gene from a snp vcf file.
    output: `chrom start end gene snp_num gene_length heter_rate`

    ** Total heterozygosity is output in stdout.
    Examples:
        %(prog)s sample.vcf sample.gff -o out.tsv 
"""

from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys

import multiprocessing
import pandas as pd 
import numpy as np 

logging.basicConfig(level=logging.DEBUG)
def filter_vcf(args):
    vcf_df, _type, minDP, remove_snp, remove_indel, out = args
    #allow_type = ('0/1', '1/1', '0/0', 'heter', 'homo', 'variant', 'all')
    
    vcf_df.dropna(inplace=True)
    vcf_df = vcf_df[vcf_df.SAMPLE.map(lambda x: x.split(":")[0]).isin(_type)]
    if remove_snp:
        vcf_df = vcf_df[(vcf_df['REF'].map(lambda x: len(x)) > 1) | \
                        (vcf_df['ALT'].map(lambda x: len(x)) > 1)]
    if remove_indel:
        vcf_df = vcf_df[(vcf_df['REF'].map(lambda x: len(x)) == 1) & \
                        (vcf_df['ALT'].map(lambda x: len(x)) == 1)]
    def deal_info(info):
        return dict(map(lambda x: x.split("="), info.split(";")))
    if minDP:
        vcf_df = vcf_df[vcf_df.INFO.map(lambda x: int(deal_info(x)['DP']) >= minDP)]
    
    if out:
        try:
            vcf_df.to_csv(out, sep='\t', header=None)
        except:
            logging.warning('Warning: Failed to output filter vcf result...')
      
    return vcf_df


def import_vcf(vcf_path,
               _type='0/1',
               minDP=0,
               chunksize=1e6,
               threads=12, 
               remove_snp=False, 
               remove_indel=True,
               out=None):
    standard_type = ('0/1', '1/1', '0/0')
    standard_type_dict = {'heter': '0/1', 'homo': '1/1'}
    if _type not in standard_type:
        if _type == 'all':
            _type = ['0/1', '1/1', '0/0', './.']
        elif _type == 'variant':
            _type = ['0/1', '1/1']
        elif _type == 'heter':
            _type = ['0/1']
        elif __type == 'homo':
            _type = ['1/1']
        else:
            logging.error('Error input the vcf type `{}`'.format(_type))
    else:
        _type = [_type] 
    
    compression = 'gzip' if vcf_path.endswith('.gz') else 'infer'
    vcf_df = pd.read_csv(vcf_path, sep='\t', header=None, 
                         index_col=0,comment="#",
                         chunksize=chunksize,
                         compression=compression,
                         names=('CHROM', 'POS', 
                              'ID', 'REF', 'ALT',
                              'QUAL', 'FILTER', 'INFO',
                              'FORMAT', 'SAMPLE'))
    
    args_list = [(i, _type, minDP,  remove_snp, 
                    remove_indel, out) for i in vcf_df]
    with multiprocessing.Pool(threads) as pool:
        res = pool.map(filter_vcf, args_list)
  
    logging.debug('Loaded vcf file of `{}`'.format(vcf_path))
    return pd.concat(res)


def import_gff(gff_path, _type=None):
    compression = 'gzip' if gff_path.endswith('.gz') else 'infer'
    gff_df = pd.read_csv(gff_path, sep='\t', 
                         header=None, 
                         index_col=0, 
                         comment="#", 
                         compression = compression,
                         names=('chrom', 'source',
                                'type', 'start', 
                                'end', 'score', 
                                'strand','phase',
                                'attributes'))
    if _type in set(gff_df.head(100)['type']):
        if _type: 
            gff_df = gff_df[gff_df['type'] == _type]
    else:
        logging.warning('Warning: Failed to filter data by' 
              'type, input type is not correct')
    
    logging.debug('Loaded gff file of `{}`'.format(gff_path))
    return gff_df


def get_gene(attributes):
    db = dict(map(lambda x: x.split('='), 
                  [i for i in attributes.split(';') if i]))
    return db['ID']


def calc_heter_by_chrom(args):
    chrom, tmp_vcf_df, genes, starts, ends = args
    
    res = []

    for gene, (start, end) in zip(genes ,zip(starts, ends)):
        
        nums = len(tmp_vcf_df[(tmp_vcf_df.POS >= start) & \
                            (tmp_vcf_df.POS <= end)])
        length = end - start + 1
        res.append([chrom, start, end, gene, nums, 
                    length, nums*1.0 / length])
    logging.debug('successful run `{}`'.format(chrom))
    
    return res
    

def calc_heter_all(vcf_path, 
                gff_path, 
                out,
                threads=12,
                minDP=0,
                _type='gene',
                include_indel=False):
    
    remove_indel = False if include_indel else True
    
    gff_df = import_gff(gff_path, _type=_type)
    vcf_df = import_vcf(vcf_path, 
                        threads=threads,
                        minDP=minDP,
                        remove_indel=remove_indel)
    
    vcf_chroms = [i for i in set(vcf_df.index.to_list())
                if i in set(gff_df.index)]
    args_list = []
    for chrom in vcf_chroms:
        try:
            args_list.append((chrom, vcf_df.loc[chrom],
                   gff_df.loc[chrom].attributes.map(get_gene),
                   gff_df.loc[chrom].start,
                   gff_df.loc[chrom].end))

        except AttributeError:
            genes = [get_gene(gff_df.loc[chrom].attributes)]
            args_list.append((chrom, vcf_df.loc[chrom],
                                genes, [gff_df.loc[chrom].start],
                                [gff_df.loc[chrom].end]))

    #args_list = [ (chrom, vcf_df.loc[chrom],
    #               gff_df.loc[chrom].attributes.map(get_gene),
     #              gff_df.loc[chrom].start,
     #              gff_df.loc[chrom].end)
    #               for chrom in vcf_chroms]
    with multiprocessing.Pool(threads) as pool:
        res = pool.map(calc_heter_by_chrom, args_list)
    
    res = np.concatenate(res)
    res_df = pd.DataFrame(res, columns=['chrom', 'start', 
                                'end', 'gene','nums', 
                                'length', 'rate'])
    ref_df = res_df.sort_values(by=['chrom', 'start'])
    res_df.to_csv(out, sep='\t', header=None, index=None)
    cmd = 'sort -V {0} -o {0}'.format(out)
    os.system(cmd)
    logging.debug('Successful, per chrom heter result is in `{}`'.format(out))
    nums = res_df.nums.apply(int).sum()
    length = res_df.length.apply(int).sum()
    rate = 1.0 * nums / length
    
    stat_out = sys.stdout if out is sys.stdout \
                    else open(out.rsplit('.', 1)[0] + '.stat', 'w')
    
    print('\t'.join(('SNP_numbers', 'Gene_length', 'Heterozygosity')),
            file=stat_out)
    print('\t'.join(map(str, (nums, length, rate))), 
            file=stat_out)
    logging.debug('Successful, total stat is in `{}`'.format(stat_out.name))



if __name__ == "__main__":
    p = argparse.ArgumentParser(prog=__file__,
                        description=__doc__,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('vcf', help='snp vcf file')
    pReq.add_argument('gff', help='annotation gff file')
    pOpt.add_argument('-t', '--threads', type=int, default=12,
            help='threads of program [default: %(default)s]')
    pOpt.add_argument('--minDP', type=int, default=0,
            help='minimum DP value [default: %(default)s]')
    pOpt.add_argument('--include_indel', action='store_true',
            default=False, 
            help='include indel to calculation [default: %(default)s]')
    pOpt.add_argument('-o', '--out', default=sys.stdout, 
            help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args()

    calc_heter_all(args.vcf,
                args.gff,
                args.out,
                threads=args.threads,
                minDP=args.minDP,
                include_indel=args.include_indel)