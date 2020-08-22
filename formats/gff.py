#!/usr/bin/env python
# -*- coding:utf-8 -*-

"""
Util tools for gff formats.
"""


from __future__ import print_function

import argparse
import logging
import os
import os.path as op
import sys
import pandas as pd
import numpy as np
import gffutils

from collections import OrderedDict

from TDGP.apps.base import ActionDispatcher
from TDGP.apps.base import check_file_exists, debug
from TDGP.apps.base import listify


def main():

    actions = (
            ("renameAttributes", "rename attributes by a two column list file"),
        )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


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

def get_gene_id(attributes):
    db = dict(map(lambda x: x.split('='), 
                  [i for i in attributes.split(';') if i]))
    return db['ID']


def createDB(gff, dbfn=None):
    """
    create a db file for gffutils

    Params:
    --------
    gff: `str` gff file
    dbfn: `str` db name [default: gff.db]

    Returns:
    --------
    db: `gffutils.FeatureDB` db for gffutils

    Examples:
    ---------
    >>> createDB('in.gff3')
    <FeatureDB ...>
    """
    ## create gff db
    if not dbfn: 
        db_name = gff + ".db"
    else:
        db_name = dbfn
    if not op.exists(db_name):
        logging.debug('No such database file of `{}`, creating ...'.format(db_name))
        gffutils.create_db(gff, dbfn=db_name, keep_order=True)
    else:
        logging.debug('Already exists DB file of `{}`, skip.'.format(db_name))
    db = gffutils.FeatureDB(db_name)
    
    return db


## outside command 
def renameAttributes(args):
    """
    %(prog)s <in.gff3> <rename.list> [Options]

        Change the attributes within the gff3
    """
    p = p=argparse.ArgumentParser(prog=renameAttributes.__name__,
                        description=renameAttributes.__doc__,
                        formatter_class=argparse.RawTextHelpFormatter,
                        conflict_handler='resolve')
    pReq = p.add_argument_group('Required arguments')
    pOpt = p.add_argument_group('Optional arguments')
    pReq.add_argument('gff', 
            help='gff file ')
    pReq.add_argument('renameList', 
            help='rename list, two columns <old_gene_name\tnew_gene_name>')
    pOpt.add_argument('-d', '--dbname', default=None,
            help='gff database name [default: gff3.db')
    pOpt.add_argument('-o', '--output', type=argparse.FileType('w'),
            default=sys.stdout, help='output file [default: stdout]')
    pOpt.add_argument('-h', '--help', action='help',
            help='show help message and exit.')
    
    args = p.parse_args(args)
    output = args.output
    db = createDB(args.gff, args.dbname)
    rename_db = OrderedDict(i.strip().split() for i in open(args.renameList)
                        if i.strip())
    
    for ID in rename_db:
        print(str(db[ID]).replace(ID, rename_db[ID]), file=output)
        for feature in db.children(ID, order_by='start'):
            print(str(feature).replace(ID, rename_db[ID]), file=output)
        print("", file=output)
    logging.debug("Successful. Output file is in `{}`".format(output.name))



if __name__ == "__main__":
    main()