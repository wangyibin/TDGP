#!/usr/bin/env python

from __future__ import print_function

import gffutils
import sqlite3
import sys

def gff2bed(gfffile, dbname, feature='gene', target='ID'):
    # try to create database, but if it already exists,
    # it will use the existing one.
    try:
        db = gffutils.create_db(gfffile, dbname)
    except sqlite3.OperationalError:
        db = gffutils.FeatureDB(dbname)

    for item in db.features_of_type(feature):
        tmp, = item[target]
        out = [item.seqid, item.start, item.end, tmp]
        print("\t".join(map(str, out)))

if __name__ == "__main__":
    from optparse import OptionParser
    p = OptionParser(__doc__)
    p.add_option('-d','--dbname', dest='dbname', default='*gff.db',
                 help='the database of the gff file, if existing it will'
                      'use it')
    p.add_option("--feature", dest='feature', default='gene',
          help='extract feature of annotations file [default: %default')
    p.add_option("--target", dest='target', default='ID',
          help="the target of feature [default: %default]")

    opts, args = p.parse_args()
    if len(args) != 1:
        sys.exit(p.print_help())

    gfffile, = args
    dbname = opts.dbname
    if dbname == "*gff.db":
        dbname = gfffile.replace(".gz", ".db")
    gff2bed(gfffile, dbname=dbname,
            feature=opts.feature, target=opts.target)
