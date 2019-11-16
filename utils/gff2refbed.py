#!/usr/bin/env python
# -*- coding:utf-8 -*-


import gzip
import sys

from collections import OrderedDict


def deal_gff_info(info):
    return dict(map(lambda x: x.split("="), info.split(";")))


def read_gff(gfffile):
    if gfffile[-2:] == "gz":
        openfile = gzip.open(gfffile)
    else:
        openfile = open(gfffile)
    name = ''
    db = OrderedDict()
    flag = 1
    with openfile as fp:
        for line in fp:
            if line[0] == "#":
                continue

            line_list = line.strip().split("\t")
            if line_list[2] == "mRNA":
                flag = 1
                name = deal_gff_info(line_list[8])["Name"]
                db[name] =[line_list[0],  line_list[3], line_list[4],
                         line_list[3], line_list[4], line_list[6], name, name, "", [], []]
           # elif line_list[2] == 'mRNA' and flag == 1:
           #     db[name][3] = line_list[3]
           #     db[name][4] = line_list[4]
            elif line_list[2] == 'CDS':
                db[name][9].append(line_list[3])
                db[name][10].append(line_list[4])

    return db

def out_result(gfffile):
    db = read_gff(gfffile)
    for gene in db:
        db[gene][9] = ",".join(db[gene][9])
        db[gene][10] = ",".join(db[gene][10])
        print("\t".join(db[gene]))


if __name__ == "__main__":
    out_result(sys.argv[1])
    

