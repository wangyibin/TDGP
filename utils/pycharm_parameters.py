#!/usr/bin/env python
# -*- coding:utf-8 -*-

# @Time: 2019/6/9 15:48

import sys
import os.path as op




if __name__ == '__main__':
    from optparse import OptionParser
    p = OptionParser(__doc__)
    p.add_option("-d",dest="d")
    opts, args = p.parse_args()
    print(args)
    print(opts)
    s = [i.strip() for i in open(opts.d)]
    print(s)