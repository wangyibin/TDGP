#!/usr/bin/env python
# -*- coding:utf-8 -*-

import os
import sys

import cooler
import numpy as np


class CalObsExp(object):
    def __init__(self, matrix_source, outpre):
        cool = cooler.Cooler(matrix_source)
        for c in cool.chromnames:
            raw = cool.matrix(balance=False, sparse=False).fetch(c)
            raw[np.isnan(raw)] = 0
            expected = self.expected_matrix(raw)
            obs_exp = raw / expected
            obs_exp[expected == 0] = 0
            outfil = outpre + '.{0}.npy'.format(c)
            np.save(outfil, obs_exp)

    def expected_matrix(self, raw):
        tmp = raw.sum(axis=0) != 0  # valid rows or columns
        n = raw.shape[0]
        expected = np.zeros_like(raw)
        idx = np.arange(n)
        for i in idx:
            if i > 0:
                valid = tmp[:-i] * tmp[i:]
            else:
                valid = tmp
            current = raw.diagonal(i)[valid]
            if current.size > 0:
                v = current.mean()
                if i > 0:
                    expected[idx[:-i], idx[i:]] = v
                    expected[idx[i:], idx[:-i]] = v
                else:
                    expected[idx, idx] = v
        return expected


if __name__ == '__main__':
    datasets = sys.argv[1]
    outpre = sys.argv[2]
    work = CalObsExp(datasets, outpre)
