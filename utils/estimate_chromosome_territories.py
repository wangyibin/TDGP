#!/usr/bin/env python
# -*- coding:utf-8 -*-

import warnings 
warnings.simplefilter(action="ignore", category=RuntimeWarning)
warnings.simplefilter(action="ignore", category=PendingDeprecationWarning)

from math import log2
from scipy.sparse import csr_matrix
import numpy as np
import sys

from hicmatrix import HiCMatrix as hm
from hicexplorer.utilities import obs_exp_matrix, expected_interactions
from hicexplorer.utilities import convertNansToZeros, convertInfsToZeros

import logging
log = logging.getLogger(__name__)


def _obs_exp(pSubmatrix):
    obs_exp_matrix_ = obs_exp_matrix(pSubmatrix)
    obs_exp_matrix_ = convertNansToZeros(csr_matrix(obs_exp_matrix_))
    obs_exp_matrix_ = convertInfsToZeros(csr_matrix(obs_exp_matrix_)).todense()
    return obs_exp_matrix_


def get_expected_matrix(pSubmatrix):
    expected_interactions_in_distance = expected_interactions(pSubmatrix)
    row, col = pSubmatrix.nonzero()
    distance = np.ceil(np.absolute(row - col) / 2).astype(np.int32)
    expected = expected_interactions_in_distance[distance]
    pSubmatrix.data = expected
    pSubmatrix = convertNansToZeros(csr_matrix(pSubmatrix))
    pSubmatrix = convertInfsToZeros(csr_matrix(pSubmatrix)).todense()

    return pSubmatrix


def get_observed_matrix(pSubmatrix):
    return pSubmatrix.todense()
    

def get_whole_chrom_matrix(coolfile):
    hic_ma = hm.hiCMatrix(pMatrixFile=coolfile)
    #submatrix = _obs_exp(hic_ma.matrix)
    print(hic_ma.matrix.todense())
    obsmatrix = get_observed_matrix(hic_ma.matrix)
    print(obsmatrix)
    expmatrix = get_expected_matrix(hic_ma.matrix)
    print(expmatrix)
    chrom_list = []
    for chrom in hic_ma.getChrNames():
        if chrom[:3] == 'Chr':
            chrom_list.append(chrom)
    whole_chrom_matrix = np.zeros((len(chrom_list), len(chrom_list)))
    for i, chrom1 in enumerate(chrom_list):
        for j, chrom2 in enumerate(chrom_list):
            if chrom1 != chrom2:
                chrom1_bin_range1, chrom1_bin_range2 = hic_ma.getChrBinRange(chrom1)
                chrom2_bin_range1, chrom2_bin_range2 = hic_ma.getChrBinRange(chrom2)
                #chrom1_chrom2 = submatrix[chrom1_bin_ranges[0]: chrom1_bin_ranges[1]+1, chrom2_bin_ranges[0]: chrom2_bin_ranges[1]+1]
                #size = chrom1_chrom2.shape[0]*chrom1_chrom2.shape[1]
                obs = obsmatrix[chrom1_bin_range1: chrom1_bin_range2+1, chrom2_bin_range1: chrom2_bin_range2+1].sum()
                exp = expmatrix[chrom1_bin_range1: chrom1_bin_range2+1, chrom2_bin_range1: chrom2_bin_range2+1].sum()
                whole_chrom_matrix[i, j] = log2(obs/exp)

    return whole_chrom_matrix


def plot_heatmap(matix):
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots()
    im = ax.imshow(matrix, cmap='coolwarm', vmax=3, vmin=-3)
    fig.colorbar(im, pad=0.01)
    plt.savefig('heatmap.pdf', dpi=300)


if __name__ == "__main__":
    matrix = get_whole_chrom_matrix(sys.argv[1])
    plot_heatmap(matrix)
