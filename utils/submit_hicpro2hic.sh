#!/bin/bash


cmd="hicpro2juicebox.sh -i *.allValidPairs -g chrom_reference.sizes -j /public1/home/stu_wangyibin/software/juicer/scripts/juicer_tools.1.7.6_jcuda.0.8.jar -r *I.bed "

echo $cmd > run_hicpro2hic.sh
qsub -pe mpi 1 -j y -cwd -S /bin/bash run_hicpro2hic.sh
