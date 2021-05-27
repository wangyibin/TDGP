#!/bin/bash 

# calculate a feature density of a genome in particular window

bed=$1
size=$2
window=$3

if [[ -z $bed || -z $size || -z $window ]]; then
    echo 
    echo -e "Usage: `basename $0` <bed> <chrom.sizes> <window>"
    echo -e "           calculate density for a bed file in window size"
    exit;
fi

PID=$!

bedtools makewindows -g ${size} -w ${window} > ${PID}.${window}.window
bedtools intersect -a ${PID}.${window}.window -b ${bed} -c > ${bed%%.bed}_${window}.bg 
bedGraphToBigWig  ${bed%%.bed}_${window}.bg ${size}  ${bed%%.bed}_${window}.bw


rm ${PID}.${window}.window