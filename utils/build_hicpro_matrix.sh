#!/bin/bash


allpairs=$1
binsize=$2
chrom_size=$3
outprefix=$4

if [[ -z $allpairs || -z $binsize || -z $chrom_size || -z $outprefix ]]; then
        echo
        echo "Usage: $0 allValidPairs 100000 chrom_reference.sizes sample/raw/100000/sample_100000"
        echo
        exit;
fi

mkdir -p `dirname $outprefix`
mkdir -p `dirname $outprefix | sed 's/raw/iced/g'`
iced_outprefix=`echo $outprefix | sed 's/raw/iced/g'`
#cat $allpairs | /public1/home/stu_wangyibin/software/HiC-Pro_2.11.1/scripts/build_matrix --matrix-format upper --binsize $binsize --chrsizes $chrom_size -ifile /dev/stdin --oprefix $outprefix

python /public1/home/stu_wangyibin/software/HiC-Pro_2.11.1/scripts/ice --results_filename ${iced_outprefix}_iced.matrix --filter_low_counts_perc 0.02 --filter_high_counts_perc 0 --max_iter 100 --eps 0.1 --remove-all-zeros-loci --output-bias 1 --verbos 1 ${outprefix}.matrix


