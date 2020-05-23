#!/bin/bash


bed=$1
matrix=$2
genome=$3
chromsizes=$4
if [[ -z $bed || -z $matrix || -z $genome || -z $chromsizes ]];then
        echo
        echo "Usage: `basename $0` rice_100000_abs.bed rice_100000_iced.matrix rice chrom.sizes"
        echo
        exit;
fi


if [ ! -f ${cworld_dir}/${genome}.refseq.txt ];then
        echo "No such file of $genome.refseq.txt in $cworld_dif"
        echo "to execute command of :"
        echo "
                    ## extract gene density bed from gff file
                    python ~/code/TDGP/utils/gff2bed.py Osativa_323_v7.0.gene.gff3.gz > rice.refseq.txt
                    ## move above txt to cowrld dictory
                    cp rice.refseq.txt $cworld_dir
        "
        exit;
fi



echo "starting sparse to dense"
sparseToDense.py -b $bed $matrix -o ${bed%%_abs.bed}.dense.matrix  -c
echo "staring create cworld header"
cworld_header.py $bed $genome -c

name=${bed%%_abs.bed}
for dense in *${name}*dense.matrix; do
        echo ${dense%%.dense.matrix}
done | parallel -j 12 "addMatrixHeaders.pl -i {}.dense.matrix --xhf {}_abs.bed.header --yhf {}_abs.bed.header -v -o {}"

for matrix in *${name}*addedHeaders.matrix.gz;do
        echo ${matrix}
done | parallel -j 12 "matrix2compartment.pl -i {}"

cat *${name}*addedHeaders.zScore.eigen1.bedGraph |sort -V |grep -v track > ${matrix%%.matrix}_all_eigen1.bg


bedtools intersect -a ${matrix%%.matrix}_all_eigen1.bg -b $cworld_dir/${genome}.refseq.txt -c > ${matrix%%.matrix}_all_eigen1_gene_density.bg
cut -f 1-3,5 ${matrix%%.matrix}_all_eigen1_gene_density.bg > ${matrix%%.matrix}_gene_density.bg

python -m TDGP.analysis.ab quickPlot ${matrix%%.matrix}_all_eigen1.bg $chromsizes -g ${matrix%%.matrix}_gene_density.bg | parallel -j 12 {}