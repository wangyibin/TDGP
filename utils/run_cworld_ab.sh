#!/bin/bash


bed=$1
matrix=$2
genome=$3

if [[ -z $bed || -z $matrix || -z $genome ]];then
        echo
        echo "Usage: `basename $0` rice_100000_abs.bed rice_100000_iced.matrix rice"
        echo
        exit;
fi

cworld_dir=/public1/home/stu_wangyibin/software/cworld-dekker/scripts/util/geneDensity/
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
python /public1/home/stu_wangyibin/code/TDGP/utils/cworld_header.py $bed $genome -c


for dense in *dense.matrix; do
        echo ${dense%%.dense.matrix}
done | parallel -j 12 "addMatrixHeaders.pl -i {}.dense.matrix --xhf {}_abs.bed.header --yhf {}_abs.bed.header -v -o {}"

for matrix in *addedHeaders.matrix.gz;do
        echo ${matrix}
done | parallel -j 12 "matrix2compartment.pl -i {}"

cat *addedHeaders.zScore.eigen1.bedGraph |sort -V |grep -v track > ${matrix%%matrix}all_eigen1.bg

bedtools intersect -a ${matrix%%matrix}all_eigen1.bg -b $cworld_path/${genome}.refseq.txt -c > ${matrix%%.matrix}_all_eigen1_gene_density.bg
