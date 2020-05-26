#!/bin/bash


bed=$1
matrix=$2
genome=$3
chromsizes=$4
outdir=$5
thread=$6
if [[ -z $bed || -z $matrix || -z $genome || -z $chromsizes ]];then
        echo
        echo "Usage: `basename $0` rice_100000_abs.bed rice_100000_iced.matrix rice chrom.sizes"
        echo
        exit;
fi

if [[ -z $outdir ]]; then
        outdir="./"
fi

if [[ -z $thread ]]; then
        thread=12
fi

if [ ! -f ${cworld_dir}/${genome}.refseq.txt ];then
        echo "No such file of $genome.refseq.txt in $cworld_dir"
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

name=`basename ${bed} | sed 's/_abs.bed//g'`
prefix=`basename ${matrix} | sed 's/.matrix//g'`
sparseToDense.py -b $bed $matrix -o ${name}.dense.matrix  -c


mv *${name}.dense.matrix ${outdir}
echo "staring create cworld header"
cworld_header.py $bed $genome -o ${outdir} -c


for dense in ${outdir}/*${name}*dense.matrix; do
        echo ${dense%%.dense.matrix}
done | parallel -j ${thread} "addMatrixHeaders.pl -i {}.dense.matrix --xhf {}_abs.bed.header --yhf {}_abs.bed.header -v -o {}" && \
for matrix in ${outdir}/*${name}*addedHeaders.matrix.gz;do
        echo ${matrix}
done | parallel -j ${thread} "matrix2compartment.pl -i {}" && \
mv *${name}.addedHeaders* ${outdir} && \
cat ${outdir}/*${name}*addedHeaders.zScore.eigen1.bedGraph |sort -V |grep -v track > ${outdir}/${name}_all_eigen1.bg && \
bedtools intersect -a ${outdir}/${name}_all_eigen1.bg -b ${cworld_dir}/${genome}.refseq.txt -c > ${outdir}/${prefix}_all_eigen1_gene_density.bg && \
cut -f 1-3,5 ${outdir}/${prefix}_all_eigen1_gene_density.bg > ${outdir}/${prefix}_gene_density.bg && \
ab_boxplot.py ${outdir}/${prefix}_all_eigen1_gene_density.bg --xlabel 'Gene' --ylabel 'Count' -o ${outdir}/${prefix}_all_eigen1_gene_density.pdf &
python -m TDGP.analysis.ab quickPlot ${outdir}/${name}_all_eigen1.bg $chromsizes -g ${outdir}/${name}_gene_density.bg -o ${outdir}/quickPlot | parallel -j ${N_CPU} {} &
wait