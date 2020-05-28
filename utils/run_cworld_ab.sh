#!/bin/bash


bed=$1
matrix=$2
genome=$3
chromsizes=$4
outdir=$5
thread=$6
TE_DATA=$7

resolution=`basename $matrix | perl -lne '/_(\d+)[_.]/ && print $1'`
RES_FILE_NAME=`basename ${bed} | sed 's/_abs.bed//g' | sed "s/_${resulotion}//g"`
name=`basename ${bed} | sed 's/_abs.bed//g'`
prefix=`basename ${matrix} | sed 's/.matrix//g'`
if [[ -z $bed || -z $matrix || -z $genome || -z $chromsizes ]];then
        echo
        echo "Usage: `basename $0` rice_100000_abs.bed rice_100000_iced.matrix rice chrom.sizes"
        echo "Usage: `basename $0` rice_100000_abs.bed rice_100000_iced.matrix rice chrom.sizes outdir thread TE_DATA"
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


sparseToDense(){
        echo "starting sparse to dense"
        sparseToDense.py -b $bed $matrix -o ${name}.dense.matrix  -c

        if [[ ${outdir} != "./" ]]; then mv *${name}.dense.matrix ${outdir}/cworld_results; fi
}

cworldheader(){
        echo "staring create cworld header"
        cworld_header.py $bed $genome -o ${outdir}/cworld_results -c
        echo "created header"
}


addheader(){
        echo "starting analysis compartments"
        for dense in ${outdir}/cworld_results/*${name}*dense.matrix; do
                echo ${dense%%.dense.matrix}
        done | parallel -j ${thread} "addMatrixHeaders.pl -i {}.dense.matrix --xhf {}_abs.bed.header --yhf {}_abs.bed.header -v -o {}" 
}
matrix2compartment(){
        cd  ${outdir}/cworld_results
        for mtrx in *${name}*addedHeaders.matrix.gz;do
                echo ${mtrx}
        done | parallel -j ${thread} "matrix2compartment.pl -i {} " && \
        cd -
        echo "matrix2compartments done"
} 

mergeBedGraph(){
        echo "merge all chromosome results"
        cat ${outdir}/cworld_results/*${name}*addedHeaders.zScore.eigen1.bedGraph |sort -V |grep -v track > ${outdir}/${prefix}_all_eigen1.bg && \
        bedGraphToBigWig ${outdir}/${prefix}_all_eigen1.bg ${chromsizes} ${outdir}/${prefix}_all_eigen1.bw
        echo "Merge results done"
}


TE_analysis(){
        local te_data=$1
        echo "Starting analysis the TE density"
        #bedtools makewindows -g ${chromsizes} -w ${bsize} > ${outdir}/${RES_FILE_NAME}_${bsize}.window && \
        #bedtools intersect -a ${outdir}/${RES_FILE_NAME}_${bsize}.window -b ${te_data} -c > ${outdir}/${RES_FILE_NAME}_${bsize}_TE.bg
        bedtools intersect -a ${outdir}/${prefix}_all_eigen1.bg -b ${te_data} -c > ${outdir}/${prefix}_all_eigen1_TE_density.bg && \
        cut -f 1-3,5 ${outdir}/${prefix}_all_eigen1_TE_density.bg > ${outdir}/${prefix}_TE_density.bg && \
        bedGraphToBigWig ${outdir}/${prefix}_TE_density.bg ${chromsizes} ${outdir}/${prefix}_TE_density.bw
        ab_boxplot.py ${outdir}/${prefix}_all_eigen1_TE_density.bg --xlabel 'TE' --ylabel 'Count' -o ${outdir}/${prefix}_all_eigen1_TE_density.pdf
        ab_dotplot.py ${outdir}/${prefix}_all_eigen1_gene_density.bg --xlabel 'Gene' --ylabel 'Count' -o ${outdir}/${prefix}_all_eigen1_TE_density_dotplot.pdf
        echo "TE analysis done"
}

GENE_analysis(){
        echo "Start analysis the gene density"
        bedtools intersect -a ${outdir}/${prefix}_all_eigen1.bg -b ${cworld_dir}/${genome}.refseq.txt -c > ${outdir}/${prefix}_all_eigen1_gene_density.bg && \
        cut -f 1-3,5 ${outdir}/${prefix}_all_eigen1_gene_density.bg > ${outdir}/${prefix}_gene_density.bg && \
        bedGraphToBigWig ${outdir}/${prefix}_gene_density.bg ${chromsizes} ${outdir}/${prefix}_gene_density.bw
        ab_boxplot.py ${outdir}/${prefix}_all_eigen1_gene_density.bg --xlabel 'Gene' --ylabel 'Count' -o ${outdir}/${prefix}_all_eigen1_gene_density.pdf &
        ab_dotplot.py ${outdir}/${prefix}_all_eigen1_gene_density.bg --xlabel 'Gene' --ylabel 'Count' -o ${outdir}/${prefix}_all_eigen1_gene_density_dotplot.pdf
        echo "Gene analysis done"
}




mkdir -p ${outdir}/cworld_results

sparseToDense 
cworldheader
addheader

matrix2compartment 
mergeBedGraph
GENE_analysis &

if [[ -e ${TE_DATA} ]]; then
TE_analysis ${TE_DATA} &  
TE_suffix="--TE ${outdir}/${prefix}_TE_density.bg"   
fi
wait

#python -m TDGP.analysis.ab quickPlot ${outdir}/${prefix}_all_eigen1.bg $chromsizes -g ${outdir}/${prefix}_all_eigen1_gene_density.bg  ${TE_suffix} -o ${outdir}/quickPlot | parallel -j ${thread} {} &

wait
