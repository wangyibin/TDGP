#!/bin/bash


bed=$1
matrix=$2
genome=$3
chromsizes=$4
outdir=$5
thread=$6
RNA_DATA=$7
#TE DATA BED FILE
TE_DATA=$8
if [[ -z $bed || -z $matrix || -z $genome || -z $chromsizes ]];then
        echo
        echo "Usage: `basename $0` rice_100000_abs.bed rice_100000_iced.matrix rice chrom.sizes"
        echo "Usage: `basename $0` rice_100000_abs.bed rice_100000_iced.matrix rice chrom.sizes outdir thread RNA_DATA TE_DATA"
        echo "      TE_DATA: four column bed forth column is the type of TE: LTR DNA...."
        echo "      RNA_DATA: bam file of RNA-Seq"
        exit;
fi
resolution=`basename $matrix | perl -lne '/_(\d+)[_.]/ && print $1'`
RES_FILE_NAME=`basename ${bed} | sed 's/_abs.bed//g' | sed "s/_${resulotion}//g"`
name=`basename ${bed} | sed 's/_abs.bed//g'`
prefix=`basename ${matrix} | sed 's/.matrix//g'`


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


RNA_analysis(){
        local rna_data=$1
        echo "Staring analysis the RNA-Seq data"
        rna_name=`basename ${rna_data}`
        bigwig=${rna_name%.bam}_RNA_${resolution}.bw
        bedgraph=${bigwig%.bw}.bg
        bamCoverage -p ${thread} -b ${rna_data} --normalizeUsing RPKM -bs ${resolution} -o ${outdir}/${bigwig}
        bigWigToBedGraph ${outdir}/${bigwig} ${outdir}/${bedgraph}
        bedtools intersect -a ${outdir}/${prefix}_all_eigen1.bg -b ${outdir}/${bedgraph} -wao | awk '{print $1,$2,$3,$4,$8}' OFS='\t'> ${outdir}/${prefix}_all_eigen1_RNA_density.bg
        cut -f 1-3,5 ${outdir}/${prefix}_all_eigen1_RNA_density.bg > ${outdir}/${prefix}_RNA_density.bg
        python -c "import pandas as pd;import numpy as np;df=pd.read_csv('${outdir}/${prefix}_RNA_density.bg', sep='\t',header=None, index_col=None);df[3] = np.log1p(df[3]); df.to_csv('${outdir}/${prefix}_RNA_log1p_density.bg', header=None,index=None, sep='\t');"
        bedGraphToBigWig ${outdir}/${prefix}_RNA_log1p_density.bg ${chromsizes} ${outdir}/${prefix}_RNA_log1p_density.bw 
        ab_boxplot.py ${outdir}/${prefix}_all_eigen1_RNA_density.bg --xlabel 'RNA' --ylabel 'RPKM' -o ${outdir}/${prefix}_all_eigen1_RNA_density_barplot.pdf
        ab_dotplot.py ${outdir}/${prefix}_all_eigen1_RNA_density.bg --log1p --topA 10 --topB 10  --title 'RNA' --xlabel 'PC1' --ylabel 'RPKM' -o ${outdir}/${prefix}_all_eigen1_RNA_density_dotplot.pdf
}
TE_analysis(){
        local te_data=$1
        echo "Starting analysis the TE density"
        #bedtools makewindows -g ${chromsizes} -w ${bsize} > ${outdir}/${RES_FILE_NAME}_${bsize}.window && \
        #bedtools intersect -a ${outdir}/${RES_FILE_NAME}_${bsize}.window -b ${te_data} -c > ${outdir}/${RES_FILE_NAME}_${bsize}_TE.bg
        te_name=`basename ${te_data}`
        egrep "LTR|LINE|SINE" ${te_data} > ${outdir}/${te_name%.bed}.Retro.bed
        egrep "DNA|RC" ${te_data} > ${outdir}/${te_name%.bed}.DNA.bed
        bedtools intersect -a ${outdir}/${prefix}_all_eigen1.bg -b ${outdir}/${te_name%.bed}.Retro.bed -c > ${outdir}/${prefix}_all_eigen1_Retro_density.bg && \
        cut -f 1-3,5 ${outdir}/${prefix}_all_eigen1_Retro_density.bg > ${outdir}/${prefix}_Retro_density.bg && \
        bedGraphToBigWig ${outdir}/${prefix}_Retro_density.bg ${chromsizes} ${outdir}/${prefix}_Retro_density.bw
        ab_boxplot.py ${outdir}/${prefix}_all_eigen1_Retro_density.bg --xlabel 'Retro-TE' --ylabel 'Count' -o ${outdir}/${prefix}_all_eigen1_Retro_density_barplot.pdf
        ab_dotplot.py ${outdir}/${prefix}_all_eigen1_Retro_density.bg --title 'Retro-TE' --xlabel 'PC1' --ylabel 'Count' -o ${outdir}/${prefix}_all_eigen1_Retro_density_dotplot.pdf

        bedtools intersect -a ${outdir}/${prefix}_all_eigen1.bg -b ${outdir}/${te_name%.bed}.DNA.bed -c > ${outdir}/${prefix}_all_eigen1_DNA_density.bg && \
        cut -f 1-3,5 ${outdir}/${prefix}_all_eigen1_DNA_density.bg > ${outdir}/${prefix}_DNA_density.bg && \
        bedGraphToBigWig ${outdir}/${prefix}_DNA_density.bg ${chromsizes} ${outdir}/${prefix}_DNA_density.bw
        ab_boxplot.py ${outdir}/${prefix}_all_eigen1_DNA_density.bg --xlabel 'DNA-TE' --ylabel 'Count' -o ${outdir}/${prefix}_all_eigen1_DNA_density_barplot.pdf
        ab_dotplot.py ${outdir}/${prefix}_all_eigen1_DNA_density.bg --title 'DNA-TE' --xlabel 'PC1' --ylabel 'Count' -o ${outdir}/${prefix}_all_eigen1_DNA_density_dotplot.pdf
        echo "TE analysis done"
}

GENE_analysis(){
        echo "Start analysis the gene density"
        bedtools intersect -a ${outdir}/${prefix}_all_eigen1.bg -b ${cworld_dir}/${genome}.refseq.txt -c > ${outdir}/${prefix}_all_eigen1_gene_density.bg && \
        cut -f 1-3,5 ${outdir}/${prefix}_all_eigen1_gene_density.bg > ${outdir}/${prefix}_gene_density.bg && \
        bedGraphToBigWig ${outdir}/${prefix}_gene_density.bg ${chromsizes} ${outdir}/${prefix}_gene_density.bw
        ab_boxplot.py ${outdir}/${prefix}_all_eigen1_gene_density.bg --xlabel 'Gene' --ylabel 'Count' -o ${outdir}/${prefix}_all_eigen1_gene_density.pdf &
        ab_dotplot.py ${outdir}/${prefix}_all_eigen1_gene_density.bg --title 'Gene' --xlabel 'PC1' --ylabel 'Count' -o ${outdir}/${prefix}_all_eigen1_gene_density_dotplot.pdf
        echo "Gene analysis done"
}

mkdir -p ${outdir}/cworld_results

sparseToDense 
cworldheader
addheader

matrix2compartment 
mergeBedGraph
GENE_analysis &

if [[ -e ${RNA_DATA} ]]; then
RNA_analysis ${RNA_DATA} &
RNA_suffix="--RNA ${outdir}/${prefix}_RNA_density.bg"
fi

if [[ -e ${TE_DATA} ]]; then
TE_analysis ${TE_DATA} &  
TE_suffix="--Retro ${outdir}/${prefix}_Retro_density.bg --DNA ${outdir}/${prefix}_DNA_density.bg"   
fi
wait

python -m TDGP.analysis.ab quickPlot ${outdir}/${prefix}_all_eigen1.bg $chromsizes \
        -g ${outdir}/${prefix}_all_eigen1_gene_density.bg  ${TE_suffix} ${RNA_suffix} \
        -o ${outdir}/quickPlot | parallel -j ${thread} {}