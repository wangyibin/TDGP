#!/bin/sh


cool=$1
bwPath=$2

if [[ -z ${cool} && ! -f ${cool} ]]; then
    echo "Please input correct cool file path"
    echo "Uasge: `basename $0` <sample_100000_iced.cool> [/path/to/bigwig]"
    echo "Usage: `basename $0` sample_100000_iced.cool"
    echo "Usage: `basename $0` ample_100000_iced.cool ../../../../100000"
    exit
fi


#plotting without bw, set bwPaht to empty

resolution=`basename $cool | perl -lne '/_(\d+)[_.]/ && print $1'`
prefix=`basename $cool | sed 's/.cool//g'`

suffix=""
if [ ! -z $bwPath ]; then
    suffix=""" --bigwig ${bwPath}/${resolution}/${prefix}_all_eigen1.bw \
    ${bwPath}/${resolution}/${prefix}_gene_density.bw \
    ${bwPath}/${resolution}/${prefix}_RNA_log1p_density.bw \
    ${bwPath}/${resolution}/${prefix}_Retro_density.bw \
    ${bwPath}/${resolution}/${prefix}_DNA_density.bw \
    --pyLabel 'Compartments' 'Gene' 'RNA' 'Retro-TE' 'DNA-TE'
    """
fi
for i in {1..8}; do 
echo "hicPlotMatrix --matrix ${cool} --dpi 300 --log1p \
    --chromosomeOrder Chr${i}A Chr${i}B Chr${i}C Chr${i}D \
    -o ${prefix}_Chr${i}.pdf --clearMaskedBins \
    ${suffix}"
done | tee ${PWD}/`basename ${0%%.sh}`_pdf.commands | parallel -j 8 {}

for i in {1..8}; do 
echo "hicPlotMatrix --matrix ${cool} --dpi 300 --log1p \
    --chromosomeOrder Chr${i}A Chr${i}B Chr${i}C Chr${i}D \
    -o ${prefix}_Chr${i}.png --clearMaskedBins \
    ${suffix}"
done | tee ${PWD}/`basename ${0%%.sh}`_png.commands | parallel -j 8 {}