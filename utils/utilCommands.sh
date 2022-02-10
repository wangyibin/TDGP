#!/bin/bash



## my usualy use command library

why(){
        echo "^o^"
}


linkHiCProMatrix(){

        ## link hicpro matrix and bed file to this directory
        dir=$1
        res=$2
        if [[ -z $dir || -z $res ]];then
                echo
                echo "Usage: linkHiCProMatrix directory resolution"
                echo "      link hicpro matrix and bed file to this directory"
                return 1
        fi
        ln -s $dir/hic_results/matrix/*/iced/$res/*matrix .
        ln -s $dir/hic_results/matrix/*/raw/$res/*bed .
        ln -s $dir/hic_results/matrix/*/raw/$res/*matrix .
}


hicpro2cool1() {
        ## Convert hicpro to cool format use hicexplorer
        bedfile=$1
        matrix=$2

        if [[ -z $bedfile || -z $matrix ]]; then
                echo
                echo "Usage: hicpro2cool rice_5000_iced.matrix rice_5000_abs.bed"
                echo
                return 1
        fi
        outname=`basename ${matrix%%.matrix}`
        hicConvertFormat -m $matrix --bedFileHicpro $bedfile --inputFormat hicpro --outputFormat cool -o ${outname}.cool
}


 qssh() {
        
        ## qrlogin node

        node=$1
        if [[ -z $node ]];then
                echo
                echo "Usage: qssh node1"
                echo
                return 1
        fi
        
        qrsh -q all.q@$node
}


sra2fq() {
        ## use fastq-dump to convert sra to fastq
        sra_num=$1

        if [ -z $sra_num ]; then
                echo
                echo "Usage: sra2fq SRA000000"
                echo 
                return 1
        fi
        nohup fastq-dump --gzip --split-3 $sra_num.sra -O srr${sra_num##SRR} &

}

ascp_sra() {
        ## use ascp to download sra 
        sra_num=$1
        if [ -z $sra_num ]; then
                echo
                echo "Usage: ascp_sra SRA000000"
                echo
                return 1
        fi

        ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T 1 -l 200m anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${sra_num:0:6}/$sra_num/${sra_num}.sra
}


wget_sra() {
        ## prefetch sra through wget

        sra_num=$1
        if [ -z $sra_num ]; then
                echo
                echo "Usage: wget_sra SRA000000"
                echo
                return 1
        fi

        wget -cb ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/${sra_num:0:6}/$sra_num/${sra_num}.sra 
}

dump_fq() {
        ## fastq-dump sra
        sra_file=$1

        if [ -z $sra_file ]; then
                echo
                echo "Usage: dump_fq SRR000000"
                echo
                return 1
        fi

        fastq-dump --gzip --split-3 $sra_file -O srr${sra_file##SRR} &
}

sra_fastp() {
        ## qsub submit all sra fastq file to use fastp control quality
        for sra in *_1.fastq.gz; do s=${sra%%_1.fastq.gz};echo "fastp -i ${s}_1.fastq.gz -o ${s}_R1.fastq.gz -I ${s}_2.fastq.gz -O ${s}_R2.fastq.gz -j ${s}.json -h ${s}.html"> run_${s}_fastp.sh;qsub -pe mpi 4 run_${s}_fastp.sh; done
}

sra_fastp_PE() {
        ## qsub submit all PE sra fastq file to use fastp control quality
        for sra in  *.fastq.gz; do s=${sra%%.fastq.gz};echo "fastp -i ${s}.fastq.gz -o ${s}_fp_.fastq.gz -j ${s}.json -h ${s}.html"> run_${s}_fastp.sh;qsub -pe mpi 4 run_${s}_fastp.sh; done
}


add() {
        awk '{s+=$1} END {print s}'
}

pwc() {
        if [[ -z $@ ]];then
                echo
                echo -e "Usage: cat sample.txt | $0 -l  "
                echo -e "\tparallel wc"
                echo
                return;
        fi
        parallel --block=10M --pipe wc $@ 2>/dev/null | awk '{s+=$1} END {print s}'
}


clean_lastdb(){
        rm *.{tis,suf,ssp,sds,prj,des,bck}
}

