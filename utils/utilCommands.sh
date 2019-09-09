#!/bin/bash



## my usualy use command library


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
        fastq-dump --gzip --split-3 $sra_num.sra -O srr${sra_num##SRR} &

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
