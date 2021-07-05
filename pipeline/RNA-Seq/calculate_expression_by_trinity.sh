#!/bin/bash

##################################
#             configure          #
##################################


transcripts="../data/F5.cds"
dataPath="../data"
sampleList="F5.list"
R1="_R1.fq.gz"
R2="_R2.fq.gz"
thread=4


for sample in `cat ${sampleList}`; do
left="${dataPath}/${sample}${R1}"
right="${dataPath}/${sample}${R2}"
cmd=" 
align_and_estimate_abundance.pl --transcripts ${transcripts} \
--seqType fq --left ${left} \
--right ${right} --est_method RSEM \
--output_dir ${sample} --aln_method bowtie \
--prep_reference --thread_count ${thread} 
"
clusterHeader -t $thread -s "$cmd" > run_align_$sample.sh
done 
