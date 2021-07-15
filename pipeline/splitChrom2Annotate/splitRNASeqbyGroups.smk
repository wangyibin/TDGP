## 
## split RNA-Seq data by groups


import os
import os.path as op
import sys
import glob

from collections import OrderedDict


groups_db_file = config['groups']
groups_db = OrderedDict()
with open(groups_db_file) as fp:
    for line in fp:
        group, chrom = line.strip().split()
        if group not in groups_db:
            groups_db[group] = []
        groups_db[group].append(chrom)
groups_list = list(groups_db)
genome = config['genome']
ncpus = config['ncpus']
ext = config['ext']
R1, R2 = ext
suffix = config['suffix']
hisat2_k = 4

if not op.exists('data/'):
    print("[Error]: No such file of data")
    sys.exit()
SAMPLES = list(map(op.basename, glob.glob("data/*{}.{}".format(R1, suffix))))
SAMPLES = list(map(lambda x: x.replace("."+suffix, "").replace("_"+R1, ""), SAMPLES))

rule all:
    input:
        expand("{group}/{group}_{ext}.fastq.gz", group=groups_list, ext=ext),
 

rule hisat2build:
    input:
        genome
    output:
        expand("{genome}.{idx}.ht2l", genome=(genome, ), idx=range(1,9))
    log:
        f"logs/{genome}.hisat2-build.log"
    threads: ncpus
    shell:
        "hisat2-build -p {threads} {input} {input} 1>{log} 2>&1"

rule getBed:
    input:
        "{group}/{group}.fasta"
    output:
        "{group}/{group}.bed"
    shell:
        "samtools faidx {input} && "
        "awk '{{print $1,0,$2}}' OFS='\t' {input}.fai > {output}"

rule hisat2:
    input:
        index = expand(f"{genome}.{{idx}}.ht2l", idx=range(1,9)),
        left = f"data/{{sample}}_{R1}.{suffix}",
        right = f"data/{{sample}}_{R2}.{suffix}"
    output:
        protected("hisat2_result/{sample}.sorted.bam")
    log:
        "logs/{sample}.hisat2.log"
    threads: ncpus 
    params:
        k = hisat2_k
    shell:
        "hisat2 -p {threads} -k {params.k} -x {genome} -1 {input.left} "
        "-2 {input.right} 2>{log} | samtools sort -T hisat2_result/ -@ 4 | " 
        "samtools view -bhS > {output}"

rule extractFastqFromBam:
    input:
        fq = f"data/{{sample}}_{{ext}}.{suffix}",
        bam = "hisat2_result/{sample}.sorted.bam",
        bed = "{group}/{group}.bed"
    output:
        temp("{group}/{sample}.{group}_{ext}.fastq")
    threads: 4
    shell:
        "samtools view -@ {threads} {input.bam} -L {input.bed} | "
        "parallel --pipe -k cut -f 1 | "
        "seqkit grep -f - {input.fq} > {output}"
    
rule mergeFastq:
    input:
        expand("{{group}}/{sample}.{{group}}_{{ext}}.fastq", sample=SAMPLES)
    output:
        "{group}/{group}_{ext}.fastq.gz"
    threads: 4
    shell:
        "cat {input} | pigz -p {threads} > {output}"


