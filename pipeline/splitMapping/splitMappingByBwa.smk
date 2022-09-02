import os
import os.path as op
import glob
import sys


__version__ = "v1.0"



FA = config['genome']
SAMPLE = config['sample'][0]
ncpus = config['ncpus']
fq_suffix = config['fq_suffix']
tag1 = config['tag'][0]
tag2 = config['tag'][0]

SAMPLES = list(map(op.basename, glob.glob("data/*{}.{}".format(tag1, fq_suffix))))
SAMPLES = list(map(lambda x: x.replace("."+fq_suffix, "").replace("_"+tag1, ""), SAMPLES))

rule all:
    input:
        expand("{sample}.sorted.bam", sample=(SAMPLE,))

rule fa_index:
    input:
        FA
    output:
        "{input}.bwt",
        "{input}.fai"
    log:
        "logs/{input}.log"
    shell:
        "samtools faidx {input} & bwa index {input} 2>{log}"

rule bwa_mem:
    input:
        fabwt = lambda wildcards: "{}.bwt".format(FA),
        fq1 = expand("data/{{sample}}_{tag}.{fq_suffix}", tag=[tag1], fq_suffix=[fq_suffix]),
        fq2 = expand("data/{{sample}}_{tag}.{fq_suffix}", tag=[tag2], fq_suffix=[fq_suffix])
    output:
        temp("bwa_result/{sample}.sorted.bam")
    log:
        "logs/{sample}_bwa_mem.log"
    threads: ncpus
    shell:
        "bwa mem -t {threads} {FA} {input.fq1} {input.fq2} | "
        "samtools sort -@ {threads} > {output} 2>{log}"

rule bam_merge:
    input:
        expand("bwa_result/{sample}.sorted.bam", sample=SAMPLES)
    output:
        expand("{name}.sorted.bam", name=(SAMPLE,))
    log:
        expand("logs/{name}.merged.log", name=config['sample'])
    threads: ncpus
    shell:
        "samtools merge -@ {threads} {output} {input} 2>{log}"