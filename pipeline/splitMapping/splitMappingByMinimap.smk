import os
import os.path as op
import glob
import sys


__version__ = "v1.0"



FA = config['genome']
SAMPLE = config['sample'][0]
ncpus = config['ncpus']
fq_suffix = config['fq_suffix']

SAMPLES = list(map(op.basename, glob.glob("data/*.{}".format(fq_suffix))))
SAMPLES = list(map(lambda x: x.replace("."+fq_suffix, ""), SAMPLES))

rule all:
    input:
        expand("{sample}.sorted.bam", sample=(SAMPLE,))

rule fa_index:
    input:
        FA
    output:
        "{input}.fai"
    log:
        "logs/{input}.log"
    shell:
        "samtools faidx {input} 2>{log}"

rule minimap2:
    input:
        fa = FA,
        fai = f'{FA}.fai',
        fq = expand("data/{{sample}}.{fq_suffix}", fq_suffix=[fq_suffix])
    output:
        temp("minimap2_result/{sample}.sorted.bam")
    log:
        "logs/{sample}_minimap2.log"
    params: mm2_param = config['mm2_param']
    threads: ncpus
    shell:
        "minimap2 -t {threads} -a {params.mm2_param} {input.fa} {input.fq} | "
        "samtools view -@ {threads} -bt {input.fai} | "
        "samtools sort -@ {threads} > {output} 2>{log}"

rule bam_merge:
    input:
        expand("minimap2_result/{sample}.sorted.bam", sample=SAMPLES)
    output:
        expand("{name}.sorted.bam", name=(SAMPLE,))
    log:
        expand("logs/{name}.merged.log", name=config['sample'])
    threads: ncpus
    shell:
        "samtools merge -@ {threads} {output} {input} 2>{log}"