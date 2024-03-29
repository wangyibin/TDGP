##
## allhic pipeline
##          snakemake pipeline for executing ALLHiC diploid
## Examples:
##          snakemake -s allhic_diploid_pipeline.smk \
##              --configfile config_allhic_diploid_pipeline.yaml \
##              -j 12 \
##              --cluster "qsub -l select=1:ncpus={threads} -q workq -j oe"
##---------------------------------------------------------------------------
## 
## genome:
##         tig.fasta
##
## sample:
##         - AD
##
## ncpus:
##         20
##
## enzyme:
##         HINDIII
##
## cluster_N:
##         1
## Allele:
##          Allele.ctg.table
## tag:
##         - R1
##         - R2

import os
import os.path as op
import glob
import sys


__version__ = "v2.0"


enzyme_repo = {"MBOI":"GATC", "HINDIII":"AAGCTT"}


FA = config['genome']
SAMPLE = config['sample'][0]
enzyme = config['enzyme'].upper()
enzyme_base = enzyme_repo[enzyme]
N = config["cluster_N"]
ncpus = config['ncpus']
bin_sizes = config['bin_sizes']
fq_suffix = config['fq_suffix']
tag1 = config['tag'][0]

Allele = config['Allele']

SAMPLES = list(map(op.basename, glob.glob("data/*{}.{}".format(tag1, fq_suffix))))
SAMPLES = list(map(lambda x: x.replace("."+fq_suffix, "").replace("_"+tag1, ""), SAMPLES))

if config['plot']:
    rule all:
        input:
            # lambda wildcards: "process/{}.bwa_aln.REduced.paired_only.bam".format(config['sample'][0])
            #"allhic_result/groups.asm.fasta"
            whole_pdf = expand("allhic_result/pdf/{bin_size}_Whole_genome.pdf", bin_size=bin_sizes),
            chroms_pdf = expand("allhic_result/pdf/{bin_size}_group{n}.pdf", bin_size=bin_sizes, n=range(1,N+1))
            
else:
    rule all:
        input:
            # lambda wildcards: "process/{}.bwa_aln.REduced.paired_only.bam".format(config['sample'][0])
            "allhic_result/groups.asm.fasta"

rule fa_index:
    input:
        FA
    output:
        "{input}.bwt",
        "{input}.fai"
    log:
        "logs/{input}.log"
    shell:
        "samtools faidx {input} & bwa index -a bwtsw {input} 2>{log}"


rule bwa_aln:
    input:
        fabwt = lambda wildcards: "{}.bwt".format(FA),
        fq = expand("data/{{sample}}.{fq_suffix}", fq_suffix=[fq_suffix])
    output:
        temp("bwa_result/{sample}.sai")
    log:
        "logs/{sample}.log"
    threads: 4
    shell:
        "bwa aln -t {threads} {FA} {input.fq} > {output} 2>{log}"


rule bwa_sampe:
    input:
        fq = expand("data/{{sample}}_{tag}.fastq.gz", tag=config["tag"]),
        sai = expand("bwa_result/{{sample}}_{tag}.sai", tag=config["tag"]),
        
    output:
        temp("bwa_result/{sample}.bwa_sampe.bam")
    log:
        "logs/{sample}.bwa_sampe.log"
    shell:
        "bwa sampe {FA} {input.sai} {input.fq} 2>{log}| samtools view -bhS -> {output}"

rule preprocess_bam:
    input:
        bam = "bwa_result/{sample}.bwa_sampe.bam"
    output:
        temp("bwa_result/{sample}.bwa_sampe.REduced.paired_only.bam")
    log:
        "logs/{sample}.preprocess.log"
    params:
       enzyme = config["enzyme"]
    threads: 4
    shell:
        "PreprocessSAMs.pl {input} {FA} {params.enzyme} 4 2>&1 1>{log}"


rule filter_bam:
    input:
        "bwa_result/{sample}.bwa_sampe.REduced.paired_only.bam"
    output:
        temp("bwa_result/{sample}.clean.sam")
    log:
        "logs/{sample}.filter.log"
    shell:
        "filterBAM_forHiC.pl {input} {output} 1>{log} 2>&1"

rule sam2bam:
    input:
        "bwa_result/{sample}.clean.sam"
    output:
        temp("bwa_result/{sample}.clean.bam")
    log:
        "logs/{sample}.sam2bam.log"
    threads: ncpus
    shell:
        "samtools view -@ {threads} -bt {FA}.fai {input} > {output} 2>{log}"

rule merge_bam:
    input:
        bam = expand("bwa_result/{sample}.clean.bam", sample=SAMPLES)
    output:
        "allhic_result/{SAMPLE}.merged.bam"
    log:
        "logs/{SAMPLE}.merged.log"
    threads: ncpus
    shell:
        "samtools merge -@ {threads} {output} {input.bam} 2>{log}"

rule prune:
    input:
        lambda wildcards: "allhic_result/{}.merged.bam".format(SAMPLE)
    output:
        "prunning.bam"
    log:
        "logs/{}.prunning.log".format(SAMPLE)
    threads: ncpus
    shell:
        "ALLHiC_prune -i {Allele} -b {input} -r {FA} 2>{log} "

rule partition:
    input:
        "prunning.bam"
    output:
        #expand("allhic_result/{sample}.pruning.counts_{enzyme_base}.{N}g{n}.txt",
        #    sample=(SAMPLE,), enzyme_base=(enzyme_base,), N=(N,), n=range(1,N+1)),
        "prunning.clusters.txt",
        "prunning.counts_{}.txt".format(enzyme_base),
        # "allhic_result/prunning.clm"
    log:
        "logs/{}.partition.log".format(SAMPLE)
    shell:
        "ALLHiC_partition -b {input} -r {FA} -e {enzyme_base} -k {N} 1>{log} 2>&1 "
        

rule rescue:
    input:
        bam = lambda wildcards: "allhic_result/{}.merged.bam".format(SAMPLE),
        cluster = "prunning.clusters.txt",
        counts = lambda wildcards: "prunning.counts_{}.txt".format(enzyme_base),
        # clm = "allhic_result/prunning.clm"
    output:
        expand("allhic_result/group{n}.txt", n=range(1, N+1)),
    log:
        "logs/{}.rescue.log".format(SAMPLE)
    threads: ncpus
    shell:
        "ALLHiC_rescue -b {input.bam} -r {FA} -c {input.cluster} -i {input.counts} 2>{log} && "
        "mv group*txt prunning* allhic_result"
   
rule extract:
    input:
        bam = lambda wildcards: "allhic_result/{}.merged.bam".format(SAMPLE),
        #txt = "allhic_result/{sample}.clusters.txt"
        
    output:
        "allhic_result/{SAMPLE}.merged.clm",
    log:
        "logs/{SAMPLE}_extract.log"
    shell:
        "allhic extract {input.bam} {FA} --RE {enzyme_base} 2>{log}"

rule optimize:
    input:
        clm = lambda wildcards: "allhic_result/{}.merged.clm".format(SAMPLE),
        txt = "allhic_result/group{n}.txt"
    output:
        "allhic_result/group{n}.tour"
    log:
        "logs/group{n}_optimize.log"
    threads: 4
    shell:
        "allhic optimize {input.txt} {input.clm} 1>{log} 2>&1"

rule build:
    input:
        expand("allhic_result/group{n}.tour", n=range(1,N+1))
    output:
        fasta = "allhic_result/groups.asm.fasta",
        agp = "allhic_result/groups.agp"
    log:
        "logs/allhic.build.log"
    shell:
        "cd allhic_result && ALLHiC_build ../{FA} 2 > ../{log}"

rule getSizes:
    input:
        "allhic_result/groups.asm.fasta"
    output:
        "allhic_result/chrom.sizes"
    shell:
        "getChrLength.py {input} | grep -v tig | grep group > {output}"

rule sort_bam:
    input:
        expand("allhic_result/{name}.merged.bam", name=config['sample'])
    output:
        expand("allhic_result/{name}.sorted.bam", name=config['sample'])
    log:
        expand("logs/{name}.sort_bam.log", name=config['sample'])
    threads: ncpus
    shell:
        "samtools sort -@ {threads} {input} -o {output} 2>{log}"

rule bam_index:
    input:
        expand("allhic_result/{name}.sorted.bam", name=config['sample'])
    output:
        expand("allhic_result/{name}.sorted.bam.bai", name=config['sample'])
    log:
        expand("logs/{name}.bam_index.log", name=config['sample'])
    threads: ncpus
    shell:
        "samtools index -@ {threads} {input}  2>{log}"

if config['plot']:
    rule plot:
        input:
            bam = expand("allhic_result/{sample}.sorted.bam", sample=(SAMPLE,)),
            bai = expand("allhic_result/{sample}.sorted.bam.bai", sample=(SAMPLE, )),
            agp = "allhic_result/groups.agp",
            sizes = "allhic_result/chrom.sizes"
        output:
            whole_pdf = expand("allhic_result/pdf/{bin_size}_Whole_genome.pdf", bin_size=bin_sizes),
            chroms_pdf = expand("allhic_result/pdf/{bin_size}_group{n}.pdf", bin_size=bin_sizes, n=range(1,N+1))
        params:
            ",".join(bin_sizes)
        threads: ncpus
        shell:
            "ALLHiC_plot3 {input.bam} {input.agp} {input.sizes} -t {threads} --bin_size {params} && mv *.pdf allhic_result/pdf"

    