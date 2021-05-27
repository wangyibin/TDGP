##
## allhic pipeline
##          snakemake pipeline for executing ALLHiC diploid with correct
## Examples:
##          snakemake -s allhic_diploid_withCorrect_pipeline.smk \
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
## tag:
##         - R1
##         - R2

import os
import os.path as op
import glob
import sys


__version__ = "v2.0"


enzyme_repo = {"MBOI":"GATC", "HINDIII":"AAGCTT"}

defaults = {'plot': True}

FA = config['genome']
SAMPLE = config['sample'][0]
enzyme = config['enzyme'].upper()
enzyme_base = enzyme_repo[enzyme]
N = config["cluster_N"]
ncpus = config['ncpus']
bin_sizes = config['bin_sizes']
fq_suffix = config['fq_suffix']
tag1 = config['tag'][0]
tag2 = config['tag'][0]


SAMPLES = list(map(op.basename, glob.glob("data/*{}.{}".format(tag1, fq_suffix))))
SAMPLES = list(map(lambda x: x.replace("."+fq_suffix, "").replace("_"+tag1, ""), SAMPLES))

if config['plot']:
    rule all:
        input:
            # lambda wildcards: "process/{}.bwa_aln.REduced.paired_only.bam".format(config['sample'][0])
            #"allhic_result/groups.asm.fasta"
            whole_pdf = expand("allhic_result/pdf/{bin_size}_Whole_genome.pdf", bin_size=bin_sizes),
            chroms_pdf = expand("allhic_result/pdf/{bin_size}_{sample}.merged.counts_{enzyme_base}.{N}g{n}.pdf",
                bin_size=bin_sizes, sample=(SAMPLE,), enzyme_base=(enzyme_base,), N=(N,), n=range(1,N+1))
else:
    rule all:
        input:
            #lambda wildcards: "process/{}.bwa_aln.REduced.paired_only.bam".format(config['sample'][0]),
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

rule round1_bwa_mem:
    input:
        fabwt = lambda wildcards: "{}.bwt".format(FA),
        fq1 = expand("data/{{sample}}_{tag}.{fq_suffix}", tag=[tag1], fq_suffix=[fq_suffix]),
        fq2 = expand("data/{{sample}}_{tag}.{fq_suffix}", tag=[tag2], fq_suffix=[fq_suffix])
    output:
        temp("1.bwa_result/{sample}_1st.sorted.bam")
    log:
        "logs/{sample}_bwa_mem_1st.log"
    threads: 4
    shell:
        "bwa mem -t {threads} {FA} {input.fq1} {input.fq2} 2>{log} | "
        "samtools sort -@ {threads} > {output} 2>{log}"

rule round1_bam_merge:
    input:
        expand("1.bwa_result/{sample}_1st.sorted.bam", sample=SAMPLES)
    output:
        expand("1.bwa_result/{name}_1st.sorted.bam", name=config['sample'])
    log:
        expand("logs/{name}_1st.merged.log", name=config['sample'])
    threads: ncpus
    shell:
        "samtools merge -@ {threads} {output} {input} 2>{log}"

rule round1_bam_index:
    input:
        "1.bwa_result/{sample}_1st.sorted.bam"
    output:
        "1.bwa_result/{sample}_1st.sorted.bam.bai"
    log:
        "logs/{sample}_1st.sorted.index.log"
    threads: ncpus
    shell:
        "samtools index -@ {threads} {input} 2>{log}"

rule correct:
    input:
        fasta = FA,
        bam = "1.bwa_result/{sample}_1st.sorted.bam",
        bai = "1.bwa_result/{sample}_1st.sorted.bam.bai"
    output:
       "{sample}.HiCcorrected.fasta"
    log:
        "logs/{sample}.correct.log"
    threads: ncpus
    shell:
        "ALLHiC_corrector -m {input.bam} -r {input.fasta} -o {output} -t {threads} 2>&1 1>{log}"

rule corrected_fa_index:
    input:
        expand("{sample}.HiCcorrected.fasta", sample=(SAMPLE, ))
        
    output:
        expand("{sample}.HiCcorrected.fasta.bwt", sample=(SAMPLE, )),
        expand("{sample}.HiCcorrected.fasta.fai", sample=(SAMPLE, ))
    log:
        expand("logs/{sample}_index.log", sample=(SAMPLE, ))
    shell:
        "samtools faidx {input} & bwa index -a bwtsw {input} 2>{log}"

rule round2_bwa_mem:
    input:
        fabwt = lambda wildcards: "{}.HiCcorrected.fasta.bwt".format(SAMPLE),
        fq1 = expand("data/{{sample}}_{tag}.{fq_suffix}", tag=[tag1], fq_suffix=[fq_suffix]),
        fq2 = expand("data/{{sample}}_{tag}.{fq_suffix}", tag=[tag2], fq_suffix=[fq_suffix])
    output:
        "2.bwa_result/{sample}.sorted.bam"
    log:
        "logs/{sample}_bwa_mem_2st.log"
    threads: 4
    shell:
        "bwa mem -t {threads} {FA} {input.fq1} {input.fq2} 2>{log} | "
        "samtools sort -@ {threads} - > {output} 2>{log}"

rule round2_bam_merge:
    input:
        expand("2.bwa_result/{sample}.sorted.bam", sample=SAMPLES)
    output:
        temp(expand("allhic_result/{name}.merged.bam", name=config['sample']))
    log:
        expand("logs/{name}_2st.merged.log", name=config['sample'])
    threads: ncpus
    shell:
        "samtools merge -@ {threads} {output} {input} 2>{log}"

rule preprocess_bam:
    input:
        bam = expand("allhic_result/{name}.merged.bam", name=SAMPLE),
        fasta = lambda wildcards: "{}.HiCcorrected.fasta".format(SAMPLE)
    output:
        expand("allhic_result/{name}.merged.REduced.paired_only.bam", name=SAMPLE)
    log:
        lambda wildcards: "logs/{}.preprocess.log".format(SAMPLE)
    params:
       enzyme = config["enzyme"]
    threads: ncpus
    shell:
        "PreprocessSAMs.pl {input.bam} {input.fasta} {params.enzyme} {threads} 2>&1 1>{log}"

rule partition:
    input:
        bam =  expand("allhic_result/{name}.merged.REduced.paired_only.bam", name=SAMPLE),
        fasta = lambda wildcards: "{}.HiCcorrected.fasta".format(SAMPLE)
    output:
        expand("allhic_result/{sample}.merged.REduced.paired_only.counts_{enzyme_base}.{N}g{n}.txt",
            sample=(SAMPLE,), enzyme_base=(enzyme_base,), N=(N,), n=range(1,N+1)),
        expand("allhic_result/{sample}..REduced.paired_only.clm", sample=(SAMPLE, ))
    log:
        "logs/{}.partition.log".format(SAMPLE)
    shell:
        "ALLHiC_partition -b {input.bam} -r {input.fasta}  -e {enzyme_base} -k {N} 1>{log} 2>&1"

# rule extract:
#     input:
#         bam = "allhic_result/{sample}.merged.bam",
#         #txt = "allhic_result/{sample}.clusters.txt"
#         txt = expand("allhic_result/{sample}.merged.counts_{enzyme_base}.{N}g{n}.txt",
#             sample=(SAMPLE,), enzyme_base=(enzyme_base,), N=(N,), n=range(1,N+1))

#     output:
#         "allhic_result/{sample}.merged.clm",
#     log:
#         "logs/{sample}_extract.log"
#     shell:
#         "allhic extract {input.bam} {FA} --RE {enzyme_base} 2>{log}"

rule optimize:
    input:
        clm = "allhic_result/{sample}.merged.clm",
        txt = "allhic_result/{sample}.merged.counts_{enzyme_base}.{N}g{n}.txt"
    output:
        "allhic_result/{sample}.merged.counts_{enzyme_base}.{N}g{n}.tour"
    log:
        "logs/{sample}.{enzyme_base}.{N}g{n}_optimize.log"
    threads: ncpus
    shell:
        "allhic optimize {input.txt} {input.clm} 2>{log}"

rule build:
    input:
        expand("allhic_result/{sample}.merged.counts_{enzyme_base}.{N}g{n}.tour",
            sample=(SAMPLE,), enzyme_base=(enzyme_base,), N=(N,), n=range(1,N+1))
    output:
        fasta = "allhic_result/groups.asm.fasta",
        agp = "allhic_result/groups.agp"
    log:
        "logs/allhic.build.log"
    shell:
        "cd allhic_result/ && ALLHiC_build ../{FA} 2 > ../{log}"


rule getSizes:
    input:
        "allhic_result/groups.asm.fasta"
    output:
        "allhic_result/chrom.sizes"
    shell:
        "getChrLength.py {input} | grep -v tig > {output}"

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
            chroms_pdf = expand("allhic_result/pdf/{bin_size}_{sample}.merged.counts_{enzyme_base}.{N}g{n}.pdf",
                bin_size=bin_sizes, sample=(SAMPLE,), enzyme_base=(enzyme_base,), N=(N,), n=range(1,N+1))
        params:
            ",".join(bin_sizes)
        threads: ncpus
        shell:
            "ALLHiC_plot3 {input.bam} {input.agp} {input.sizes} -t {threads} --bin_size {params} && mv *.pdf allhic_result/pdf"

        