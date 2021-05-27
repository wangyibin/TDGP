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

Allele = config['Allele']
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
        fasta = lambda wildcards: "{}.HiCcorrected.fasta".format(SAMPLE),
        fabwt = lambda wildcards: "{}.HiCcorrected.fasta.bwt".format(SAMPLE),
        fq1 = expand("data/{{sample}}_{tag}.{fq_suffix}", tag=[tag1], fq_suffix=[fq_suffix]),
        fq2 = expand("data/{{sample}}_{tag}.{fq_suffix}", tag=[tag2], fq_suffix=[fq_suffix])
    output:
        "2.bwa_result/{sample}.sorted.bam"
    log:
        "logs/{sample}_bwa_mem_2st.log"
    threads: 4
    shell:
        "bwa mem -t {threads} {input.fasta} {input.fq1} {input.fq2} 2>{log} | "
        "samtools sort -@ {threads} - > {output} 2>{log}"


rule make_bed_around_RE_site:
    input:
        expand("{sample}.HiCcorrected.fasta", sample=(SAMPLE, ))
    output:
        f"{SAMPLE}.HiCcorrected.fasta.near_{enzyme_base}.500.bed",
        f"{SAMPLE}.HiCcorrected.fasta.pos_of_{enzyme_base}.txt",
    log: lambda wildcards: "logs/{}.make_re_site.log".format(SAMPLE)
    params:
       enzyme = enzyme_base
    shell:
        "make_bed_around_RE_site.pl {input} {params.enzyme} 500 2>{log}"

rule round2_bam_preprocess:
    input:
        lambda wildcards: "{}.HiCcorrected.fasta.pos_of_{}.txt".format(SAMPLE, enzyme_base),
        bed_re_file = lambda wildcards: "{}.HiCcorrected.fasta.near_{}.500.bed".format(SAMPLE, enzyme_base),
        bam = "2.bwa_result/{sample}.sorted.bam"
    output:
        reduced_bam = temp("2.bwa_result/{sample}.reduced.bam"),
        pair_only_bam = temp("2.bwa_result/{sample}.reduced.paired_only.bam")
    log: "logs/{sample}.2st.bam.preprocess.log"
    threads: 4
    shell:
        "bedtools intersect -abam {input.bam} -b {input.bed_re_file} > {output.reduced_bam} 2>{log} &&"
        "samtools view -@ {threads} -F12 {output.reduced_bam} -b -o {output.pair_only_bam} 2>{log}"   

rule round2_bam_merge:
    input:
        expand("2.bwa_result/{sample}.reduced.paired_only.bam", sample=SAMPLES)
    output:
        temp(expand("allhic_result/{name}.clean.bam", name=config['sample']))
    log:
        expand("logs/{name}_2st.merged.log", name=config['sample'])
    threads: ncpus
    shell:
        "samtools merge -@ {threads} {output} {input} 2>{log}"


rule prune:
    input: 
        bam =   expand("allhic_result/{name}.clean.bam", name=SAMPLE),
        fasta = lambda wildcards: "{}.HiCcorrected.fasta".format(SAMPLE)
    output:
        "prunning.bam"
    log:
        "logs/{}.prunning.log".format(SAMPLE)
    threads: ncpus
    shell:
        "ALLHiC_prune -i {Allele} -b {input.bam} -r {input.fasta} 2>{log}"

rule partition:
    input:
        bam =  "prunning.bam",
        fasta = lambda wildcards: "{}.HiCcorrected.fasta".format(SAMPLE)
    output:
        "prunning.clusters.txt",
        "prunning.counts_{}.txt".format(enzyme_base)
    log:
        "logs/{}.partition.log".format(SAMPLE)
    shell:
        "ALLHiC_partition -b {input.bam} -r {input.fasta}  -e {enzyme_base} -k {N} 1>{log} 2>&1"

rule rescue:
    input:
        fasta = lambda wildcards: "{}.HiCcorrected.fasta".format(SAMPLE),
        bam =  expand("allhic_result/{name}.clean.bam", name=SAMPLE),
        cluster = "prunning.clusters.txt",
        counts = lambda wildcards: "prunning.counts_{}.txt".format(enzyme_base),
        # clm = "allhic_result/prunning.clm"
    output:
        expand("allhic_result/group{n}.txt", n=range(1, N+1)),
    log:
        "logs/{}.rescue.log".format(SAMPLE)
    threads: ncpus
    shell:
        "ALLHiC_rescue -b {input.bam} -r {input.fasta} -c {input.cluster} -i {input.counts} 2>{log} && "
        "mv group*txt prunning* allhic_result"

rule extract:
    input:
        bam =  expand("allhic_result/{name}.clean.bam", name=(SAMPLE, )),
        fasta = lambda wildcards: "{}.HiCcorrected.fasta".format(SAMPLE),
        
    output:
        "allhic_result/{SAMPLE}.clean.clm",
    log:
        "logs/{SAMPLE}_extract.log"
    shell:
        "allhic extract {input.bam} {input.fasta} --RE {enzyme_base} 2>{log}"

rule optimize:
    input:
        clm = lambda wildcards: "allhic_result/{}.clean.clm".format(SAMPLE),
        txt = "allhic_result/group{n}.txt"
    output:
        "allhic_result/group{n}.tour"
    log:
        "logs/group{n}_optimize.log"
    threads: ncpus
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
        "getChrLength.py {input} | grep -v tig > {output}"

rule sort_bam:
    input:
        expand("allhic_result/{name}.clean.bam", name=config['sample'])
    output:
        expand("allhic_result/{name}.sorted.bam", name=config['sample'])
    log:
        expand("logs/{name}.sort_bam.log", name=config['sample'])
    threads: ncpus
    shell:
        "samtools sort -@ {threads} {input} -o {output} 2>{log}"

rule bam_index:
    input:
        expand("allhic_result/{name}.sorted.bam", name=SAMPLE),
    output:
        expand("allhic_result/{name}.sorted.bam.bai", name=SAMPLE),
    log:
        expand("logs/{name}.bam_index.log", name=config['sample'])
    threads: ncpus
    shell:
        "samtools index -@ {threads} {input}  2>{log}"

if config['plot']:
    rule plot:
        input:
            bam =  expand("allhic_result/{name}.sorted.bam", name=SAMPLE),
            bai = expand("allhic_result/{name}.sorted.bam.bai", name=SAMPLE),
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

        