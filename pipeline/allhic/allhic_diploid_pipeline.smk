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
## tag:
##         - R1
##         - R2



__version__ = "v1.0"


enzyme_repo = {"MBOI":"GATC","HINDIII":"AAGCTT"}


FA = config['genome']
SAMPLE = config['sample'][0]
enzyme = config['enzyme'].upper()
enzyme_base = enzyme_repo[enzyme]
N = config["cluster_N"]
ncpus = config['ncpus']

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
        fq = "data/{sample}_{tag}.fastq.gz"
    output:
        "bwa_result/{sample}_{tag}.sai"
    log:
        "logs/{sample}_{tag}.log"
    threads: ncpus
    shell:
        "bwa aln -t {threads} {FA} {input.fq} > {output} 2>{log}"


rule bwa_sampe:
    input:
        fq = expand("data/{sample}_{tag}.fq.gz", sample=config['sample'], tag=config["tag"]),
        sai = expand("bwa_result/{sample}_{tag}.sai",sample=config["sample"],tag=config["tag"])
    output:
        "bwa_result/{sample}.bwa_sampe.bam"
    log:
        "logs/{sample}.bwa_sampe.log"
    shell:
        "bwa sampe {FA} {input.sai} {input.fq} | samtools view -F 4 -bS > {output}"
 


rule preprocess_bam:
    input:
        bam = "bwa_result/{sample}.bwa_sampe.bam"
    output:
        "bwa_result/{sample}.bwa_aln.REduced.paired_only.bam"
    log:
        "logs/{sample}_preprocess.log"
    params:
       enzyme = config["enzyme"]
    shell:
        "PreprocessSAMs.pl {input} {FA} {params.enzyme}"


rule filter_bam:
    input:
        "bwa_result/{sample}.bwa_aln.REduced.paired_only.bam"
    output:
        "allhic_result/{sample}.clean.sam"
    log:
        "logs/{sample}_filter.log"
    shell:
        "filterBAM_forHiC.pl {input} {output} 2>{log}"


rule sam2bam:
    input:
        "allhic_result/{sample}.clean.sam"
    output:
        "allhic_result/{sample}.clean.bam"
    threads: ncpus
    shell:
        "samtools view -@ {threads} -bt {FA}.fai {input} > {output}"

rule partition:
    input:
        lambda wildcards: "allhic_result/{}.clean.bam".format(SAMPLE)
    output:
        expand("allhic_result/{sample}.clean.counts_{enzyme_base}.{N}g{n}.txt",
            sample=(SAMPLE,), enzyme_base=(enzyme_base,), N=(N,), n=range(1,N+1))
    #log:
    #    lambda wildcards: "logs/{}_partition.log".format(SAMPLE)
    shell:
        "ALLHiC_partition -b {input} -r {FA} -e {enzyme_base} -k {N} 2>{log}"

rule extract:
    input:
        bam = "allhic_result/{sample}.clean.bam",
        #txt = "allhic_result/{sample}.clusters.txt"
        txt = expand("allhic_result/{sample}.clean.counts_{enzyme_base}.{N}g{n}.txt",
            sample=(SAMPLE,), enzyme_base=(enzyme_base,), N=(N,), n=range(1,N+1))

    output:
        "allhic_result/{sample}.clean.clm",
    log:
        "logs/{sample}_extract.log"
    shell:
        "allhic extract {input.bam} {FA} --RE {enzyme_base} 2>{log}"

rule optimize:
    input:
        clm = "allhic_result/{sample}.clean.clm",
        txt = "allhic_result/{sample}.clean.counts_{enzyme_base}.{N}g{n}.txt"
    output:
        "allhic_result/{sample}.clean.counts_{enzyme_base}.{N}g{n}.tour"
    log:
        "logs/{sample}_{enzyme_base}.{N}g{n}_optimize.log"
    shell:
        "allhic optimize {input.txt} {input.clm} 2>{log}"

rule build:
    input:
        expand("allhic_result/{sample}.clean.counts_{enzyme_base}.{N}g{n}.tour",
            sample=(SAMPLE,), enzyme_base=(enzyme_base,), N=(N,), n=range(1,N+1))
    output:
        "allhic_result/groups.asm.fasta"
    log:
        "logs/allhic_build.log"
    shell:
        "cd allhic_result/ && ALLHiC_build ../{FA} 2 > ../{log}"

