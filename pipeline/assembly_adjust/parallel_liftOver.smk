
import glob
import os
import os.path as op
import sys

script_path = "~/code/TDGP/pipeline/assembly_adjust"

SAMPLE = config['sample']
SAMPLES = list(map(op.basename, glob.glob('data/*')))
ncpus = config['ncpus']
old_agp = config['old_agp']
new_agp = config['new_agp']
try:
    os.makedirs("tmp_smk")
except:
    pass

rule all:
    input:
        "out.merged.ValidPairs"

rule liftOver:
    input:
        "data/{sample}"
    output:
        temp("split_liftover_results/{sample}.new.validpairs")
    log: "logs/{sample}.liftover.logs"
    threads: ncpus
    params:
        chunksize = 1000000,
        tmp = "tmp_smk/{sample}"
    shell:
        "python {script_path}/liftOverValidPairs.py "
                "{input} {old_agp} {new_agp} -t {threads} -o {output} "
                "-c {params.chunksize} -T {params.tmp}/ 2>{log}"
        

rule merge:
    input:
        expand("split_liftover_results/{sample}.new.validpairs", sample=SAMPLES)
    output:
        "out.merged.ValidPairs"
    shell:
        "cat {input} > {output}"

