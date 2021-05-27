import os.path as op


SAMPLES = config['samples']
memory = 10
try:
    ncpus = config['ncpus']
except:
    ncpus = 20

try:
    kmers = config['kmers']
except:
    kmers = [23, 25, 27]

sample1, sample2, sample3 = SAMPLES
KMC_SUFFIX = ('kmc_suf', 'kmc_pre')
CONBINES = ('100', '110', '101', '010',
            '011', '001', '111')
EXPRESSION = {'100': f'{sample1} - {sample2} - {sample3}', 
               '110': f'{sample1} * {sample2} - {sample3}', 
               '101': f'{sample1} * {sample3} - {sample2}',
               '010': f'{sample2} - {sample1} - {sample3}', 
               '011': f'{sample2} * {sample3} - {sample1}', 
               '001': f'{sample3} - {sample1} - {sample2}', 
               '111': f'{sample1} * {sample2} * {sample3}'
}


rule all:
    input:
        expand(f"{sample1}-{sample2}-{sample3}.k{{k}}.venn.{{fmt}}",
            k=kmers, fmt=('png', 'pdf', 'svg'))

rule kmer_count:
    input:
        expand("data/{{sample}}_{suffix}" ,suffix=("R1.fastq.gz", "R2.fastq.gz"))
    output:
        "{sample}.k{k}.data.list",
        expand("k{{k}}/{{sample}}.k{{k}}.{suffix}", suffix=KMC_SUFFIX) 
    threads: ncpus
    log: "logs/{sample}.k{k}.kmc.log"
    resources: 
        mem_mb = 10000000
    params:
        mem = 10000000/1000000
    shell:
        "mkdir -p tmp_{wildcards.sample}_{wildcards.k} && "
        "echo -e '{input[0]}\\n{input[1]}' > {output[0]} &&"
        "kmc -k{wildcards.k} -m{params.mem} -t{threads} "
        "@{output[0]} k{wildcards.k}/{wildcards.sample}.k{wildcards.k} "
        "./tmp_{wildcards.sample}_{wildcards.k} 2>{log} && "
        "rm -rf tmp_{wildcards.sample}_{wildcards.k}"

rule configure_complex:
    output:
        expand('k{{k}}/{conbine}.k{{k}}.config', conbine=CONBINES)
    run:
        s1 = f"k{wildcards.k}/{sample1}.k{wildcards.k}"
        s2 = f"k{wildcards.k}/{sample2}.k{wildcards.k}"
        s3 = f"k{wildcards.k}/{sample3}.k{wildcards.k}"
        for conbine, expression in EXPRESSION.items():
            
            CONFIGSTRING = f"""INPUT:
{sample1}={s1}
{sample2}={s2}
{sample3}={s3}
OUTPUT:
k{wildcards.k}/{conbine}={expression}
"""
            with open(f'k{wildcards.k}/{conbine}.k{wildcards.k}.config', 'w') as out:
                out.write(CONFIGSTRING)

rule kmc_complex:
    input:
        'k{k}/{conbine}.k{k}.config',
        expand("k{{k}}/{sample}.k{{k}}.{suffix}", sample=SAMPLES, suffix=KMC_SUFFIX)
    output:
        expand("k{{k}}/{{conbine}}.{suffix}", suffix=KMC_SUFFIX)
    resources:
        mem_mb = 20000000
    log: "logs/{conbine}.k{k}.complex.log"
    shell:
        "kmc_tools complex {input[0]} 2>{log}"

rule kmc_dump:
    input:
        expand("k{{k}}/{{conbine}}.{suffix}", suffix=KMC_SUFFIX)
    output:
        "k{k}/{conbine}.k{k}.res"
    shell:
        "kmc_dump k{wildcards.k}/{wildcards.conbine} /dev/stdout | wc -l > {output}"

rule concat_res:
    input:
        expand("k{{k}}/{conbine}.k{{k}}.res", conbine=CONBINES)
    output:
        f"k{{k}}/{sample1}-{sample2}-{sample3}.k{{k}}.res"
    run:
        out = open(output[0], 'w')

        for sample_res in input:
            s = op.basename(sample_res).rsplit(".")[0]
            with open(sample_res, 'r') as fp:
                number = int(fp.read().strip())
                out.write(f"{s}\t{number}\n")
        out.close()
        
rule plot_venn3:
    input:
        f"k{{k}}/{sample1}-{sample2}-{sample3}.k{{k}}.res"
    output:
        [f"{sample1}-{sample2}-{sample3}.k{{k}}.venn.{fmt}" 
                for fmt in ('png', 'pdf', 'svg')]
    run:
        import matplotlib as mpl 
        mpl.use('Agg')
        import matplotlib.pyplot as plt 
        from matplotlib_venn import venn3, venn3_circles
        
        res_db = {}
        with open(input[0], 'r') as fp:
            for line in fp:
                conbine, number = line.strip().split()
                number = int(number)
                res_db[conbine] = number
        
        venn3(res_db, set_labels=SAMPLES, subset_label_formatter=lambda x: f"{x:,}",
             alpha=0.2)
        venn3_circles(res_db, linewidth=1, alpha=0.6)
        for out in output:
            plt.savefig(out, dpi=300)
