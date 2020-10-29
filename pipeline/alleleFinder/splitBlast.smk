


fasta = config['fasta']
query = config['query']
dbtype = config['dbtype']
nparts = config['nparts']
ncpus = config['ncpus']
out = config['out']

evalue = config['evalue']
outfmt = config['outfmt']
num_alignments = config['num_alignments']

query_string_list = query.split(".")
query_string_list.insert(-1, "part_{}")
query_format = ".".join(query_string_list)
split_query_list = [ query_format.format("{:0>3}".format(i)) 
                                for i in range(1, nparts + 1)]

if dbtype == 'nucl':
    suffix = ['nin', 'nhr', 'nsq']
    PROG = 'blastn'
elif dbtype == 'prot':
    suffix = ['pin', 'phr', 'psq']
    PROG = 'blastp'

rule all:
    input:
        out

rule splitQuery:
    input:
        query
    output:
        expand("{{query}}.split/{fa}", fa=split_query_list)
    params:
        nparts = nparts
    shell:
        "seqkit split -p {params.nparts} {input} 2> /dev/null"

rule makeblastdb:
    input:
        fasta
    output:
        expand("{fasta}.{suffix}", fasta=(fasta, ), suffix=suffix)
    log:
        "logs/makeblastdb.log"
    params:
        dbtype=dbtype
    shell:
        "makeblastdb -in {input} -out {input} -dbtype {params.dbtype} 2>{log}"
    
rule splitBlast:
    input:
        fa = f"{query}.split/{{sample}}",
        index = expand("{fasta}.{suffix}", fasta=(query, ), suffix=suffix)
    output:
        temp("{sample}.tmp.blast")
    log: 
        "logs/{sample}.blast.log"
    params:
        evalue = evalue,
        outfmt = outfmt,
        num_alignments = num_alignments
    threads:
        ncpus
    shell:
        "{PROG} -query {input.fa} -db {fasta} -out {output} "
        "-evalue {params.evalue} -num_threads {threads} "
        "-outfmt {params.outfmt} -num_alignments {params.num_alignments} 2>{log}"

rule mergeBlast:
    input:
        expand("{sample}.tmp.blast", sample=split_query_list)
    output:
        out
    shell:
        "cat {input} > {output}"