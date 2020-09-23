##
## split homo protein by groups using tblastn


import os
import os.path as op
import sys

from Bio import SeqIO
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

homo_proteins = config['homo_proteins']
homo_part_number = config['homo_part_number']
prefix, suffix = op.splitext(homo_proteins)

numbers = [ "{:03}".format(i) for i in range(1, 1+homo_part_number)]
ncpus = config['ncpus']

evalue = config['evalue']
num_alignments = config['num_alignments']
outfmt = 6

rule all:
    input:
        expand("{group}/{group}.homo.fasta", group=groups_list)

rule splitHomoProtein:
    input:
        homo = homo_proteins,

    output:
        expand("{fasta}.split/{prefix}.part_{number}{suffix}", 
                fasta=(homo_proteins, ), number=numbers, 
                prefix=(prefix, ), suffix=(suffix, )) 
    log:
        "logs/splitHomoProtein.log"
    params:
        num = homo_part_number
    shell:
        "seqkit split -p {params.num} {input.homo} 2>{log}"

rule makeblastdb:
    input:
        "{group}/{group}.fasta"
    output:
        expand("{{group}}/{{group}}.fasta.db.{ext}", ext=['nsq', 'nin', 'nhr'])
    log:
        "logs/{group}.makeblastdb.log"
    params:
        dbtype = "nucl"
    shell:
        "makeblastdb -in {input} -out {input}.db -dbtype {params.dbtype} 2>{log}"

rule tblastn:
    input:
        target_fasta = "{group}/{group}.fasta",
        db =  expand("{{group}}/{{group}}.fasta.db.{ext}", ext=['nsq', 'nin', 'nhr']),
        fasta = f"{homo_proteins}.split/{prefix}.part_{{number}}{suffix}"
                
    output:
        f"{{group}}/results/{prefix}.part_{{number}}{suffix}.tsv"
    log:
        "logs/tblastn_{group}.log"
    threads: ncpus
    params:
        evalue = evalue,
        num_alignments = num_alignments,
        outfmt = outfmt
    shell:
        "tblastn -query {input.fasta} -db {input.target_fasta}.db -evalue {params.evalue} "
        "-num_threads {threads} -num_alignments {params.num_alignments} "
        "-outfmt {params.outfmt} -out {output}"

rule mergeBlastResult:
    input:
        expand("{{group}}/results/{prefix}.part_{number}{suffix}.tsv", 
                fasta=(homo_proteins, ), prefix=(prefix, ), 
                number=numbers, suffix=(suffix, )) 
    output:
        "{group}/{group}.tsv"
    shell:
        "cat {input} > {output}"

rule blast2fasta:
    input:
        "{group}/{group}.tsv"
    output:
        homoFasta = "{group}/{group}.homo.fasta"
    run:
        gene_set = set()
        with open(input[0]) as fp:
        
            for line in fp:
                if line.startswith("#"):
                    continue
                gene = line.strip().split()[0]
                gene_set.add(gene)

        with open(output.homoFasta, 'w') as out:
            fp = open(homo_proteins)
            fa = SeqIO.parse(fp, 'fasta')
            for record in fa:
                if record.id in gene_set:
                    SeqIO.write(record, out, 'fasta')