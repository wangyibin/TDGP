##
## split homo protein by groups using tblastn


import os
import os.path as op
import sys

from Bio import SeqIO
from collections import OrderedDict
from pytools.persistent_dict import PersistentDict

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


homo_protein_list = [ i.strip() for i in \
    os.popen(f"seqkit seq -i -n {homo_proteins}") ]
homo_protein_list = list(set(homo_protein_list))
homo_protein_path_list = list(map(lambda x:f"{homo_proteins}.split/{prefix}.id_{x}{suffix}", homo_protein_list))
rule splitHomoProtein:
    input:
        homo = homo_proteins   

    output:
        f'{homo_proteins}.split.ok'
    log:
        "logs/splitHomoProtein.log"
    threads: ncpus
    shell:
        "seqkit split -j {threads} -i {input.homo} 2>{log} && touch {output}"

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

rule write_tblastn_cmd:
    input:
        f'{homo_proteins}.split.ok',
        target_fasta = "{group}/{group}.fasta",
        db =  expand("{{group}}/{{group,}}.fasta.db.{ext}", ext=['nsq', 'nin', 'nhr']),
        
    output:
        f"{{group}}/cmd.list"

    params:
        evalue = evalue,
        num_alignments = num_alignments,
        outfmt = outfmt,
        threads = ncpus
    run:
        if not op.exists(f'{wildcards.group}/results'):
            os.makedirs(f'{wildcards.group}/results')
        with open(output[0], 'w') as outfile:
            for fasta in homo_protein_path_list:
                base_fasta = op.basename(fasta)
                print(f"tblastn -query {fasta} -db {input.target_fasta}.db -evalue {params.evalue} "
                        f"-num_threads {params.threads} -num_alignments {params.num_alignments} "
                        f"-outfmt {params.outfmt} -out {wildcards.group}/results/{fasta}.tsv", file=outfile)

rule ParaBlast:
    input:
        "{group}/cmd.list"
    
    output:
        "{group}/ParaBlast.ok",
        # expand("{{group}}/results/{prefix}.id_{ID}{suffix}.tsv", 
        #         fasta=(homo_proteins, ), prefix=(prefix, ), 
        #         ID=homo_protein_list, suffix=(suffix, )) 
    log:
        "logs/tblastn_{group}.log"
    threads: ncpus
    shell:
        'ParaFly -c {input} -CPU {threads} -failed_cmds {wildcards.group}/FailedCommands 2>{log} &&'
        'if [[ ! -s {wildcards.group}/FailedCommands ]]; then touch {output[0]}; else {wildcards.group}/ParaBlast.failed; fi' 

# rule rescueFailedBlast:
#     input:
#         "{group}/ParaBlast.failed"
#     output:
#         "{group}/ParaBlast.ok"
#     log:
#         "logs/rescueFailedBlast_{group}.log"
#     threads:
#         ncpus
#     shell:
#         'if [[ ParaFly -c {wildcards.group}/FaileCommands -CPU {threads} -failed_cmds {wildcards.group}/FailedCommands.rescue'


rule mergeBlastResult:
    input:
        "{group}/ParaBlast.ok",
        # expand("{{group}}/results/{prefix}.id_{ID}{suffix}.tsv", 
        #         fasta=(homo_proteins, ), prefix=(prefix, ), 
        #         ID=homo_protein_list, suffix=(suffix, )) 
    output:
        "{group}/{group}.tsv"
    shell:
        "cat {wildcards.group}/results/*tsv > {output}"

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