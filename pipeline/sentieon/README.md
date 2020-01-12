## Snakemake Pipeline of Sentieon

### Tutorial
1. method-1
- configure file
```bash
cp pipeline/sentieon/default.yaml .
```
and modify 
```yaml
## all sample with correct extense will be run in this directory
fastq_folder: 'fastq'

## fasta file 
fasta: 'fasta'
```
- run snakemake
```bash
snakemake -s pipeline/sentieon/Snakefile -c default.yaml 
```
2. method-2
- without modify yaml, add config parameters in command line
```bash
snakemake -s pipeline/sentieon/Snakefile --config fasta=fasta fastq_folder=fastq
```


### Advance
- execute with SGE cluster
```bash 
snakemake -s pipeline/sentieon/Snakefile -c default.yaml -j 5 --cluster "qsub -V -cwd -j y -q all.q -pe mpi {threads}"
```

