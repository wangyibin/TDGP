## kmer data venn plot
count kmer in a resequencing data and plot venn picture
### Dependency
- [`kmc`](https://github.com/refresh-bio/KMC)
- [`snakemake`](https://snakemake.readthedocs.io/en/stable/)
- [`matplotlib-venn`](https://pypi.org/project/matplotlib-venn/)
- [`matplotlib`](https://matplotlib.org/)

### Usage

- pre_data
    ```bash
    data/
        |-- HD_R1.fastq.gz
        |-- HD_R2.fastq.gz
        |-- JGY_R1.fastq.gz 
        |-- JGY_R2.fastq.gz
        |-- TGY_R1.fastq.gz
        `-- TGY_R2.fastq.gz
    ```
- run pipeline
    ```bash
    snakemake -s ~/code/TDGP/pipeline/kmerVenn/kmerVenn3.smk \
        --config samples="['TGY', 'HD', 'JGY']" -j 10 \
        --cluster 'qsub -l nodes=1:ppn={threads} -j oe -V -q workq' -p
    ```

