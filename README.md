# Introduction

A nextflow-based workflow for amplicon sequencing data analysis.

# Install

```sh
# Install Qiime2
wget https://data.qiime2.org/distro/core/qiime2-2020.11-py36-linux-conda.yml
conda env create -n qiime2-2020.11 --file qiime2-2020.11-py36-linux-conda.yml
conda activate qiime2-2020.11
conda install -c bioconda bioconductor-phyloseq r-optparse -y
conda install -c rujinlong r-qiime2r
```

# Usage

```sh
nextflow run ampseq.nf -profile hpc_slurm -resume --mode "remove_adapters"

# Check p02_primer_trimmed.qzv and modify nextflow.config
qiime tools view p02_primer_trimmed.qzv

nextflow run ampseq.nf -profile hpc_slurm -resume --mode "all" --classifier false
```
