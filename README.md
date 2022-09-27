QUILT-wrap Wrapper around QUILT using Snakemake
=============================================

This folder contains scripts to facilitate imputing samples using snakemake, targetting performing large scale analysis on a cluster

A workflow could look like what follows, based on parameters specified in main.smk

Note that the below uses `run.sh`, a lightweight wrapper to facilitate running snakemake, for example
```
./run.sh all local 1 --dry-run
./run.sh all cluster 100 --dry-run
```
where the four arguments are the keyword of what to run (the target), whether to run locally or using a cluster (local vs target), how many cores to use, and additional optional arguments to snakemake (for example `--dry-run`).

## Install snakemake

For instance, using conda. See the main QUILT README for basic information about how to install a program using conda. 

```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
conda install --name snakemake snakemake
source activate snakemake
```

## Activate local environment

Here in the activate script we set variables specific to this analysis, like where we want the output to go, and the path to the analysis ready BAM files
```
source activate snakemake
. ./activate
```

## Download reference information

Install liftOver

```
cd ${ANALYSIS_DIR}
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod +x liftOver
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

On the cluster I use, only the head node has internet access, so we download files there

```
# Download recombination rate (IMPORTANT if relevant, change population, see main.smk parameters)
# Also download 1000 Genomes reference data
./run.sh download local 1

# Convert per-chormosome files 
./run.sh convert cluster 1000

```

## Prepare reference information

```
## Determine per-chromosome chunks
## this could be done using Snakemake programatically, but we do using R for simplicity
R -f determine_chunks.R

## Prepare per-region imputation files
./run.sh prep cluster 1000
```

## Impute samples
```
## Requires bamlist.txt to be defined in ANALYSIS_DIR
./run.sh impute cluster 1000 --dryrun
```



