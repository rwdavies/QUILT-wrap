#!/usr/bin/env bash

set -e 

## a lightweight wrapper around snakemake

## e.g. ./run.sh all local 1 --dry-run
## e.g. ./run.sh all cluster 100 --dry-run

what=$1
what="${what:-all}"
where=$2
where="${where:-local}"
nCores=${3}
nCores="${nCores:-1}"
other="${4}"

if [ "${ANALYSIS_DIR}" == "" ]
then
    echo "Variable ANALYSIS_DIR not set. Did you forget to activate?"
    exit 1
fi

## put somewhere else, no space on home
LOG_DIR="${ANALYSIS_DIR}/logs/"

mkdir -p ${ANALYSIS_DIR}
mkdir -p ${LOG_DIR}
mkdir -p "${ANALYSIS_DIR}/results"

rsync -a ${BAMLIST} ${ANALYSIS_DIR}/bamlist.txt ## if not already there

cd ${ANALYSIS_DIR}
SCRIPTPATH="${QUILT_WRAP_HOME}"
SNAKEFILE="${SCRIPTPATH}main.smk"
SNAKEMAKE="snakemake"

## puts QUILT.R in PATH
export PATH=${QUILT_HOME}:${PATH}

if [ $where == "cluster" ]
then
    ## {params.queue}
    ## --max-status-checks-per-second 0.01 ## once per 100 seconds
    ${SNAKEMAKE} --max-status-checks-per-second 0.01 --snakefile ${SNAKEFILE} -w 30 --cluster "qsub -cwd -V -N {params.N} -pe shmem {params.threads} -q short.qc -P davies.prjc -j Y -o ${LOG_DIR}" --jobs ${nCores}  ${other} ${what}
elif [ $where == "local" ]
then
    ${SNAKEMAKE} --snakefile ${SNAKEFILE} ${other} --cores ${nCores} ${what}
    ## --dag > dag.input    
    ## dot -v -Tpdf dag.input  > dag.pdf
    echo done
else
    echo bad where: ${where}
    exit 1
fi

