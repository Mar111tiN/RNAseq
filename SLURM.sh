#!/bin/bash

# can be overruled on CLI with -J <NAME>
#SBATCH --job-name=RNAseq

# Set the file to write the stdout and stderr to (if -e is not set; -o or --output).
#SBATCH --output=slogs/%x-%j.log

# Set the number of cores (-n or --ntasks).
#SBATCH --ntasks=2

# Force allocation of the two cores on ONE node.
#SBATCH --nodes=1

# Set the memory per CPU. Units can be given in T|G|M|K.
#SBATCH --mem-per-cpu=2500M

# Set the partition to be used (-p or --partition).
#SBATCH --partition=medium

# Set the expected running time of your job (-t or --time).
# Formats are MM:SS, HH:MM:SS, Days-HH, Days-HH:MM, Days-HH:MM:SS
#SBATCH --time=20:00:00

SNAKE_HOME=$(pwd);

export LOGDIR=${HOME}/scratch/slogs/${SLURM_JOB_NAME}-${SLURM_JOB_ID}
export TMPDIR=/fast/users/${USER}/scratch/tmp;
mkdir -p $LOGDIR;

# make conda available
eval "$($(which conda) shell.bash hook)"
# activate snakemake env
conda activate snake-env;
echo $CONDA_PREFIX "activated";

set -x;

# !!! leading white space is important
SLURM_CLUSTER="sbatch -p {cluster.partition} -t {cluster.t} --mem-per-cpu={cluster.mem} -J {cluster.name} --nodes={cluster.nodes} -n {cluster.threads}";
SLURM_CLUSTER="$SLURM_CLUSTER -o ${LOGDIR}/{rule}-%j.log"
snakemake --unlock --rerun-incomplete;
snakemake --dag | awk '$0 ~ "digraph" {p=1} p' | dot -Tsvg > dax/dag.svg;
snakemake --use-conda --rerun-incomplete --cluster-config config/cluster/RNAseq-cluster.json --cluster "$SLURM_CLUSTER" -prk -j 1000;
# -k ..keep going if job fails
# -p ..print out shell commands
