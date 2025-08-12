#!/bin/bash
#SBATCH --job-name=snakemake
#SBATCH --partition=RM-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3-00:00:00

###ADJUST FOLLOWING TWO LINES TO ACTIVATE YOUR SNAKEMAKE ENVIRONMENT###
###source /jet/home/agomez3/.bashrc
###source activate /jet/home/agomez3/miniforge3/envs/snakemake

snakemake --unlock
snakemake --use-conda --rerun-incomplete --workflow-profile profile/
