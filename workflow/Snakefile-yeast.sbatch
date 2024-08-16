#!/usr/bin/env bash
#SBATCH --job-name=yeast
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=7-00:00:00
#SBATCH --output=results/yeast.output.txt
source $HOME/.bashrc
module load eth_proxy
conda activate af2genomics-envA
cd -P $HOME/af2genomics/
date; time snakemake --profile smk-simple-slurm-eu --snakefile workflow/Snakefile-yeast.smk yeast_all --dry-run
