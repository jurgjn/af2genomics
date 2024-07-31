#!/usr/bin/env bash
#SBATCH --job-name=yeast
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G
#SBATCH --time=3-00:00:00
#SBATCH --output=results/yeast.stdout.txt
#SBATCH --error=results/yeast.stderr.txt
source /cluster/home/jjaenes/.bashrc
module load eth_proxy
conda activate af2genomics-envA
cd /cluster/work/beltrao/jjaenes/24.06.10_af2genomics/
date; time snakemake all --local-cores $SLURM_NTASKS
