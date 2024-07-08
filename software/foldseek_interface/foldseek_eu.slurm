#!/usr/bin/env bash
# Compile foldseek interface on oldest-generation Euler currently use for compatibility (by specifying --constraint=ibfabric6)..
#SBATCH --job-name=foldseek-eu
#SBATCH --ntasks=1
#SBATCH --constraint=ibfabric6
#SBATCH --mem-per-cpu=32G
#SBATCH --time=0-1:00:00
#SBATCH --output=foldseek_eu.stdout.txt
#SBATCH --error=foldseek_eu.stderr.txt
source /cluster/home/jjaenes/.bashrc
conda activate af2genomics-env
cd /cluster/work/beltrao/jjaenes/24.06.10_af2genomics/software/foldseek_interface
time jupyter nbconvert --to notebook --inplace --execute foldseek_eu.ipynb
