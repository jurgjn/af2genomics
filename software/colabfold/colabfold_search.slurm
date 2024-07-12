#!/usr/bin/env bash
# Setup databases as in: https://colabfold.mmseqs.com
#SBATCH --job-name=colabfold_search
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=20G # Batch searches with many sequences against the ColabFoldDB quires a machine with approx. 128GB RAM - https://colabfold.mmseqs.com
#SBATCH --time=1-00:00:00
#SBATCH --output=colabfold_search.stdout.txt
#SBATCH --error=colabfold_search.stderr.txt
source /cluster/home/jjaenes/.bashrc
cd /cluster/work/beltrao/jjaenes/24.06.10_af2genomics/software/colabfold
module load eth_proxy
conda activate colabfold-env
time jupyter nbconvert --to notebook --inplace --execute colabfold_search.ipynb
