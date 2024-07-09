#!/usr/bin/env bash
# Setup databases as in: https://colabfold.mmseqs.com
#SBATCH --job-name=setup_databases
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=200G # Batch searches with many sequences against the ColabFoldDB quires a machine with approx. 128GB RAM - https://colabfold.mmseqs.com
#SBATCH --time=2-00:00:00
#SBATCH --output=setup_databases_eu.stdout.txt
#SBATCH --error=setup_databases_eu.stderr.txt
source /cluster/home/jjaenes/.bashrc
conda activate colabfold-env
cd /cluster/work/beltrao/jjaenes/24.06.10_af2genomics/software/colabfold
mkdir -p database
date; time ./setup_databases.sh database/
