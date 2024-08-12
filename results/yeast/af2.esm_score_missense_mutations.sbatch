#!/usr/bin/env bash
#SBATCH --job-name=esm_score_missense_mutations
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=10G
#SBATCH --gpus=rtx_4090:1
#SBATCH --gres=gpumem:16g
#SBATCH --tmp=16834
#SBATCH --time=7-00:00:00
#SBATCH --output=results/yeast/af2.esm_score_missense_mutations_output.txt
source $HOME/.bashrc
module load eth_proxy
module load stack/2024-05 gcc/13.2.0 cuda/12.2.1
cd /cluster/work/beltrao/jjaenes/24.06.10_af2genomics/software/esm-variants
date; time conda run -n esm-variants-env python3 esm_score_missense_mutations.py \
    --input-fasta-file ../../results/yeast/af2.esm_score_missense_mutations.fasta \
    --output-csv-file ../../results/yeast/af2.esm_score_missense_mutations.csv
