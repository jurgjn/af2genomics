#!/usr/bin/env bash
#SBATCH --job-name=multimer
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=32G
#SBATCH --time=0-23:00:00
#SBATCH --gpus=rtx_4090:1
#SBATCH --gres=gpumem:16g
#SBATCH --tmp=16834
#SBATCH --output=multimer.stdout.txt
#SBATCH --error=multimer.stderr.txt
source /cluster/home/jjaenes/.bashrc
conda activate af2genomics-env
cd /cluster/work/beltrao/jjaenes/24.06.10_af2genomics/examples/openfold
time jupyter nbconvert --to notebook --execute openfold-run-multimer.ipynb