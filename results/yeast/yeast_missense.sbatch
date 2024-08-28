#!/usr/bin/env bash
#SBATCH --job-name=yeast_missense
#SBATCH --ntasks=8
#SBATCH --mem-per-cpu=32G
#SBATCH --tmp=16834
#SBATCH --time=0-12:00:00
#SBATCH --output=results/yeast/yeast_missense_output.txt
source $HOME/.bashrc
conda activate af2genomics-envA
cd -P ~/af2genomics
jupyter nbconvert --to notebook --inplace --execute results/yeast/yeast_missense.ipynb
myjobs -j $SLURM_JOB_ID
