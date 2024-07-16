#!/usr/bin/env bash
#SBATCH --job-name=run_multimer
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=10G
#SBATCH --gpus=rtx_4090:1
#SBATCH --gres=gpumem:16g
#SBATCH --tmp=16834
#SBATCH --time=3-00:00:00
#SBATCH --output=results/openfold/run_multimer.stdout.txt
#SBATCH --error=results/openfold/run_multimer.stderr.txt
source /cluster/home/jjaenes/.bashrc
module load eth_proxy stack/2024-05 gcc/13.2.0 cuda/12.2.1
conda activate openfold_env
cd /cluster/work/beltrao/jjaenes/24.06.10_af2genomics
profile_euler/run_local run_pretrained_openfold_multimer --snakefile workflow/targets/openfold.smk
