#!/usr/bin/env bash
#SBATCH --job-name=openfold.slurm
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20G
#SBATCH --time=3-00:00:00
#SBATCH --output=results/openfold.stdout.txt
#SBATCH --error=results/openfold.stderr.txt
source /cluster/home/jjaenes/.bashrc
module load eth_proxy
conda activate openfold_env
cd /cluster/work/beltrao/jjaenes/24.06.10_af2genomics/
smk sbatch precompute_alignments_ --snakefile workflow/targets/openfold.smk
