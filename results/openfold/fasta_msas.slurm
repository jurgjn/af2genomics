#!/usr/bin/env bash
#SBATCH --job-name=precompute_alignments
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=10G
#SBATCH --tmp=16834
#SBATCH --time=3-00:00:00
#SBATCH --output=results/openfold/fasta_msas.stdout.txt
#SBATCH --error=results/openfold/fasta_msas.stderr.txt
source /cluster/home/jjaenes/.bashrc
conda activate openfold_env
cd /cluster/work/beltrao/jjaenes/24.06.10_af2genomics
profile_euler/run_local precompute_alignments --snakefile workflow/targets/openfold.smk
