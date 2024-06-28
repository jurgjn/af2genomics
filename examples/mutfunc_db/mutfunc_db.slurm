#!/usr/bin/env bash
#SBATCH --job-name=mutfunc_db
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=32G
#SBATCH --time=0-23:00:00
#SBATCH --output=mutfunc_db.stdout.txt
#SBATCH --error=mutfunc_db.stderr.txt
source /cluster/home/jjaenes/.bashrc
conda activate af2genomics-env
cd /cluster/work/beltrao/jjaenes/24.06.10_af2genomics/examples/mutfunc_db
time jupyter nbconvert --to notebook --execute mutfunc_db.ipynb