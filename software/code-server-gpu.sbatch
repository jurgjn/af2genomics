#!/usr/bin/env bash
#SBATCH --job-name=code-server-gpu
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=4G
#SBATCH --gpus=rtx_4090:1
#SBATCH --gres=gpumem:16g
#SBATCH --tmp=16834
#SBATCH --time=0-12:00:00
#source /cluster/home/jjaenes/.bashrc
module load eth_proxy
module load cuda/12.2.1
conda activate af2genomics-env
VSC_IP_REMOTE="$(hostname -i)"
echo $VSC_IP_REMOTE
cd
.local/bin/code-server --bind-addr=${VSC_IP_REMOTE}:8898
