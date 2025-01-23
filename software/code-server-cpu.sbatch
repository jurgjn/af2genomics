#!/usr/bin/env bash
#SBATCH --job-name=code-server-cpu
#SBATCH --ntasks=16
#SBATCH --mem-per-cpu=8G
#SBATCH --tmp=16834
#SBATCH --time=0-12:00:00
source /cluster/home/jjaenes/.bashrc
module load eth_proxy
VSC_IP_REMOTE="$(hostname -i)"
echo $VSC_IP_REMOTE
cd
.local/bin/code-server --bind-addr=${VSC_IP_REMOTE}:8897
