#!/usr/bin/env bash
#SBATCH --job-name=code-server
#SBATCH --ntasks=2
#SBATCH --mem-per-cpu=4G
# Time format: d-hh:mm:ss
#SBATCH --time=0-12:00:00
#source /cluster/home/jjaenes/.bashrc
module load eth_proxy
#conda activate code-server-env
VSC_IP_REMOTE="$(hostname -i)"
echo $VSC_IP_REMOTE
#code-server --bind-addr=${VSC_IP_REMOTE}:8898
cd
.local/bin/code-server --bind-addr=${VSC_IP_REMOTE}:8898