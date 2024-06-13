#!/usr/bin/env -S bash -l
# Add -l to have module subsystem working on the new Ubuntu nodes: https://scicomp.ethz.ch/wiki/Migrating_to_Ubuntu
# Add -S to add an argument to the command line: https://stackoverflow.com/questions/4303128/how-to-use-multiple-arguments-for-awk-with-a-shebang-i-e
# Interactive run: srun -J pty-$(hostname) --ntasks=10 --mem-per-cpu=10G --time=1-0 --gpus=rtx_4090:1 --gres=gpumem:16g --tmp=16384 --pty bash -l
# Cleanup/re-run: rm -rf openfold openfold-setup.stderr.txt openfold-setup.stdout.txt; mamba env remove -n openfold_env -y
# GPU models: https://scicomp.ethz.ch/wiki/Using_the_batch_system#GPU
#SBATCH --job-name=openfold-setup
#SBATCH --ntasks=10
#SBATCH --mem-per-cpu=10G
#SBATCH --time=0-01:00:00
#SBATCH --gpus=rtx_4090:1
#SBATCH --gres=gpumem:16g
#SBATCH --tmp=16834
#SBATCH --output=openfold-setup.stdout.txt
#SBATCH --error=openfold-setup.stderr.txt
#source /cluster/home/jjaenes/.bashrc
module load eth_proxy
module load stack/2024-05

# /cluster/project/beltrao/jjaenes/software/miniconda3/envs/openfold_env/lib/python3.10/site-packages/torch/utils/cpp_extension.py:414: UserWarning: The detected CUDA version (12.2) has a minor version mismatch with the version that was used to compile PyTorch (12.1). Most likely this shouldn't be a problem.
#module load cuda/12.1.1
module load cuda/12.2.1

# Checkout pl_upgrades branch for CUDA12 support: https://openfold.readthedocs.io/en/latest/Installation.html
git clone https://github.com/aqlaboratory/openfold.git
cd openfold
git checkout pl_upgrades

# Adapted from: https://github.com/aqlaboratory/openfold/blob/pl_upgrades/environment.yml
# Adds latest compatible/matching versions for gcc, cudatoolkit-dev
echo "Download CUTLASS, required for Deepspeed Evoformer attention kernel"
git clone https://github.com/NVIDIA/cutlass --depth 1

mamba env create -n openfold_env -f ../../workflow/envs/openfold-eu.yaml

conda activate openfold_env

# Use max number of cores allocated by SLURM
export MAX_JOBS=$SLURM_NTASKS

##scripts/install_third_party_dependencies.sh
# Download folding resources
wget -N --no-check-certificate -P openfold/resources \
    https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

# Certain tests need access to this file
mkdir -p tests/test_data/alphafold/common
ln -rs openfold/resources/stereo_chemical_props.txt tests/test_data/alphafold/common

# Decompress test data
gunzip -c tests/test_data/sample_feats.pickle.gz > tests/test_data/sample_feats.pickle

python setup.py install

#echo "Download CUTLASS, required for Deepspeed Evoformer attention kernel"
#git clone https://github.com/NVIDIA/cutlass --depth 1
conda env config vars set CUTLASS_PATH=$PWD/cutlass

# This setting is used to fix a worker assignment issue during data loading
conda env config vars set KMP_AFFINITY=none

export LIBRARY_PATH=$CONDA_PREFIX/lib:$LIBRARY_PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

##
conda deactivate
conda activate openfold_env
# Use max number of cores allocated by SLURM
export MAX_JOBS=$SLURM_NTASKS

export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

scripts/download_alphafold_params.sh openfold/resources
scripts/run_unit_tests.sh
# run with conda run -n instead to test?
