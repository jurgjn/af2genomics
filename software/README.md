```
conda env create --file AlphaCutter/environment.yml
```

```
sbatch openfold-setup.slurm
```

```
wget https://github.com/gnina/gnina/releases/download/v1.1/gnina
chmod u+x gnina
```

```
pip install bash_kernel
python -m bash_kernel.install
```

## conda
* `nvcc` is in `cudatoolkit` and/or `cudatoolkit-dev`, e.g. [Problem with openfold installation.  · Issue \#217 · gcorso/DiffDock](https://github.com/gcorso/DiffDock/issues/217)
* `gcc` , `g++`, ... are in `gcc_linux-64` and/or `gxx_linux-64`, e.g. [Installing Difficult Conda Packages on O2 - HMS IT RC O2 - Confluence](https://harvardmed.atlassian.net/wiki/spaces/O2/pages/2557575188/Installing+Difficult+Conda+Packages+on+O2)
