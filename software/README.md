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

```
.local/bin/code-server --install-extension ms-python.python
.local/bin/code-server --install-extension ms-toolsai.jupyter
.local/bin/code-server --install-extension snakemake.snakemake-lang
.local/bin/code-server --install-extension anwar.papyrus-pdf
.local/bin/code-server --install-extension mechatroner.rainbow-csv
#.local/bin/code-server --install-extension reageyao.biosyntax # auto-switches theme??
```

