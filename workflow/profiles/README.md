The default profile maps snakemake [standard resources](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#standard-resources)
to [Euler equivalents](https://scicomp.ethz.ch/wiki/Using_the_batch_system#Resource_requirements)
by adapting [jdblischak/smk-simple-slurm](https://github.com/jdblischak/smk-simple-slurm).

- `threads` sets number of cores; mapped to `--ntasks`
- `mem_mb` sets total RAM; mapped to `--mem-per-cpu` by accounting for cores: `expr {resources.mem_mb} / {threads}`
- `disk_mb` sets Euler's [local scratch space](https://scicomp.ethz.ch/wiki/Using_local_scratch); mapped to `--tmp`
- `runtime` sets max wall-clock time; mapped to `--time` and converted for `sbatch` if needed, e.g. `1d` becomes `'1-00:00:00'`
- `xtra_args` are directly added to `sbatch` to specify [conditional resources](https://github.com/jdblischak/smk-simple-slurm/tree/main/examples/conditional-resource) such as [GPUs](https://scicomp.ethz.ch/wiki/Using_the_batch_system#GPU)
