
def yeast_af2():
    df_ = pd.read_csv('results/yeast/af2.tsv', names=['pdb_gz'])
    df_['struct_id'] = df_['pdb_gz'].map(lambda pdb_gz: os.path.basename(pdb_gz).removesuffix('.pdb.gz'))
    return df_.head(2).tail(1).struct_id.tolist()

rule all:
    # eu-login-39 $ srun -J pty-$(hostname) --ntasks=1 --mem-per-cpu=16G --time=0-12 --tmp=16384 --pty bash
    # $ cd -P af2genomics
    # $ conda activate af2genomics-env
    # $ profile_euler/run_local all --dry-run
    # $ profile_euler/run_slurm all --dry-run
    # $ profile_euler/run_slurm all --keep-going --dry-run
    # $ profile_euler/run_slurm all --delete-temp-output --dry-run
    input:
        [ pfile(struct_id=struct_id, step='af2.trim_bf', suffix='.pdb', base='results/yeast') for struct_id in yeast_af2()[:1000] ],
 