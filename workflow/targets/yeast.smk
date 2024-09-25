
def yeast_af2():
    df_ = pd.read_csv('results/yeast/af2.tsv', names=['pdb_gz'])
    df_['struct_id'] = df_['pdb_gz'].map(lambda pdb_gz: os.path.basename(pdb_gz).removesuffix('.pdb.gz'))
    #return df_.head(100).struct_id.tolist()
    return df_.struct_id.tolist()

def yeast_af2_trim_bf():
    df_ = pd.read_csv('results/yeast/af2.trim_bf.tsv', sep='\t')
    printlen(df_, 'raw structures')
    df_ = df_.query('nresid_trim_bf >= 16')
    printlen(df_, 'after filtering')
    return df_.struct_id.tolist()#[:200]

rule all:
    # eu-login-39 $ srun -J pty-$(hostname) --ntasks=1 --mem-per-cpu=16G --time=0-12 --tmp=16384 --pty bash
    # $ cd -P af2genomics
    # $ conda activate af2genomics-env
    # $ profile_euler/run_local all --dry-run
    # $ profile_euler/run_slurm all --dry-run
    # $ profile_euler/run_slurm all --keep-going --dry-run
    # $ profile_euler/run_slurm all --delete-temp-output --dry-run
    # $ sbatch workflow/targets/yeast.slurm
    input:
        #[ pfile(struct_id=struct_id, step='af2.trim_bf', suffix='.pdb', base='results/yeast') for struct_id in yeast_af2() ],
        [ pfile(struct_id=struct_id, step='af2.trim_bf.repairpdb.pssm', suffix='.tsv', base='results/yeast') for struct_id in yeast_af2_trim_bf() ],
