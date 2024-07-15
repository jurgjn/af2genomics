
def some_human():
    #EGFR/P00533
    #ERBB2/P04626
    #TP53/P04637
    #MED1/Q15648
    #ATOX1/O00244
    #UBIAD1/Q9Y5Z9
    return ['P00533-F1', 'P04626-F1', 'P04637-F1', 'Q15648-F1', 'O00244-F1', 'Q9Y5Z9-F1', ]

rule all:
    # eu-login-39 $ srun -J pty-$(hostname) --ntasks=1 --mem-per-cpu=16G --time=0-12 --tmp=16384 --pty bash
    # $ cd -P af2genomics
    # $ conda activate af2genomics-env
    # $ profile_euler/run_local all --dry-run
    # $ profile_euler/run_slurm all --dry-run
    # $ profile_euler/run_slurm all --keep-going --dry-run
    # $ profile_euler/run_slurm all --delete-temp-output --dry-run
    input:
        [ pfile(struct_id=struct_id, step='af2.trim_bf.repairpdb.pssm', suffix='.tsv', base='results/human') for struct_id in some_human() ], # Test on a handful of human structures
        #[ pfile(struct_id=struct_id, step='af2.trim_bf.repairpdb.pssm', suffix='.zip', base='results/yeast') for struct_id in all_yeast() ],
 
