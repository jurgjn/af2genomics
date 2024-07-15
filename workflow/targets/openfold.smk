
rule run_unit_tests:
    """
    openfold rules need snakemake started in openfold_env
        $ module load stack/2024-05 gcc/13.2.0 cuda/12.2.1
        $ conda activate openfold_env
        $ profile_euler/run_local run_unit_tests --snakefile workflow/targets/openfold.smk --dry-run
    """
    params:
        openfold_dir = '/cluster/work/beltrao/jjaenes/24.06.10_af2genomics/software/openfold',
    shell: """
        module list
        cd {params.openfold_dir}
        date; time scripts/run_unit_tests.sh
    """

rule fasta_dir:
    # eu-login-39 $ srun -J pty-$(hostname) --ntasks=1 --mem-per-cpu=16G --time=0-12 --tmp=16384 --pty bash
    # $ cd -P af2genomics
    # $ conda activate af2genomics-env
    # $ profile_euler/run_local fasta_dir --dry-run
    # $ profile_euler/run_slurm fasta_dir --dry-run
    # $ profile_euler/run_slurm fasta_dir --keep-going --dry-run
    # $ profile_euler/run_slurm fasta_dir --delete-temp-output --dry-run
    output:
        fasta_dir = directory('results/openfold/fasta_dir'),
    shell: """
        mkdir -p {output.fasta_dir}
        workflow/scripts/curl_uniprot Q71U36 > {output.fasta_dir}/Q71U36.fasta
        workflow/scripts/curl_uniprot P07437 > {output.fasta_dir}/P07437.fasta
        workflow/scripts/curl_uniprot O43663 > {output.fasta_dir}/O43663.fasta
    """

rule precompute_alignments:
    # $ rm -rf results/openfold/fasta_msas
    # $ profile_euler/run_local precompute_alignments --snakefile workflow/targets/openfold.smk --dry-run
    # $ profile_euler/run_sbatch precompute_alignments --snakefile workflow/targets/openfold.smk --dry-run
    input:
        dir = 'results/openfold/fasta_dir',
    output:
        dir = directory('results/openfold/fasta_msas'),
    params:
        openfold_dir = '/cluster/work/beltrao/jjaenes/24.06.10_af2genomics/software/openfold',
        workdir = '/cluster/work/beltrao/jjaenes/24.06.10_af2genomics',
    shell: """
        mkdir -p {output.dir}
        module list
        cd {params.openfold_dir}
        # Turn off templates --max_template_date '1950-01-01' from: https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1995177985/Using+AlphaFold+on+O2
        export BASE_DATA_DIR='/cluster/project/alphafold'
        export TEMPLATE_MMCIF_DIR=$BASE_DATA_DIR/pdb_mmcif/mmcif_files/
        date; time python3 scripts/precompute_alignments.py \
            --uniref90_database_path $BASE_DATA_DIR/uniref90/uniref90.fasta \
            --mgnify_database_path $BASE_DATA_DIR/mgnify/mgy_clusters_2018_12.fa \
            --pdb_seqres_database_path $BASE_DATA_DIR/pdb_seqres/pdb_seqres.txt \
            --uniref30_database_path $BASE_DATA_DIR/uniref30/UniRef30_2021_03 \
            --uniprot_database_path $BASE_DATA_DIR/uniprot/uniprot.fasta \
            --bfd_database_path $BASE_DATA_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
            --jackhmmer_binary_path jackhmmer \
            --hhblits_binary_path hhblits \
            --hmmsearch_binary_path hmmsearch \
            --hmmbuild_binary_path hmmbuild \
            --kalign_binary_path kalign \
            --max_template_date '1950-01-01' \
            --cpus_per_task 8 \
            {params.workdir}/{input.dir} {params.workdir}/{output.dir}
        cd -
    """

'''
rule run_pretrained_openfold:
# Turn off templates --max_template_date '1950-01-01' from: https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1995177985/Using+AlphaFold+on+O2
export BASE_DATA_DIR='/cluster/project/alphafold'
export TEMPLATE_MMCIF_DIR=$BASE_DATA_DIR/pdb_mmcif/mmcif_files
export INPUT_FASTA_DIR='/cluster/work/beltrao/jjaenes/24.06.10_af2genomics/examples/openfold/multimer_input_dir'
export OUTPUT_DIR='/cluster/work/beltrao/jjaenes/24.06.10_af2genomics/examples/openfold/multimer_output_dir'
date; time python3 run_pretrained_openfold.py $INPUT_FASTA_DIR $TEMPLATE_MMCIF_DIR --output_dir $OUTPUT_DIR \
    --uniref90_database_path $BASE_DATA_DIR/uniref90/uniref90.fasta \
    --mgnify_database_path $BASE_DATA_DIR/mgnify/mgy_clusters_2018_12.fa \
    --pdb_seqres_database_path $BASE_DATA_DIR/pdb_seqres/pdb_seqres.txt \
    --uniref30_database_path $BASE_DATA_DIR/uniref30/UniRef30_2021_03 \
    --uniprot_database_path $BASE_DATA_DIR/uniprot/uniprot.fasta \
    --bfd_database_path $BASE_DATA_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt \
    --jackhmmer_binary_path jackhmmer \
    --hhblits_binary_path hhblits \
    --hmmsearch_binary_path hmmsearch \
    --hmmbuild_binary_path hmmbuild \
    --kalign_binary_path kalign \
    --config_preset model_1_multimer_v3 \
    --max_template_date '1950-01-01' \
    --model_device 'cuda:0'
'''

# Caveats:
# - Set variables during environment initialisation
# - Create symlink to local AF2 weights
# - Move?: Emitting ninja build file /cluster/home/jjaenes/.cache/torch_extensions/py310_cu121/evoformer_attn/build.ninja...

## Run monomer
#export BASE_DATA_DIR='/cluster/project/alphafold'
#export TEMPLATE_MMCIF_DIR=$BASE_DATA_DIR/pdb_mmcif/mmcif_files

# /cluster/home/jjaenes/work/24.06.11_openfold_build
#export INPUT_DIR='/cluster/home/jjaenes/work/24.06.11_openfold_build/example_inp'
#export OUTPUT_DIR='/cluster/home/jjaenes/work/24.06.11_openfold_build/example_out'
# wget https://rest.uniprot.org/uniprotkb/O00244.fasta

# Needs to run in the OpenFold directory to find the AF2 weights
#python3 run_pretrained_openfold.py $INPUT_DIR $TEMPLATE_MMCIF_DIR --output_dir $OUTPUT_DIR --config_preset model_1_ptm --uniref90_database_path $BASE_DATA_DIR/uniref90/uniref90.fasta --mgnify_database_path $BASE_DATA_DIR/mgnify/mgy_clusters_2018_12.fa --pdb70_database_path $BASE_DATA_DIR/pdb70/pdb70 --uniclust30_database_path $BASE_DATA_DIR/uniclust30/uniclust30_2018_08/uniclust30_2018_08 --bfd_database_path $BASE_DATA_DIR/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt --model_device "cuda:0"


#source /cluster/home/jjaenes/.bashrc
#conda activate af2genomics-env
#cd /cluster/home/jjaenes/work/23.12.06_ppi_reselect
#time jupyter nbconvert --to notebook --execute af2-models-split-ifresid.ipynb
