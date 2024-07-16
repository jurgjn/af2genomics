
rule fasta_file:
    output:
        fasta = 'results/openfold/fasta_dir/{uniprot_id}.fasta',
    shell: """
        workflow/scripts/curl_uniprot {wildcards.uniprot_id} > {output.fasta}
    """

def fasta_dir_input():
    import pandas as pd, af2genomics
    return ['Q71U36', 'P07437'] + pd.read_excel(af2genomics.workpath('24.07.08_tubb/Protein list for test run.xlsx'), names=['uniprot_id'])['uniprot_id'].tolist()

rule fasta_dir:
    # $ conda activate af2genomics-env
    # $ profile_euler/run_local fasta_dir --snakefile workflow/targets/openfold.smk --dry-run
    input:
        expand('results/openfold/fasta_dir/{uniprot_id}.fasta', uniprot_id=fasta_dir_input()),

rule precompute_alignments:
    # $ rm -rf results/openfold/fasta_msas
    # $ conda activate openfold_env
    # $ profile_euler/run_local precompute_alignments --snakefile workflow/targets/openfold.smk --dry-run
    # $ sbatch results/openfold/fasta_msas.slurm
    input:
        dir = 'results/openfold/fasta_dir',
    output:
        dir = directory('results/openfold/fasta_msas'),
        sstat = 'results/openfold/fasta_msas/.sstat.tsv',
    threads: 8
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
            --cpus_per_task {threads} \
            {params.workdir}/{input.dir} {params.workdir}/{output.dir}
        cd -
        workflow/scripts/sstat_eu > {output.sstat}
    """

rule run_pretrained_openfold_multimer:
    # mkdir -p fasta_dir_multimer
    # cat fasta_dir/O43663.fasta fasta_dir/P07437.fasta fasta_dir/Q71U36.fasta > fasta_dir_multimer/O43663_P07437_Q71U36.fasta
    # profile_euler/run_local run_pretrained_openfold_multimer --snakefile workflow/targets/openfold.smk --dry-run
    # https://openfold.readthedocs.io/en/latest/Inference.html#advanced-options-for-increasing-efficiency
    input:
        fasta_dir = 'results/openfold/fasta_dir_multimer',
        use_precomputed_alignments = 'results/openfold/fasta_msas',
    output:
        dir = directory('results/openfold/run_multimer'),
        sstat = 'results/openfold/run_multimer/.sstat.tsv',
    threads: 8
    params:
        openfold_dir = '/cluster/work/beltrao/jjaenes/24.06.10_af2genomics/software/openfold',
        workdir = '/cluster/work/beltrao/jjaenes/24.06.10_af2genomics',
    shell: """
        mkdir -p {output.dir}
        module list
        cd {params.openfold_dir}
        # Turn off templates --max_template_date '1950-01-01' from: https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1995177985/Using+AlphaFold+on+O2
        export BASE_DATA_DIR='/cluster/project/alphafold'
        export TEMPLATE_MMCIF_DIR=$BASE_DATA_DIR/pdb_mmcif/mmcif_files
        date; time python3 run_pretrained_openfold.py {params.workdir}/{input.fasta_dir} $TEMPLATE_MMCIF_DIR \
            --output_dir {params.workdir}/{output.dir} \
            --use_precomputed_alignments {params.workdir}/{input.use_precomputed_alignments} \
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
            --cpus {threads} \
            --skip_relaxation \
            --use_deepspeed_evoformer_attention \
            --trace_model \
            --model_device 'cuda:0'
        cd -
        workflow/scripts/sstat_eu > {output.sstat}
    """

rule run_pretrained_openfold_help:
    # profile_euler/run_local run_pretrained_openfold_help --snakefile workflow/targets/openfold.smk --dry-run
    params:
        openfold_dir = '/cluster/work/beltrao/jjaenes/24.06.10_af2genomics/software/openfold',
        workdir = '/cluster/work/beltrao/jjaenes/24.06.10_af2genomics',
    shell: """
        cd {params.openfold_dir}
        date; time python3 run_pretrained_openfold.py --help
    """

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
