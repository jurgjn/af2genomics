
localrules: af2

rule af2:
    """
        gunzip -c resources/afdb/UP000005640_9606_HUMAN_v3/AF-{wildcards.struct_id}-F1-model_v3.pdb.gz > {output.pdb}
        tar -tf UP000005640_9606_HUMAN_v3.tar --wildcards "*.pdb.gz" > UP000005640_9606_HUMAN_v3.txt
        tar -xf UP000005640_9606_HUMAN_v3.tar --wildcards "*.pdb.gz"
        ---
        tar -xOf resources/alphafold/UP000005640_9606_HUMAN.tar AF-{wildcards.struct_id}-F1-model_v1.pdb.gz | gunzip -c > {output.pdb}
        gunzip -c resources/afdb/UP000005640_9606_HUMAN/AF-{wildcards.struct_id}-F1-model_v1.pdb.gz > {output.pdb}
        wget -O {output.pdb} https://alphafold.ebi.ac.uk/files/AF-{wildcards.struct_id}-F1-model_v1.pdb ||\
        wget -O {output.pdb} https://alphafold.ebi.ac.uk/files/AF-{wildcards.struct_id}-F1-model_v2.pdb ||\
        wget -O {output.pdb} https://alphafold.ebi.ac.uk/files/AF-{wildcards.struct_id}-F1-model_v3.pdb
    """
    output:
        pdb = pfile(struct_id='{}', step='af2', suffix='.pdb'),
    shell: """
        wget -O {output.pdb} https://alphafold.ebi.ac.uk/files/AF-{wildcards.struct_id}-model_v4.pdb
    """

localrules: trim_bf

rule trim_bf:
    """
    Moving-average pLDDT filter
    """
    input:
        pdb = pfile(struct_id='{}', step='{prev_steps}', suffix='.pdb'),
    output:
        pdb = pfile(struct_id='{}', step='{prev_steps}.trim_bf', suffix='.pdb'),
    shell: """
        pdb_trim_bf --pdbfile {input.pdb} --outpdb {output.pdb}
    """
