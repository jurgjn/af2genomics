
import io, subprocess, pandas as pd

from af2genomics.common import *

foldx_binary = '/cluster/project/beltrao/jjaenes/software/foldx/foldx_20241231'

def parse_ifresid(s):
    aa_pos = int(s[2:])
    aa_ref = s[0]
    chain = s[1]
    return chain, aa_ref, aa_pos

def AnalyseComplex(pdb):
    pdb_split = os.path.split(pdb)
    interaction_id = os.path.splitext(pdb_split[1])[0]
    output_dir = os.environ['TMPDIR']
    ps_ = subprocess.run([foldx_binary, 
        '--command=AnalyseComplex', '--analyseComplexChains=A,B',
        f'--pdb-dir={pdb_split[0]}',
        f'--pdb={pdb_split[1]}',
        f'--output-dir={output_dir}'
    ], capture_output=True)
    #print('foldx output:')
    #print(ps_.stdout.decode('ascii'))

    Indiv_energies_ = os.path.join(output_dir, f'Indiv_energies_{interaction_id}_AC.fxout')
    Interaction_ = os.path.join(output_dir, f'Interaction_{interaction_id}_AC.fxout')
    Interface_Residues_ = os.path.join(output_dir, f'Interface_Residues_{interaction_id}_AC.fxout')
    Summary_ = os.path.join(output_dir, f'Summary_{interaction_id}_AC.fxout')

    assert os.path.isfile(Indiv_energies_) and os.path.isfile(Interaction_) and os.path.isfile(Interface_Residues_) and os.path.isfile(Summary_)
    #!cat {Indiv_energies_}
    #!cat {Interaction_}
    #!cat {Interface_Residues_}
    #!cat {Summary_}

    df_Indiv_energies_ = pd.read_csv(Indiv_energies_, sep='\t', skiprows=8)
    df_Interaction_ = pd.read_csv(Interaction_, sep='\t', skiprows=8)
    df_Interface_Residues_ = pd.read_csv(Interface_Residues_, skiprows=9)
    df_Summary_ = pd.read_csv(Summary_, sep='\t', skiprows=8)

    # Raw DataFrames:
    #return df_Indiv_energies_, df_Interaction_, df_Interface_Residues_, df_Summary_

    # Interface residues:
    resid_ = df_Interface_Residues_.squeeze()
    #l_resid = list(filter(len, df_Interface_Residues_.squeeze().split('\t')))
    #return l_resid

    # Interaction energy
    interaction_energy = df_Interaction_['Interaction Energy'].squeeze()
    return interaction_energy, resid_
