
from af2genomics.common import *

import pandas as pd, openbabel, openbabel.pybel

def smi_largest_component(smiles):
    #https://github.com/dkoes/rdkit-scripts/blob/master/rdconf.py#L80
    if '.' in smiles:
        return max(smiles.split('.'), key=len)
    else:
        return smiles

def smi_strip_chemaxon(smiles):
    # Remove Chemaxon extensions from SMILES string
    # https://chemistry.stackexchange.com/questions/47242/what-is-this-smiles-notation-csicc-11
    if ' |' in smiles:
        smiles = smiles.split()[0]
    return smiles

def read_gnina(fp):
    def get_(mol_, key_):
        return float(mol_.data[key_]) if key_ in mol_.data else float('NaN')
    cols_ = ['minimizedAffinity', 'CNNscore', 'CNNaffinity', 'CNN_VS']
    df_ = pd.DataFrame.from_records([[
        get_(mol, 'minimizedAffinity'),
        get_(mol, 'CNNscore'),
        get_(mol, 'CNNaffinity'),
        get_(mol, 'CNN_VS'),
    ] for mol in openbabel.pybel.readfile('sdf', fp) ], columns=cols_)
    return df_

def read_sdf(fp, drop_=['MOL Chiral Flag', 'OpenBabel Symmetry Classes']):
    cols_ = ['ligand_id'] + next(openbabel.pybel.readfile('sdf', fp)).data.keys()
    df_ = pd.DataFrame.from_records(([mol.title] + mol.data.values() for mol in openbabel.pybel.readfile('sdf', fp)), columns=cols_)
    df_ = df_.drop(drop_, axis=1)
    df_.insert(loc=1, column='mode_id', value=df_.groupby('ligand_id', as_index=False).cumcount() + 1)
    return df_

def read_Corsello2020_compounds():
    fp_ = workpath('23.07.25_PRISM_Repurposing_23Q2/43018_2019_18_MOESM2_ESM.xlsx')
    df_ = pd.read_excel(fp_, sheet_name=1, header=2).rename({'name': 'compound_id'}, axis=1)[['compound_id', 'smiles']].query('smiles == smiles')
    printlen(df_, 'raw records')
    df_['compound_id'] = df_['compound_id'].str.lower()
    df_['smiles'] = df_['smiles'].str.split(',')
    df_ = df_.explode('smiles')
    printlen(df_, 'after parsing smiles')
    df_['smiles'] = df_['smiles'].str.lstrip().str.rstrip()
    df_ = df_.query('smiles != ""').reset_index(drop=True)

    df_['compound_id_penc'] = df_['compound_id'].map(penc)
    df_ = df_.drop_duplicates(subset=['compound_id_penc'], keep='first').reset_index(drop=True)
    printlen(df_, 'after de-duplicating by smiles')

    df_ = pd.concat([df_, pd.DataFrame.from_dict(list(df_['smiles'].map(lambda smiles: openbabel.pybel.readstring('smi', smiles).calcdesc())))], axis=1)
    df_ = df_.query('(150 < MW) & (MW < 1000)').reset_index(drop=True)
    printlen(df_, 'after filtering for molecular weight')
    df_['smiles'] = [* map(lambda smiles: smi_largest_component(smi_strip_chemaxon(smiles)), df_['smiles']) ]

    printlen(df_, 'rows with', uf(df_['compound_id_penc'].nunique()), 'unique values in compound_id_penc')
    return df_

def read_Piazza2020_Fig3():
    fp_ = '/cluster/project/beltrao/jjaenes/22.12_pocketomes/adhoc/23.02_Piazza2020_remap/Piazza2020_Fig3.tab'
    df_ = pd.read_csv(fp_, sep='\t', skiprows=2)
    df_[['uniprot_id', 'cofactor_id']] = df_['alphafill_id'].str.split('_', expand=True)
    return df_
