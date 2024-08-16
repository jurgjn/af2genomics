
import pandas as pd

from af2genomics.common import *

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

def read_samples():
    fp_ = workpath('23.07.25_Corsello2017/repurposing_samples_20240610.txt')
    df_ = pd.read_csv(fp_, comment='!', sep='\t')#.rename({'pert_iname': 'compound_id'}, axis=1)
    printlen(df_, 'raw samples')
    printlen(df_.drop_duplicates(subset=['pert_iname']), 'unique compounds')

    df_['smiles'] = [* map(lambda smiles: smi_largest_component(smi_strip_chemaxon(smiles)), df_['smiles']) ]
    #df_['n_unique_smiles'] = df_.groupby('pert_iname')['smiles'].transform(lambda x: len(set(x)))
    #printlen(df_.query('n_unique_smiles > 1').drop_duplicates(subset=['pert_iname']), 'compounds with non-unique smiles')

    #df_ = df_.drop_duplicates(subset=['pert_iname'], keep='first').reset_index(drop=True)
    #printlen(df_, 'compounds/SMILES after dedup')
    return df_
