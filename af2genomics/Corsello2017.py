
import pandas as pd

from af2genomics.common import *

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

def Corsello2017_smiles():
    return read_samples()[['broad_id', 'smiles']].drop_duplicates(keep='first').reset_index(drop=True)

def Corsello2020_smiles():
    df_ = pd.read_excel(workpath('23.07.25_PRISM_Repurposing_23Q2/43018_2019_18_MOESM2_ESM.xlsx'), sheet_name=1, header=2)[['broad_id', 'smiles']]
    df_['smiles'] = df_['smiles'].str.split(',')
    df_ = df_.explode('smiles')
    df_['smiles'] = df_['smiles'].str.lstrip().str.rstrip()
    df_['smiles'] = [* map(lambda smiles: smi_largest_component(smi_strip_chemaxon(smiles)), df_['smiles']) ]
    df_ = df_.drop_duplicates()
    return df_.query('smiles == smiles').reset_index(drop=True)