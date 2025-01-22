
import ast, collections, datetime, functools, inspect, itertools, math, os, pandas as pd, requests, sqlite3, subprocess, sys, urllib, zipfile
import numpy as np, scipy as sp, scipy.stats, scipy.stats.contingency, matplotlib, matplotlib.pyplot as plt, seaborn as sns
import sklearn as sk, sklearn.decomposition, sklearn.linear_model, sklearn.metrics, sklearn.naive_bayes, sklearn.preprocessing
import tqdm, tqdm.contrib, tqdm.contrib, tqdm.contrib.concurrent

from af2genomics.common import *

import af2genomics, af2genomics.ligands

from . import Corsello2017

'''
@functools.cache
def extended_primary_data_matrix(v=True):
    fp_ = path_('Repurposing_Public_23Q2_Extended_Primary_Data_Matrix.csv')
    df_ = pd.read_csv(fp_, index_col=0)
    #m_ = df_.index.str.endswith('_FAILED_STR')
    #if v: print(sum(m_), 'cell lines dropped due to _FAILED_STR')
    #primary_logfold_change.primary_logfold_change_ = df_[~m_]
    return df_
'''

'''
def Corsello2017_read_samples(v=True):
    df_ = pd.read_csv(workpath('23.07.25_Corsello2017/repurposing_samples_20200324.txt'), comment='!', sep='\t').query('qc_incompatible == 0')
    df_['n_unique_smiles'] = df_.groupby('pert_iname')['smiles'].transform(lambda x: len(set(x)))
    if v: print(uf(df_['pert_iname'].nunique()), 'raw compounds')
    df_ = df_.query('n_unique_smiles == 1').reset_index(drop=True)
    if v: print(uf(df_['pert_iname'].nunique()), 'compounds with consistent SMILES across samples')
    return df_ 
'''

@functools.cache
def extended_primary_compound_list(v=True):
    fp_ = workpath('23.07.25_PRISM_Repurposing_23Q2/Repurposing_Public_23Q2_Extended_Primary_Compound_List.csv')
    df_ = pd.read_csv(fp_, index_col=0)
    printlen(df_, 'raw compounds')
    return df_

    df_['Drug.Name'] = df_['Drug.Name'].str.lower()

    df_Corsello2017 = Corsello2017.read_samples()
    df_Corsello2017['broad_id'] = 'BRD:' + df_Corsello2017['broad_id']

    df_ = df_.merge(df_Corsello2017, left_on='IDs', right_on='broad_id', how='left')#.drop(['broad_id', 'pert_iname'], axis=1)
    printlen(df_, 'after merge')

    df_ = df_.query('smiles == smiles')
    printlen(df_, 'after keeping with smiles')
    df_ = df_.query('(150 < expected_mass) & (expected_mass < 1000)')
    printlen(df_, 'after filtering for mass')
    return df_

@functools.cache
def extended_primary_compound_list_adj(v=True):
    df_ = extended_primary_compound_list()
    df_['Drug.Name'] = df_['Drug.Name'].str.lower()

    df_Corsello2017 = resources.Corsello2017.read_samples()[['broad_id', 'pert_iname', 'smiles']]
    df_Corsello2017['broad_id'] = 'BRD:' + df_Corsello2017['broad_id']
    df_merge_ = df_.merge(df_Corsello2017, left_on='IDs', right_on='broad_id', how='left').drop(['broad_id', 'pert_iname'], axis=1)
    return df_merge_

@functools.cache
def extended_primary_data_matrix_adj():
    df_drugs_ = extended_primary_compound_list_adj()[['IDs', 'Drug.Name']].set_index('IDs', drop=True).drop_duplicates().rename({'Drug.Name': 'compound_name'}, axis=1)
    df_ = df_drugs_.merge(extended_primary_data_matrix(), left_index=True, right_index=True).set_index('compound_name', drop=True)
    return df_.transpose()

def read_associations(compound_id):
    fp_ = pfile(compound_id=f'prism_{penc(compound_id)}', step='prism_23Q2', suffix='.h5', base='/cluster/project/beltrao/jjaenes/22.12_pocketomes/results/23.04_bindfunc/')
    df_ = pd.read_hdf(fp_, key='summary')
    df_.insert(loc=0, column='compound_id', value=compound_id)
    return df_

def read_associations_gene_read_(ligand_id, hugo_symbol):
    return read_associations(ligand_id).query('hugo_symbol == @hugo_symbol')

@functools.cache
def read_associations_gene(hugo_symbol):
    l_ = af2genomics.ligands.read_Corsello2020_compounds().compound_id#[:10]
    return pd.concat(list(tqdm.contrib.concurrent.process_map(read_associations_gene_read_, l_, itertools.repeat(hugo_symbol, len(l_)), chunksize=10)), axis=0)

#@functools.cache
#def extended_primary_compound_list_adj(v=True):
#    df_ = extended_primary_compound_list()
#    return df_merge_

'''
@functools.cache
def extended_primary_data_matrix_adj():
    df_drugs_ = extended_primary_compound_list_adj()[['IDs', 'Drug.Name']].set_index('IDs', drop=True).drop_duplicates().rename({'Drug.Name': 'compound_name'}, axis=1)
    df_ = df_drugs_.merge(extended_primary_data_matrix(), left_index=True, right_index=True).set_index('compound_name', drop=True)
    return df_.transpose()

@functools.cache
def get_response(compound_name):
    df_ = extended_primary_data_matrix_adj()[[compound_name]]
    len_ = len(extended_primary_data_matrix_adj())
    print(f'{len_ - df_.isna().sum()} of {len_} non-NA values for {compound_name}')
    return df_

def parse_targets_(s):
    print('s')
    print(s)
    if s != s:
        return []
    return [s_i.strip() for s_i in s.split(',')]

def known_targets(compound_name):
    return parse_targets_(extended_primary_compound_list_adj().query('`Drug.Name` == @compound_name')['repurposing_target'].drop_duplicates().head(1).squeeze()) # head(1) due to small number of possibly inconsistent target annotations, e.g. betaxolol

def prism_23Q2_all():
    return extended_primary_compound_list_adj()['Drug.Name'].dropna().drop_duplicates().tolist()


@functools.lru_cache(maxsize=None)
def primary_treat_info(v=False): #https://stackoverflow.com/questions/279561/what-is-the-python-equivalent-of-static-variables-inside-a-function
    fp_ = path_('primary-screen-replicate-collapsed-treatment-info.csv')
    df_ = pd.read_csv(fp_, index_col=0)
    print(len(df_), f'records in {fp_}')
    #return resources.dup_rows(df_, subset=['name'])
    df_ = df_.drop_duplicates(subset=['name']) # TODO not clear what to do with (a small number of) duplicate records
    print(len(df_), f'records after dropping duplicates')
    assert len(df_) == 4518 # 4,518 drugs from the Drug Repurposing Hub
    return df_

@functools.lru_cache(maxsize=None)
def primary_treat_info_adj(v=False, raw=False):
    fp_ = path_('43018_2019_18_MOESM2_ESM_2_Compounds_tested.txt')
    df_ = pd.read_csv(fp_, sep='\t', skiprows=2)
    assert ~any(df_['Unnamed: 7'] == df_['Unnamed: 7'])
    df_.drop(['Unnamed: 7'], axis=1, inplace=True)
    assert df_['broad_id'].is_unique
    if raw:
        return df_l
    df_ = primary_treat_info().reset_index().merge(df_, left_on='broad_id', right_on='broad_id', suffixes=('', '_SupplTable2')).set_index('column_name')
    df_['target_count'] = [* map(lambda targetstr: len(targetstr.split(',')) if targetstr==targetstr else 0, df_['target'] ) ]
    l_secondary_name = set(secondary_treat_info()['name'])
    if v: print(len(l_secondary_name))
    df_['isin_secondary'] = df_['name'].isin(l_secondary_name)
    return df_

def parse_targets_(s):
    if s != s:
        return []
    return [s_i.strip() for s_i in s.split(',')]

def known_targets(drug_name):
    return parse_targets_(primary_treat_info().query(f'name == "{drug_name}"')['target'].values[0])

def sample_nominal_targets(n=10, skip=[], random_state=4, v=False): #https://xkcd.com/221/
    targets = pd.Series(itertools.chain.from_iterable(map(parse_targets_, primary_treat_info()['target'])))
    if v: print(len(targets), 'targets raw')
    targets = targets.sample(frac=1, random_state=random_state).reset_index(drop=True).drop_duplicates()
    if v: print(len(targets), 'targets after dedup')
    targets_drop = [target for target in targets if target not in skip]
    return targets_drop[:n]

def primary_all(names_only=True):
    df_ = primary_treat_info_adj()
    return df_['name'].dropna().tolist() if names_only else df_

def primary_targeted_launched(names_only=True):
    df_ = primary_treat_info_adj().query('drug_category == "targeted cancer" & phase == "Launched"')
    return df_['name'].tolist() if names_only else df_

def primary_targeted_launched_single(names_only=True):
    df_ = primary_treat_info_adj().query('drug_category == "targeted cancer" & phase == "Launched" & target_count == 1')
    return df_['name'].tolist() if names_only else df_

def both_targeted_launched_known(names_only=True):
    df_ = primary_treat_info_adj().query('drug_category == "targeted cancer" & phase == "Launched" & target_count >= 1 & isin_secondary')
    return df_['name'].tolist() if names_only else df_

def primary_logfold_change(v=True):
    if not hasattr(primary_logfold_change, 'primary_logfold_change_'):
        fp_ = path_('primary-screen-replicate-collapsed-logfold-change.csv')
        df_ = pd.read_csv(fp_, index_col=0)
        m_ = df_.index.str.endswith('_FAILED_STR')
        if v: print(sum(m_), 'cell lines dropped due to _FAILED_STR')
        primary_logfold_change.primary_logfold_change_ = df_[~m_]
    return primary_logfold_change.primary_logfold_change_

def primary_cell_lines():
    return primary_logfold_change().index.tolist()

def primary_logfold_change_by_name():
    return primary_treat_info().merge(primary_logfold_change().transpose(), left_index=True, right_index=True)

def get_primary_logfold_change(drug_name):
    df_ = primary_logfold_change_by_name().query(f'name == "{drug_name}"')[primary_cell_lines()].transpose()
    df_.columns = ['logfold_change']
    return df_

def secondary_treat_info(v=False):
    if not hasattr(secondary_curve_parms, 'secondary_curve_parms_'):
        fp_treat_ = path_('secondary-screen-replicate-collapsed-treatment-info.csv')
        df_treat = pd.read_csv(fp_treat_)
        if v: print(len(df_treat), 'records in secondary-screen-replicate-collapsed-treatment-info.csv')
        if v: print(len(df_treat['name'].unique()), 'unique values in name column (expected=1448)')
        #df_treat.query('name == "PAC-1"')
        #df_treat.head().transpose()

        name_MTS010 = df_treat.query('screen_id == "MTS010"')['name'].unique()
        if v: print(name_MTS010, 'drugs screened in MTS010, expected=147')
        df_treat_rm = df_treat[ ~(df_treat['name'].isin(name_MTS010) & (df_treat['screen_id'] != 'MTS010')) ]
        if v: print(df_treat_rm['name'].unique(), 'unique values in name column')#, expected=1448)
        if v: print('Drugs with more than one screen after choosing MTS010 where available:')
        df_dup_ = df_treat_rm[['name', 'screen_id']].drop_duplicates()['name'].value_counts().loc[lambda n: n > 1]
        if v: print(df_dup_)
        if v: print('Remaining duplicates come from the following screens:')
        if v: print(df_treat_rm[ df_treat_rm['name'].isin(df_dup_.index) ]['screen_id'].value_counts())
        m_ = ~((df_treat_rm['name'].isin(df_dup_.index.tolist())) & (df_treat_rm['screen_id'] == 'HTS002'))
        df_treat_uniq = df_treat_rm[m_].reset_index(drop=True)
        if v: print(df_treat_uniq[['name', 'screen_id']].drop_duplicates(),
               'drug-screen pairs after choosing MTS005, MTS006 over HTS002 (Extended Data Fig. 6)',
               len(df_treat['name'].unique()))

        secondary_curve_parms.secondary_curve_parms_ = df_treat_uniq
    return secondary_curve_parms.secondary_curve_parms_

def secondary_logfold_change(v=False):
    if not hasattr(secondary_logfold_change, 'secondary_logfold_change_'):
        fp_ = path_('secondary-screen-replicate-collapsed-logfold-change.csv')
        df_ = pd.read_csv(fp_, index_col=0)
        m_ = df_.index.str.endswith('_FAILED_STR')
        if v: print(sum(m_), 'cell lines dropped due to _FAILED_STR')
        secondary_logfold_change.secondary_logfold_change_ = df_[~m_]
    return secondary_logfold_change.secondary_logfold_change_

def secondary_curve_parms(v=False):
    if not hasattr(secondary_curve_parms, 'secondary_curve_parms_'):
        fp_ = path_('secondary-screen-dose-response-curve-parameters.csv')
        dtype_ = {'disease.area': 'string', 'indication': 'string'}
        df_ = pd.read_csv(path_(fp_), index_col=0, dtype=dtype_)
        # Remove screens that were later re-done as implied by treat_info
        qc_pass = set([(r['broad_id'], r['screen_id']) for i, r in secondary_treat_info()[['broad_id', 'screen_id']].iterrows()])
        #self.curve_parms['qc_pass'] = self.curve_parms.index.isin(qc_pass)
        df_['qc_pass'] = [* map(lambda broad_id, screen_id: (broad_id, screen_id) in qc_pass, df_.index, df_.screen_id)]
        if v: print(df_['qc_pass'].value_counts())
        if v: print(len(df_), 'before dropping NA in depmap_id')
        df_ = df_.loc[ df_['qc_pass'] ].dropna(subset=['depmap_id'])
        if v: print(len(df_), 'after dropping NA in depmap_id')
        # For now, pragmatically choose a representatives for doxycycline and U-0126 based on number of cell lines with sensitivity information
        if v: print(df_.query('name == "doxycycline"').index.value_counts())
        if v: print(df_.query('name == "U-0126"').index.value_counts())
        df_.drop('BRD-A08545410-311-03-4', inplace=True) # doxycycline -- keep BRD-A08545410-003-07-8 as latter has more cell lines
        df_.drop('BRD-K91701654-001-03-1', inplace=True) # U-0126 -- keep BRD-K18787491-001-08-6 as latter has more cell lines
        if v: print(len(df_), 'after choosing representatives for doxycycline / U-0126')
        secondary_curve_parms.secondary_curve_parms_ = df_
    return secondary_curve_parms.secondary_curve_parms_

def get_secondary_curve_parms(drug_name, sensitivity_col):
    assert sensitivity_col in ['auc', 'ec50', 'ic50'] # sensitivity_col: 'auc', 'ec50', 'ic50'
    df_curve_parms_ = secondary_curve_parms().query(f'name == "{drug_name}"').set_index('depmap_id')[[sensitivity_col]]
    return df_curve_parms_

def bimodality_coefficient(x):
    g = scipy.stats.skew(x.dropna())
    k = scipy.stats.kurtosis(x.dropna())
    n = len(x.dropna())
    return (g**2 + 1) / (k + (3*(n-1)**2) / ((n-2)*(n-3)))

def get_bic_logfold_change(drug_name, n=1, v=True):
    broad_id = primary_logfold_change_by_name().query(f'name == "{drug_name}"')['broad_id'].values[0]
    df_ = pd.concat([secondary_logfold_change(), primary_logfold_change()], axis=1, join='outer').sort_index()
    if v: print(df_.shape)
    m_ = df_.columns.str.startswith(broad_id)#'BRD-K19687926-379-07-4')
    top_screens = df_.loc[:,m_].apply(bimodality_coefficient, axis=0).sort_values(ascending=False).head(n=n).index.values
    if v: print(top_screens)
    return df_[top_screens].mean(axis=1).to_frame()

def adhoc_complex(names_only=True):
    m_ = primary_treat_info_adj()['moa'].str.contains("proteasome inhibitor|mTOR inhibitor|HDAC inhibitor").fillna(False)
    #df_ = primary_treat_info_adj()['moa'].str.contains("proteasome inhibitor|mTOR inhibitor|HDAC inhibitor").fillna(False)
    #df_ = primary_treat_info_adj().query('drug_category == "targeted cancer" & phase == "Launched" & target_count == 1')
    df_ = primary_treat_info_adj()[m_].query('target_count >= 1 & isin_secondary')
    return df_['name'].tolist() if names_only else df_

def primary_with_targets(names_only=True):
    df_ = primary_treat_info_adj().query('target_count >= 1')
    return df_['name'].tolist() if names_only else df_

def secondary_with_targets(names_only=True):
    df_ = primary_treat_info_adj().query('target_count >= 1 & isin_secondary')
    return df_['name'].tolist() if names_only else df_

def secondary_targeted_known(names_only=True):
    #df_ = primary_treat_info_adj().query('drug_category == "targeted cancer" & phase == "Launched" & target_count >= 1 & isin_secondary')
    df_ = primary_treat_info_adj().query('drug_category == "targeted cancer" & target_count >= 1 & isin_secondary')
    return df_['name'].tolist() if names_only else df_

def secondary_targeted_launched_known(names_only=True):
    df_ = primary_treat_info_adj().query('drug_category == "targeted cancer" & phase == "Launched" & target_count >= 1 & isin_secondary')
    return df_['name'].tolist() if names_only else df_

def secondary_noncancer(names_only=True):
    df_ = primary_treat_info_adj().query('drug_category == "noncancer" & isin_secondary')
    return df_['name'].tolist() if names_only else df_

def secondary_nontargets(names_only=True):
    df_ = primary_treat_info_adj().query('isin_secondary & target_count == 0')
    return df_['name'].tolist() if names_only else df_

def secondary_all(names_only=True):
    df_ = primary_treat_info_adj().query('isin_secondary')
    return df_['name'].tolist() if names_only else df_

'''