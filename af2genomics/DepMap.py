
import ast, collections, datetime, functools, inspect, itertools, math, os, pandas as pd, requests, sqlite3, subprocess, sys, urllib, zipfile
import numpy as np, scipy as sp, scipy.stats, scipy.stats.contingency, matplotlib, matplotlib.pyplot as plt, seaborn as sns
import sklearn as sk, sklearn.decomposition, sklearn.linear_model, sklearn.metrics, sklearn.naive_bayes, sklearn.preprocessing
import tqdm, tqdm.contrib, tqdm.contrib, tqdm.contrib.concurrent

from af2genomics.common import *

def path_(file_):
    #dir_ = os.path.dirname(__file__)
    #sub_ = 'DepMap_Public_21Q3'
    #sub_ = 'DepMap_Public_22Q2'
    #sub_ = 'DepMap_Public_23Q2'
    return workpath(os.path.join('23.07.26_DepMap_Public_23Q2', file_))

@functools.cache
def gene_effect():
    """
    Example:
        # https://depmap.org/portal/gene/EGFR?tab=overview
        # https://cansarblack.icr.ac.uk/target/P00533/synopsis/ligandability
        Achilles.gene_effect()[['EGFR (1956)']].hist()
    """
    return pd.read_csv(path_('CRISPRGeneEffect.csv'), index_col=0)

@functools.cache
def gene_dependency():
    return pd.read_csv(path_('CRISPRGeneDependency.csv'), index_col=0)

@functools.cache
def gene_expression(to_uniprot=False):
    df_ = pd.read_csv(path_('OmicsExpressionProteinCodingGenesTPMLogp1.csv'), index_col=0)
    if not to_uniprot:
        return df_
    else:
        d_ = pd.read_csv('results/23.04_bindfunc/prism_23Q2.depmap_to_af2.tsv', sep='\t').set_index('hugo_entrez')['uniprot_id'].to_dict()
        return resources.DepMap_Public.gene_expression().rename(d_, axis=1)

@functools.cache
def model():
    """
    TODO Achilles.sample_info()['cas9_activity'] looks like a number in most cases but isn't
    Adds cell line doubling times from Li2019: https://www.nature.com/articles/s41591-019-0404-8

    Example:
        Achilles.sample_info().head(3).transpose()
    """
    return pd.read_csv(path_('Model.csv'), index_col=0)

    '''
    if True:#not hasattr(sample_info, 'sample_info_'):
        sample_info_ = pd.read_csv(path_('sample_info.csv'), index_col=0)
        def Li2019_annotations():
            df_ = pd.read_csv('resources/Li2019/1-cell_line_annotations.txt', sep='\t')
            df_.columns = ['name', 'name_with_tissue_origins', 'classifications', 'cell_culture_media', 'gender', 'doubling_time']
            df_['growth_rate'] = math.log(2) / df_['doubling_time'] # https://www.nature.com/scitable/knowledge/library/how-populations-grow-the-exponential-and-logistic-13240157/
            return df_
        #print(len(sample_info_.merge(Li2019_annotations(), left_on='stripped_cell_line_name' , right_on='name')))
        #print(len(sample_info_.merge(Li2019_annotations(), left_on='CCLE_Name' , right_on='name_with_tissue_origins')))
        #df_ = sample_info_.merge(Li2019_annotations(), left_on='stripped_cell_line_name' , right_on='name', validate='one_to_one')
        #df_ = sample_info_.merge(Li2019_annotations(), left_on='CCLE_Name' , right_on='name_with_tissue_origins', validate='one_to_one')

        #right_ = Li2019_annotations().add_prefix('Li2019_')
        # Use reset_index() to keep preserve DepMap_ID values throughout the merge
        #merge_ = sample_info_.reset_index().merge(right_, how='left', left_on='CCLE_Name' , right_on='Li2019_name_with_tissue_origins', validate='one_to_one')
        #sample_info.sample_info_ = merge_.set_index('DepMap_ID')
        sample_info.sample_info_ = sample_info_
    return sample_info.sample_info_
    '''

def model_dummies(col, v=False, *args, **kwargs):
    """
    Example:
        Achilles.get_sample_info_dummies('culture_type')[['Adherent', 'Suspension']]

    Refs:
        https://github.com/EmanuelGoncalves/dtrace/blob/master/dtrace/Associations.py#L132-L149
        https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.get_dummies.html
    """
    covariate = model()[col]
    if v: print(covariate.value_counts(dropna=False))
    return pd.get_dummies(covariate, *args, **kwargs)

@functools.cache
def strongly_selective(min_CRISPR_LRT=100):
    fp_ = workpath('23.07.26_DepMap_Public_23Q2/strongly_selective/full_tda_summary.csv')
    strongly_selective_a_ = set(pd.read_csv(fp_).query('CRISPR_LRT > @min_CRISPR_LRT')['symbol'])
    printlen(strongly_selective_a_, 'genes in full_tda_summary.csv')

    fp_ = workpath('23.07.26_DepMap_Public_23Q2/strongly_selective/Gene Dependency Profile Summary.csv')
    strongly_selective_b_ = set(pd.read_csv(fp_).query('(Dataset == "DependencyEnum.Chronos_Combined") & (`Strongly Selective` == True)')['Gene'])
    printlen(strongly_selective_b_, 'genes in `Gene Dependency Profile Summary.csv`')

    strongly_selective_ = strongly_selective_a_ & strongly_selective_b_
    printlen(strongly_selective_, 'genes in both sets')
    return strongly_selective_

def unpack_hugo_entrez_str(hugo_entrez):
    (hugo, entrez_) = hugo_entrez.split(' ')
    entrez = entrez_.lstrip('(').rstrip(')')
    if '&' in entrez: # DEMETER2 has records of multiple genes, e.g. 399746&440888
        return (hugo, entrez)
    else:
        return (hugo, int(entrez))

def unpack_hugo_entrez(df, col='effect_name'):
    loc = df.columns.get_loc(col)
    (hugo, entrez) = zip(*df[col].map(unpack_hugo_entrez_str))
    df_unpack = df.copy()#.drop(col, axis=1)
    df_unpack.insert(loc=loc+1, column='entrez_id', value=entrez)
    df_unpack.insert(loc=loc+1, column='hugo_symbol', value=hugo)
    return df_unpack

def hugo_entrez():
    return unpack_hugo_entrez(pd.DataFrame({'hugo_entrez': gene_effect().columns}), col='hugo_entrez').astype({'entrez_id': 'int64',})