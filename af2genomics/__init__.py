
import ast, collections, datetime, functools, inspect, itertools, math, os, pandas as pd, requests, sqlite3, subprocess, sys, urllib, zipfile
import numpy as np, scipy as sp, scipy.stats, scipy.stats.contingency, matplotlib, matplotlib.pyplot as plt, seaborn as sns
import sklearn as sk, sklearn.decomposition, sklearn.linear_model, sklearn.metrics, sklearn.naive_bayes, sklearn.preprocessing

from af2genomics.common import *
__all__ = ['alphafold3', 'Corsello2017', 'DepMap', 'foldx', 'huintaf2', 'ligands', 'missense_human', 'missense_yeast', 'PRISM']

try:
    import Bio, Bio.PDB, Bio.SVDSuperimposer, Bio.SeqUtils
except ImportError:
    print('biopython not found; if needed, install with: conda install conda-forge::biopython')

try:
    import prody
except ImportError:
    print('prody not found; if needed, install with: conda install conda-forge::prody')

try:
    import IPython.display
except ImportError:
    print('Cannot import IPython.display')

__all__ += ['RANDOM_SEED']
# Fix `RANDOM_SEED` for (partial) reproducibility
RANDOM_SEED = 4 # https://xkcd.com/221

def format_pct(x):
    return '({:.2f}%)'.format(x)

def printlenq(frame, q, *args, **kwargs):
    n_q = len(frame.query(q))
    n = len(frame)
    f = n_q / n
    print(uf(n_q), 'of', uf(n), format_pct(100*f),  *args, **kwargs)

@functools.cache
def slurm_ntasks():
    try:
        return int(os.environ['SLURM_NTASKS'])
    except:
        return 4 # nproc --all?

def tqdm_map(function, iterable):
    return tqdm.contrib.concurrent.process_map(function, iterable, max_workers=slurm_ntasks(), chunksize=10)

def phead(df, n=3):
    print(uf(len(df)), 'records')
    return df.head(n).transpose()

def pmerge(left, right, *args, **kwargs):
    left_merge_right = left.merge(right, *args, **kwargs)
    print('merge:', uf(len(left_merge_right)), 'left:', uf(len(left)), 'right:', uf(len(right)))
    return left_merge_right

def read_af2_pos(add_residue_id=False, **kwargs):
    kwargs.setdefault('filepath_or_buffer', workpath('23.10.05_VEP_scores/23.10.18_af2_pos.tsv'))
    kwargs.setdefault('sep', '\t')
    df_ = pd.read_csv(**kwargs)
    if add_residue_id:
        df_.insert(loc=1, column='residue_id', value=[ f'{r.uniprot_id}/{r.pos}' for i, r in df_.iterrows() ])
    return df_

def read_human_missense_scores(add_residue_id=False, **kwargs):
    kwargs.setdefault('filepath_or_buffer', workpath('23.10.05_VEP_scores/23.10.13_human_missense_scores.tsv.gz'))
    kwargs.setdefault('dtype', {
        #'position': int, #TypeError: Cannot cast array data from dtype('float64') to dtype('int64') according to the rule 'safe'
        'uniprot_id': str, 
        'protein_variant':str,
        'am_class': str,
        'ESM1b_is_pathogenic': str,
    })
    kwargs.setdefault('sep', '\t')
    kwargs.setdefault('nrows', 10)
    df_ = pd.read_csv(**kwargs)
    if add_residue_id:
        df_.insert(loc=1, column='residue_id', value=df_['variant_id'].str[:-1])
    return df_

def read_interface_full():
    fp_ = workpath('23.09.07_dburke_af2interactions_mutfunc/interface_full.list')
    df_ppi = pd.read_csv(os.popen(f'grep -v "Error" {fp_}'), sep=' ')
    df_ppi['pair_nmodels'] = df_ppi['pair'].str.count('_') + 1
    df_ppi = df_ppi.query('pair_nmodels == 2').copy()
    df_ppi[['pairA', 'pairB']] = df_ppi['pair'].str.split('_', expand=True)
    #print(df_ppi.info())
    return df_ppi

def read_interface_full_strict():
	fp_ = workpath('23.09.07_dburke_af2interactions_mutfunc/interface_full_strict.list')
	df_ppi = pd.read_csv(os.popen(f'grep -v "Error" {fp_}'), sep=' ')
	df_ppi['pair_nmodels'] = df_ppi['pair'].str.count('_') + 1
	df_ppi = df_ppi.query('pair_nmodels == 2').copy()
	df_ppi[['pairA', 'pairB']] = df_ppi['pair'].str.split('_', expand=True)
	return df_ppi

def parse_resid(s):
    """ Return various list-of-residues representations as set-of-ints, e.g.:
    parse_resid('')
    parse_resid('1')
    parse_resid('1,2,3,4')
    parse_resid('{1,2,3,4}')
    parse_resid('[1,2,3,4]')
    """
    if isinstance(s, list) or isinstance(s, np.ndarray):
        return s
    if s == '' or s != s:
        return set()
    x = ast.literal_eval(s)
    if isinstance(x, int):
        return set([x])
    elif isinstance(x[0], int):
        return set(map(int, x))
    else:
        return set(x)

def fillna_set(s):
    # Fill NA-s with an empty set
    # https://stackoverflow.com/questions/33199193/how-to-fill-dataframe-nan-values-with-empty-list-in-pandas
    return s.apply(lambda d: d if isinstance(d, set) else set())

def g_convert(query=["CASQ2", "CASQ1", "GSTO1", "DMD", "GSTM2"], target='UNIPROTSWISSPROT_ACC', organism='hsapiens', numeric_namespace='ENTREZGENE_ACC'):
    """
    Query for HGNC gene names using g:convert (https://biit.cs.ut.ee/gprofiler/convert)
    """
    r = requests.post(url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/', json=locals())
    df_ = pd.DataFrame(r.json()['result'])
    return df_

def g_mapping(left, right, left_on, right_on, organism='hsapiens', target='ENSG', v=False):
    if v: printlen(left, left_on)
    left_d = left[[left_on]].drop_duplicates(keep='first')
    if v: printlen(left_d, 'after dedup')
    left_g = left_d.merge(g_convert(left_d[left_on].tolist(), organism=organism, target=target)[['incoming', 'converted']].rename({'incoming': left_on, 'converted': target}, axis=1), on=left_on)
    if v: printlen(left_g, f'after converting left to {target} including multi-mappings')
    left_g = left_g.query(f'{target} != "None"')
    if v: printlen(left_g, 'after removing empty mappings')

    if v: printlen(right, right_on)
    right_d = right[[right_on]].drop_duplicates(keep='first')
    if v: printlen(right_d, 'after dedup')
    right_g = right_d.merge(g_convert(right_d[right_on].tolist(), organism=organism, target=target)[['incoming', 'converted']].rename({'incoming': right_on, 'converted': target}, axis=1), on=right_on)
    if v: printlen(right_g, f'after converting right {target} including multi-mappings')
    right_g = right_g.query(f'{target} != "None"')
    if v: printlen(right_g, 'after removing empty mappings')

    merge = left_g.merge(right_g, on=target)
    if v: printlen(merge, 'after merge')

    m_dup_ = merge.duplicated(subset=left_on, keep=False) | merge.duplicated(subset=right_on, keep=False)
    merge = merge[~m_dup_]

    n_start = min(len(left_d), len(right_d))
    n_merge = len(merge)
    printlen(merge, f'of {uf(n_start)} unique mappings over {target}')
    return merge

def g_merge(left, right, left_on, right_on, organism='hsapiens', target='ENSG', v=False):
    mapping = g_mapping(left, right, left_on, right_on, organism, target, v)
    merged = left.merge(mapping, on=left_on).merge(right, on=right_on)
    printlen(merged, 'records merging')
    return merged

def merge_gconvert(frame, source_col, target, target_col, how='inner'):
    printlen(frame, 'records')
    vals_ = list(set(frame[source_col]))
    printlen(vals_, 'unique identifiers')
    map_ = g_convert(query=vals_, target=target).set_index('incoming')['converted']
    printlen(vals_, 'mapped identifiers')
    frame_map_ = frame.merge(map_.to_frame(name=target_col), left_on=source_col, right_index=True, how=how)
    printlen(frame_map_, 'after merge')
    return frame_map_

def pquery(df_, *args, **kwargs):
    df_q_ = df_.query(*args, **kwargs)
    printsrc('pquery:', uf(len(df_q_)), 'of', uf(len(df_)), 'rows selected with', args[0])
    return df_q_

def calver(timestamp=None):
    """
    Examples: https://calver.org/#other-notable-projects

    In [4]: calver(datetime.datetime.now())
    Out[4]: '23.10.31_090256'

    In [5]: calver(datetime.datetime.now() + datetime.timedelta(hours=5))
    Out[5]: '23.10.31_140258'
    """
    if timestamp is None:
        timestamp = datetime.datetime.now()
    return timestamp.strftime('%y.%m.%d_%H%M%S')

def read_af2_pLDDT(fp_):
    resseq_pLDDT = collections.OrderedDict()
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(fp_, fp_)
    for chains in structure:
        for chain in chains:
            for residue in chain:
                resname = residue.get_resname()
                hetflag, resseq, icode = residue.get_id()
                for atom in residue:
                    resseq_pLDDT[resseq] = atom.bfactor
    return resseq_pLDDT

def read_swiss_resid():
    df_swiss = pd.read_csv(workpath('23.10.31_SWISS_human_2023-09-14/SWISS-MODEL_Repository/INDEX'), comment='#', sep='\t')
    df_swiss['pos'] = [* map(lambda from_, to_: set(range(from_, to_ + 1)), df_swiss['from'], df_swiss['to']) ]
    df_resid = df_swiss.groupby(['UniProtKB_ac', 'provider']).agg(pos=('pos', lambda x: set.union(*(x_i for x_i in x)))).unstack()
    df_resid.columns = df_resid.columns.droplevel()
    df_resid.rename({'PDB': 'resid_pdb', 'SWISSMODEL': 'resid_swiss'}, axis=1)
    df_resid = df_resid.reset_index().rename({'UniProtKB_ac': 'uniprot_id', 'PDB': 'resid_pdb', 'SWISSMODEL': 'resid_swiss'}, axis=1)
    return df_resid

def nan_to_set(x):
    return x if x == x else set()

def read_Cheng2023_s5():
    # Supplementary Data S5: List of variants and AlphaMissense predictions for the ClinVar benchmark.
    df_ = pd.read_csv(workpath('23.10.05_VEP_scores/science.adg7492_data_s1_to_s9/science.adg7492_data_s5.csv'), sep=',', dtype={'label': int})
    df_['protein_variant'] = df_['protein_variant'].str.replace(':','/')
    return df_

def read_Cheng2023_s6():
    # Supplementary Data S6: List of variants and AlphaMissense predictions for the Cancer hotspot mutations benchmark.
    df_ = pd.read_csv(workpath('23.10.05_VEP_scores/science.adg7492_data_s1_to_s9/science.adg7492_data_s6.csv'), sep=',')
    df_['protein_variant'] = df_['protein_variant'].str.replace(':','/')
    return df_

def read_Cheng2023_s7():
    # Supplementary Data S7: List of variants and AlphaMissense predictions for the Deciphering Developmental Disorders benchmark.
    df_ = pd.read_csv(workpath('23.10.05_VEP_scores/science.adg7492_data_s1_to_s9/science.adg7492_data_s7.csv'), sep=',', dtype={'label': int})
    df_['protein_variant'] = df_['protein_variant'].str.replace(':','/')
    return df_

def read_Cheng2023_s8():
    # Supplementary Data S8: List of variants and AlphaMissense predictions for the ProteinGym benchmark.
    df_ = pd.read_csv(workpath('23.10.05_VEP_scores/science.adg7492_data_s1_to_s9/Supplementary_Data_S8_proteingym.csv'), sep=',')#, dtype={'label': int})
    #df_['protein_variant'] = df_['protein_variant'].str.replace(':','/')
    return df_

def read_structures():
    df_ = pd.read_csv(workpath('23.11.01_human_protein_map/structures_23.11.1.tsv'), sep='\t')
    df_['struct_resid'] = df_['n_resid'].map(lambda n_resid: set(range(1, n_resid + 1)))
    return df_

def read_missense():
    with sqlite3.connect(workpath('23.11.01_human_protein_map/missense_23.11.1.sqlite')) as db:
        df_ = pd.read_sql_query(sql='SELECT * FROM missense', con=db)
    return df_

def query_missense(variants):
    #https://stackoverflow.com/questions/28735213/pandas-read-sql-with-a-list-of-values-for-where-condition
    with sqlite3.connect(workpath('23.11.01_human_protein_map/missense_23.11.1.sqlite')) as db:
        df_ = pd.read_sql_query(sql='SELECT * FROM missense WHERE variant_id in ' + str(tuple(variants)), con=db)
    return df_.replace({
        'pred_ddg': -99,
        'pocketscore': -99,
        'pocketrank': -99,
        'interface': -99,
        'interface_strict': -99,
        'freq': -99,
    }, np.nan)

def query_missense_uniprot_id(uniprot_id):
    with sqlite3.connect(workpath('23.11.01_human_protein_map/missense_23.11.1.sqlite')) as db:
        df_ = pd.read_sql_query(sql=f'SELECT * FROM missense WHERE variant_id GLOB "{uniprot_id}*"', con=db)
    return df_

def query_missense_yeast(variants):
    with sqlite3.connect(workpath('24.08.30_yeast_missense/yeast_missense.sqlite')) as db:
        df_ = pd.read_sql_query(sql='SELECT * FROM missense WHERE variant_id in ' + str(tuple(variants)), con=db)
    return df_.set_index('variant_id').loc[variants].reset_index()

def read_pockets_resid(pocket_score):
    # Pockets, one per 
    #Q96EK7-F1      157
    cols_ = ['uniprot_id', 'struct_id', 'pocket_id', 'pocket_resid', 'pocket_score_combined_scaled']
    q_ = 'pocket_score_combined_scaled >= @pocket_score'
    return read_pockets().query(q_).explode('pocket_resid')[cols_].sort_values('pocket_score_combined_scaled', ascending=False).groupby(['uniprot_id', 'pocket_resid'])\
        .head(1)#.query('struct_id == "Q96EK7-F1" & pocket_resid == 157')

'''
@functools.cache
def read_ppi_resid(pdockq):
    df_models = read_af2_human_interactions(pdockq=pdockq)
    cols_ = ['uniprot_id', 'ifresid', 'pdockq']
    q_ne_ = 'protein1 != protein2'
    q_eq_ = 'protein1 == protein2'
    df_interfaces = pd.concat([
        df_models.query(q_ne_).rename({'protein1': 'uniprot_id', 'residues1': 'ifresid',}, axis=1)[cols_],
        df_models.query(q_eq_).rename({'protein1': 'uniprot_id', 'residues1': 'ifresid',}, axis=1)[cols_],
        df_models.query(q_ne_).rename({'protein2': 'uniprot_id', 'residues2': 'ifresid',}, axis=1)[cols_],
    ], axis=0)
    df_interfaces['ifresid'] = df_interfaces['ifresid'].map(parse_resid)
    df_ifresid = df_interfaces.explode('ifresid').sort_values('pdockq', ascending=False).groupby(['uniprot_id', 'ifresid']).head(1)#.drop_duplicates(keep='first')
    return df_ifresid.sort_values(['uniprot_id', 'ifresid']).rename({'pdockq': 'ifresid_pdockq',}, axis=1)
'''

@functools.cache
def suppl_ppi_models(pdockq, ifresid1_col=None, ifresid2_col=None):
    df_models = pd.read_csv(workpath('24.05.17_suppl/suppl_ppi_models.tsv.gz'), sep='\t', keep_default_na=False)
    df_models = df_models.query('sources != "negatome"')
    df_models = df_models.query('pdockq > @pdockq')
    if not(ifresid1_col is None) and not(ifresid2_col is None):
        fp_ifresid_ = workpath('23.12.06_ppi_reselect/af2-models-split-ifresid.tsv.gz')
        df_ifresid_ = pd.read_csv(fp_ifresid_, sep='\t', usecols=['interaction_id', ifresid1_col, ifresid2_col], nrows=None, keep_default_na=False).rename({ifresid1_col: 'ifresid1', ifresid2_col: 'ifresid2'}, axis=1)
        cols_ = df_models.columns.tolist()
        df_models = df_models.drop(['ifresid1', 'ifresid2'], axis=1).merge(df_ifresid_, on='interaction_id', how='left')[cols_]
    return df_models.reset_index(drop=True)

def suppl_ppi_residues(pdockq):
    df_models = suppl_ppi_models(pdockq)
    cols_ = ['uniprot_id', 'ifresid', 'pdockq']
    q_ne_ = 'uniprot_id1 != uniprot_id2'
    q_eq_ = 'uniprot_id1 == uniprot_id2'
    df_ppi_residues = pd.concat([
        df_models.query(q_ne_).rename({'uniprot_id1': 'uniprot_id', 'ifresid1': 'ifresid',}, axis=1)[cols_],
        df_models.query(q_eq_).rename({'uniprot_id1': 'uniprot_id', 'ifresid1': 'ifresid',}, axis=1)[cols_],
        df_models.query(q_ne_).rename({'uniprot_id2': 'uniprot_id', 'ifresid2': 'ifresid',}, axis=1)[cols_],
    ], axis=0)
    df_ppi_residues['ifresid'] = df_ppi_residues['ifresid'].map(lambda l: set(int(r[1:]) for r in l.split(',') if r != ''))
    df_ppi_residues = df_ppi_residues.explode('ifresid').sort_values('pdockq', ascending=False).groupby(['uniprot_id', 'ifresid']).head(1)
    return df_ppi_residues.sort_values(['uniprot_id', 'ifresid']).rename({'pdockq': 'ifresid_pdockq',}, axis=1)

def merge_missense(frame, variant_col, pdockq=.23, pocket_score=800, *args, **kwargs):
    cols_out = frame.columns.tolist()
    printlen(frame, 'raw records')
    frame[['_uniprot_id', '_aa_pos', '_aa_ref', '_aa_alt']] = frame.apply(lambda r: parse_varstr(r[variant_col]), axis=1, result_type='expand')

    variant_frame = query_missense(frame[variant_col]).set_index('variant_id')[['am_pathogenicity', 'am_class', 'pred_ddg']]
    frame = frame.merge(variant_frame, left_on=variant_col, right_index=True, how='left', *args, **kwargs)
    printlen(frame, 'records matched to predictions')

    frame['am_label'] = frame['am_class'] == 'pathogenic'
    frame['pred_ddg_label'] = frame['pred_ddg'].map(lambda pred_ddg: pred_ddg > 2)
    printlen(frame.query('pred_ddg_label'), 'annotated as destabilizing')
    cols_out += ['am_pathogenicity', 'am_class', 'am_label', 'pred_ddg', 'pred_ddg_label']

    #frame = frame.merge(read_ppi_resid(pdockq=pdockq), left_on=['_uniprot_id', '_aa_pos'], right_on=['uniprot_id', 'ifresid'], suffixes=('', '_ppi'), how='left').rename({'ifresid_pdockq': 'interface_pdockq'}, axis=1)
    frame = frame.merge(suppl_ppi_residues(pdockq=pdockq), left_on=['_uniprot_id', '_aa_pos'], right_on=['uniprot_id', 'ifresid'], suffixes=('', '_ppi'), how='left').rename({'ifresid_pdockq': 'interface_pdockq'}, axis=1)
    frame['interface_label'] = frame['ifresid'].map(lambda ifresid: ifresid == ifresid)
    printlen(frame.query('interface_label'), 'annotated with interfaces at pDockQ >', pdockq)
    cols_out += ['interface_pdockq', 'interface_label']

    frame = frame.merge(read_pockets_resid(pocket_score=pocket_score), left_on=['_uniprot_id', '_aa_pos'], right_on=['uniprot_id', 'pocket_resid'], how='left')
    frame['pocket_label'] = frame['pocket_score_combined_scaled'].map(lambda pocket_score_combined_scaled: pocket_score_combined_scaled == pocket_score_combined_scaled)
    printlen(frame.query('pocket_label'), 'annotated with pockets at pocket_score >', pocket_score)
    cols_out += ['pocket_label']

    def mechanistic_label_(r):
        if r.pred_ddg != r.pred_ddg:
            return 'Unassigned'
        elif r.pred_ddg_label:
            return 'Stability'
        elif r.interface_label:
            return 'Interface'
        elif r.pocket_label:
            return 'Pockets'
        else:
            return 'Unassigned'
    frame['mechanistic_label'] = frame.apply(mechanistic_label_, axis=1)
    cols_out += ['mechanistic_label']

    return frame[cols_out]

def pvalue_label(pvalue):
    # https://www.graphpad.com/support/faq/what-is-the-meaning-of--or--or--in-reports-of-statistical-significance-from-prism-or-instat/
    if pvalue <= .0001:
        return 'p≤0.0001'
    elif pvalue <= .001:
        return 'p≤0.001'
    elif pvalue <= .01:
        return 'p≤0.01'
    elif pvalue <= .05:
        return 'p≤0.05'
    else:
        return f'p={pvalue:.2f}'

def plot_merge_missense(frame, col):
    stats = pd.DataFrame.from_records([
        sp.stats.fisher_exact(pd.crosstab(frame[col], frame['am_label']).values),
        sp.stats.fisher_exact(pd.crosstab(frame[col], frame['pred_ddg_label']).values),
        sp.stats.fisher_exact(pd.crosstab(frame[col], frame['pocket_label']).values),
        sp.stats.fisher_exact(pd.crosstab(frame[col], frame['interface_label']).values),
    ], index=['Pathogenicity', 'Stability', 'Pockets', 'Interfaces'], columns=['odds_ratio', 'pvalue'])
    stats.index.name = 'predictor'
    stats['pvalue_label'] = stats['pvalue'].map(pvalue_label)
    ax = sns.barplot(data=stats.reset_index(), y='predictor', x='odds_ratio', color='tab:blue')
    ax.bar_label(ax.containers[0], labels=stats['pvalue_label'])
    return stats

def jaccard(a, b):
    return len(a & b) / len(a | b)

def interaction_id(uniprot_id_1, uniprot_id_2):
    # sorted pair of uniprot_id-s to compare interacting pairs from different sources
    return '_'.join(sorted([uniprot_id_1, uniprot_id_2]))

@functools.cache
def af2_uniprot_id():
    # AF2 human single-fragment structures (set used for pocket detection, interface modelling, etc)
    return set(read_structures()['uniprot_id'])

@functools.cache
def read_tissue_coabundance():
    fp_ = workpath('24.01.17_coabundances/DLT240108_tissue_coabundance_for_interfaces_af2_human_24.01.08.csv')
    df_ = pd.read_csv(fp_).drop(['prot1', 'prot2', 'protein1', 'protein2'], axis=1)
    cols_ = ['blood', 'brain', 'breast', 'colon', 'kidney', 'liver', 'lung', 'ovary', 'pancreas', 'stomach', 'throat']
    df_['co_abundance'] = df_[cols_].mean(axis=1)
    return df_

def flatten(l):
    return [item for sublist in l for item in sublist]

def read_summary_source(summary=False):
    df_ = pd.read_csv(workpath('23.10.02_dburke_kcl_OneDrive/summary_source.out.bz2'), na_values={'pdockq_fd': 'chain'}, nrows=None)\
        .rename({'#protdir': 'protdir'}, axis=1)
    printlen(df_, 'raw models')
    def parse_(s_):
        l_id_ = s_.split('/')[1].split('_')
        if len(l_id_) == 2:
            return interaction_id(*l_id_)
        else:
            return '.'
    df_.insert(0, 'interaction_id', [* df_['protdir'].map(parse_) ])
    df_.insert(2, 'folding_method', df_['protdir'].map(lambda s: s.split('/')[-1]))
    df_ = df_.query('interaction_id != "."').copy()
    df_[['uniprot_id_A', 'uniprot_id_B']] = df_['interaction_id'].str.split('_', expand=True)
    print(uf(len(df_)), '\t', uf(len(set(df_['interaction_id']))), 'models (interactions) after discarding non-dimers')
    df_ = df_.query('(uniprot_id_A in @af2_uniprot_id()) & (uniprot_id_B in @af2_uniprot_id())').copy()
    print(uf(len(df_)), '\t', uf(len(set(df_['interaction_id']))), 'models (interactions) after keeping human AF2 single-fragment')
    if summary:
        df_ = df_[['interaction_id', 'source']].groupby('interaction_id').agg(
            source = ('source', lambda x: ':'.join(sorted(set(flatten([ x_i.split(':') for x_i in x ] ))))),
        ).rename({'source': 'summary_source'}, axis=1)
        #.drop_duplicates(subset=['interaction_id'])\
        #.set_index('interaction_id').rename({'source': 'summary_source'}, axis=1).copy()
        df_['isin_summary_source'] = True
        #printlen(df_, 'after naive de-duplication with potential loss of information (!!)')
        return df_[['isin_summary_source', 'summary_source']]
    else:
        return df_

def read_ppi_reselect(add_source=True):
    #fp_ = workpath('23.12.06_ppi_reselect/interface_best2_p10.csv') #8A
    #fp_ = workpath('23.12.06_ppi_reselect/interface_relaxed_best2_p10.csv') # 10A
    fp_ = workpath('23.12.06_ppi_reselect/interface_strict_best2_p10.csv') #5A
    df_ = pd.read_csv(fp_, delim_whitespace=True).rename({
        'chain': 'chain1',
        'residues': 'residues1',
        'chain': 'chain1',
        'chain.1': 'chain2',
        'residues.1': 'residues2',
    }, axis=1)
    printlen(df_, 'raw records from', fp_)

    def parse_(s_):
        l_id_ = s_.split('_')
        if len(l_id_) == 2:
            return interaction_id(*l_id_)
        else:
            return '.'
    df_.insert(0, 'interaction_id', [* df_['pair'].map(parse_) ])
    df_ = df_.query('interaction_id != "."').copy()
    printlen(df_, 'after discarding non-dimers')
    df_[['pair1', 'pair2']] = df_['pair'].str.split('_', expand=True)
    #return df_
    df_ = df_.query('(protein1 in @af2_uniprot_id()) & (protein2 in @af2_uniprot_id())').copy()
    printlen(df_, 'after keeping uniprot_id-s in AF2 single fragment structures')
    df_ = df_.drop_duplicates(subset=['interaction_id'], keep='first')
    printlen(df_, 'after naive de-duplication of interaction_id')
    printlen(df_.query('pdockq > .23'), 'with pdockq > .23')
    printlen(df_.query('pdockq > .5'), 'with pdockq > .5')
    df_['pdb'] = df_.apply(lambda r: workpath(f'23.12.06_ppi_reselect/af2-models-split/{r.pair.split("_")[0]}/{r.pair}.pdb'), axis=1)
    if add_source:
        df_source = read_summary_source(summary=True)[['summary_source']]
        df_ = df_.merge(df_source, left_on='interaction_id', right_index=True)
        printlen(df_, 'after adding evidence from summary_source.out.bz2')
    return df_.reset_index(drop=True)

@functools.cache
def read_pockets():
    #fp_ = workpath('23.10.16_af2_human_pockets/23.10.16_af2_human_pocket_summary.tsv.gz')
    fp_ = workpath('23.10.16_af2_human_pockets/af2.obabel_hxr.autosite.summary.combined_score.swiss_coverage.tsv')
    df_ = pd.read_csv(fp_, sep='\t')
    df_.insert(loc=0, column='uniprot_id', value=df_['struct_id'].map(strip_fragment_id))
    df_['pocket_resid'] = df_['pocket_resid'].map(parse_resid)
    labels_ = ['0 to 200', '200 to 400', '400 to 600', '600 to 800', '800 to 1000']
    df_['pocket_score_bin'] = pd.qcut(df_['pocket_score_combined_scaled'], q=5, labels=labels_)
    return df_

__all__.append('pLDDT_bins')
pLDDT_bins = [0, 50, 70, 90, 100]

__all__.append('pLDDT_palette')
# pLDDT visualisation to match AFDB, e.g. https://alphafold.ebi.ac.uk/entry/P00533
pLDDT_palette = {
    'Very low (pLDDT ≤ 50)': '#ff7d45ff',
    'Low (70 ≥ pLDDT > 50)': '#ffdb13ff',
    'High (90 ≥ pLDDT > 70)': '#64cbf3ff',
    'Very high (pLDDT > 90)': '#0053d6ff',
}

def to_pymol_init():
    print("""delete all
bg_color white
set_color tab_blue, [31,119,180]
set_color tab_orange, [255,127,14]
set_color tab_green, [44,160,44]
set_color tab_red, [214,39,40]
set_color tab_purple, [148,103,189]
set_color tab_brown, [140,86,75]
set_color tab_pink, [227,119,194]
set_color tab_gray, [127,127,127]
set_color tab_olive, [188,189,34]
set_color tab_cyan, [23,190,207]""")

def to_pymol_struct(struct_id, base='~/euler-home/project/22.12_pocketomes/results/23.04_bindfunc'):
    print('load', pfile(struct_id=struct_id, step='af2', suffix='.pdb', base=base))
    print('color grey50,', struct_id)

def to_pymol_pocket(struct_id, pocket_id, transparency=0.5, color='tab_blue'):
    pocket_cl_file = read_pockets().query('struct_id == @struct_id & pocket_id == @pocket_id').squeeze().pocket_cl_file
    pocket_name = os.path.splitext(os.path.basename(pocket_cl_file))[0]
    pocket_path = os.path.join('~/euler-home/project/22.12_pocketomes/', pocket_cl_file)
    print(f'load {pocket_path}')
    print(f'hide everything, {pocket_name}')
    print(f'show surface, {pocket_name}')
    print(f'set transparency, {transparency}, {pocket_name}')
    print(f'color {color}, {pocket_name}')

def pymol_resid(s):
    return '+'.join(map(str, sorted(s)))

def read_chain_CA(fp_, chain_):
    # get CA coordinate matrix of a specific chain
    structure_ = Bio.PDB.PDBParser(QUIET=True).get_structure(fp_, fp_)
    return np.array([residue['CA'].coord for residue in structure_[0][chain_]]) #https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec203

"""
def rmsd_(r1_, r2_):
    # align & calculate rmsd for two structures
    coord1_ = read_chain_CA(r1_.pdb, r1_.bait_chain)
    coord2_ = read_chain_CA(r2_.pdb, r2_.bait_chain)
    #print(coord1_.shape)
    #print(coord2_.shape)

    sup = Bio.SVDSuperimposer.SVDSuperimposer()
    sup.set(coord1_, coord2_)
    sup.run()
    print(sup.get_init_rms(), sup.get_rms())
"""
def bait_rmsd(fpA, fpB):
    # RMSD of (shared) chain
    #fpA = bait_.df_bait_id.iloc[0].pdb
    #fpB = bait_.df_bait_id.iloc[9].pdb
    pdbA = prody.parsePDB(fpA)
    pdbB = prody.parsePDB(fpB)
    pdbAm, pdbBm, seqid, overlap = prody.matchChains(pdbA, pdbB)[0]
    rmsd_pre_ = prody.calcRMSD(pdbAm, pdbBm)
    pdbAm, transformation = prody.superpose(pdbAm, pdbBm)
    rmsd_post_ = prody.calcRMSD(pdbAm, pdbBm)
    #print(rmsd_pre_, rmsd_post_)
    return rmsd_post_

def read_chain_len(fp_, chain_):
    # number of residues in a chain
    structure_ = Bio.PDB.PDBParser(QUIET=True).get_structure(fp_, fp_)
    return len([residue['CA'].coord for residue in structure_[0][chain_]]) #https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec203

def parse_varstr(s):
    # df_var[['uniprot_id', 'aa_pos', 'aa_ref', 'aa_alt']] = df_var.apply(lambda r: parse_varstr(r['protein_variant']), axis=1, result_type='expand')
    uniprot_id, variant_id = s.split('/')
    aa_pos = int(variant_id[1:-1])
    aa_ref = variant_id[0]
    aa_alt = variant_id[-1]
    #print(uniprot_id, aa_pos, aa_ref, aa_alt)
    return uniprot_id, aa_pos, aa_ref, aa_alt

def parse_variant_id(variant_id):
    # df_var[['uniprot_id', 'aa_pos', 'aa_ref', 'aa_alt']] = df_var.apply(lambda r: parse_varstr(r['protein_variant']), axis=1, result_type='expand')
    aa_pos = int(variant_id[1:-1])
    aa_ref = variant_id[0]
    aa_alt = variant_id[-1]
    #print(uniprot_id, aa_pos, aa_ref, aa_alt)
    return aa_pos, aa_ref, aa_alt

__all__.append('VariantEnrichment')
class VariantEnrichment:
    # struct_data, struct_uniprot_id, struct_resid
    # var_data, var_uniprot_id, var_resid
    def __init__(self, struct_data, struct_uniprot_id, struct_resid, struct_bin, var_data, var_uniprot_id, var_resid):
        self.struct_data = struct_data[[struct_uniprot_id, struct_resid, struct_bin]].copy().rename({struct_uniprot_id: 'struct_uniprot_id', struct_resid: 'struct_resid', struct_bin: 'struct_bin'}, axis=1)
        self.struct_data['struct_resid'] = self.struct_data['struct_resid'].map(parse_resid)

        # Attach set of all residues
        #df_af2 = pd.read_csv('results/23.04_bindfunc/af2.tsv', sep='\t').query('n_frags == 1')[['uniprot_id', 'n_resid']]
        df_af2 = read_structures()
        df_af2['resid'] = df_af2['n_resid'].map(lambda n_resid: set(range(1, n_resid + 1)))
        self.struct_data = self.struct_data.merge(df_af2[['uniprot_id', 'resid']], left_on='struct_uniprot_id', right_on='uniprot_id').drop(['uniprot_id'], axis=1)

        self.var_data = var_data[[var_uniprot_id, var_resid]].copy().rename({var_uniprot_id: 'var_uniprot_id', var_resid: 'var_resid'}, axis=1)
        #self.var_data['var_resid'] = self.var_data['var_resid'].map(parse_resid)

        # merge has to start with the universe of all proteins!!
        self.merge_data = self.struct_data.merge(self.var_data, left_on='struct_uniprot_id', right_on='var_uniprot_id', how='left')
        def fillna_(x): return x if x == x else set()
        self.merge_data['var_uniprot_id'] = self.merge_data['var_uniprot_id'].map(fillna_)
        self.merge_data['var_resid'] = self.merge_data['var_resid'].map(fillna_)

        self.merge_data['n_struct_var'] = self.merge_data.apply(lambda r: len(r['struct_resid'] & r['var_resid']), axis=1)
        self.merge_data['n_struct_novar'] = self.merge_data.apply(lambda r: len(r['struct_resid'] - r['var_resid']), axis=1)
        self.merge_data['n_nostruct_var'] = self.merge_data.apply(lambda r: len((r['resid'] - r['struct_resid']) & r['var_resid']), axis=1)
        self.merge_data['n_nostruct_novar'] = self.merge_data.apply(lambda r: len((r['resid'] - r['struct_resid']) - r['var_resid']), axis=1)

    def fisher_exact_(k, l, m, n):
        print('new2')
        return sp.stats.fisher_exact([[k, l], [m, n]])

    def fisher_exact(self):
        self.merge_data_sums_ = self.merge_data[['n_struct_var', 'n_struct_novar', 'n_nostruct_var', 'n_nostruct_novar']].sum(axis=0)
        return VariantEnrichment.fisher_exact_(
            k=self.merge_data_sums_.n_struct_var,
            l=self.merge_data_sums_.n_nostruct_var,
            m=self.merge_data_sums_.n_struct_novar,
            n=self.merge_data_sums_.n_nostruct_novar,
        )

    def fisher_exact_groupby(self):
        sums_ = self.merge_data.groupby('struct_bin')[['n_struct_var', 'n_struct_novar', 'n_nostruct_var', 'n_nostruct_novar']].sum()
        sums_[['statistic', 'pvalue']] = sums_.apply(lambda r: VariantEnrichment.fisher_exact_(r.n_struct_var, r.n_nostruct_var, r.n_struct_novar, r.n_nostruct_novar), axis=1, result_type='expand')
        return sums_

def fisher_var_resid(resid, var, all):
    def agg_(s):
        return s.groupby(s.index).agg(lambda x: set.union(*x))

    var_ = agg_(var).rename('var')
    resid_ = agg_(resid).rename('resid')
    all_ = agg_(all).rename('all')

    left_ = pd.merge(all_, var_, left_index=True, right_index=True, how='left')
    merge = pd.merge(left_, resid_, left_index=True, right_index=True, how='left')

    def fillna_(x): return x if x == x else set()
    merge['var'] = merge['var'].map(fillna_)
    merge['resid'] = merge['resid'].map(fillna_)
    merge['all'] = merge['all'].map(fillna_)

    merge['n_resid_var'] = merge.apply(lambda r: len(r['resid'] & r['var']), axis=1)
    merge['n_resid_novar'] = merge.apply(lambda r: len(r['resid'] - r['var']), axis=1)
    merge['n_noresid_var'] = merge.apply(lambda r: len((r['all'] - r['resid']) & r['var']), axis=1)
    merge['n_noresid_novar'] = merge.apply(lambda r: len((r['all'] - r['resid']) - r['var']), axis=1)

    sums = merge[['n_resid_var', 'n_resid_novar', 'n_noresid_var', 'n_noresid_novar']].sum(axis=0)

    def fisher_exact_(k, l, m, n):
        return sp.stats.fisher_exact([[k, l], [m, n]])
    
    def odds_ratio_(k, l, m, n):
        return scipy.stats.contingency.odds_ratio([[k, l], [m, n]], kind='conditional').statistic
    
    (statistic, pvalue) = fisher_exact_(sums.n_resid_var, sums.n_resid_novar, sums.n_noresid_var, sums.n_noresid_novar)
    return pd.Series([statistic, pvalue], index=['odds_ratio', 'pvalue'])
    #odds_ratio = odds_ratio_(sums.n_resid_var, sums.n_resid_novar, sums.n_noresid_var, sums.n_noresid_novar)
    #return pd.Series([statistic, pvalue, odds_ratio], index=['sample_odds_ratio', 'pvalue', 'conditional_odds_ratio'])

def run_pymol(cmds):
    s_ = subprocess.run('module load gcc/6.3.0 pymol; pymol -cpQ', input='\n'.join(cmds), text=True, shell=True, capture_output=True)
    return s_

__all__.append('SharedInterfaceAnalysis')
class SharedInterfaceAnalysis():
    """
    def resid_updated_(r):
    res_A, res_B = calc_interface_residues(r.pdb, dist_threshold=5, plddt_threshold=0)
    if r.bait_chain == 'A':
        return res_A
    else:
        return res_B
    bait_.df_interactors['bait_ifresid'] = [ resid_updated_(r) for i, r in bait_.df_interactors.iterrows() ]
    """
    def __init__(self, df_interactors):
        self.uniprot_id_bait = df_interactors.head(1)['bait_id'].squeeze()
        #print(self.uniprot_id_bait, 'bait uniprot id')
        self.df_interactors = df_interactors.copy()

        pdbs_ = self.df_interactors['pdb'].copy()
        self.df_interactors['pdb'] = pdbs_.map(lambda pdb: os.path.join('/cluster/work/beltrao/jjaenes/23.12.06_ppi_reselect/af2-models-split', pdb))
        self.df_interactors['pdb_trim_bf'] = pdbs_.map(lambda pdb: os.path.join('/cluster/work/beltrao/jjaenes/23.12.06_ppi_reselect/af2-models-split-trim_bf', pdb))

        #printlen(self.df_interactors, 'interaction models sharing the bait')
        self.nresid_bait = read_chain_len(self.df_interactors.head(1).pdb.squeeze(), self.df_interactors.head(1).bait_chain.squeeze())
        #print(self.nresid_bait, 'number of bait residues')
        def parse_(s):
            return set([int(r[1:]) for r in s.split(',')]) if s != '' else set()
        self.df_interactors['bait_resseq'] = self.df_interactors['bait_ifresid'].map(parse_)
        self.df_interactors['n_bait_resseq'] = self.df_interactors['bait_resseq'].map(len)
        self.df_interactors = self.df_interactors.query('n_bait_resseq >= 10').reset_index(drop=True)

    def build_matrix(self):
        def apply_(r):
            return [resid in r.bait_resseq for resid in range(1, self.nresid_bait + 1)]
        self.df_matrix = pd.DataFrame(self.df_interactors.apply(apply_, axis=1).to_list(), index=self.df_interactors['interactor_id'], columns=range(1, self.nresid_bait + 1))
        #print(self.df_matrix.shape, 'shape matrix dimensions')
        #print(uf(sum(self.df_matrix)), 'non-zero entries in residue matrix')

    def linkage(self, criterion='distance', t=.9):
        self.linkage_ = scipy.cluster.hierarchy.linkage(self.df_matrix, method='average', metric='jaccard')
        #self.df_interactors['labels'] = scipy.cluster.hierarchy.fcluster(Z=self.linkage_, criterion='distance', t=.95)
        #self.df_interactors['labels'] = scipy.cluster.hierarchy.fcluster(Z=self.linkage_, criterion='maxclust', t=2)
        self.df_interactors['labels'] = scipy.cluster.hierarchy.fcluster(Z=self.linkage_, criterion=criterion, t=t)
        #print(self.df_interactors['labels'].value_counts())

    #def plot_clustermap_default(self):
    #    sns.clustermap(self.df_matrix, metric='jaccard', row_cluster=True, col_cluster=False, cbar_pos=None, figsize=(8,4))

    def clustermap(self, fname=None, title=None):
        row_colors_ = [* self.df_interactors['labels'].map(lambda label: matplotlib.colormaps['tab10'].colors[label - 1]) ]
        if not (fname is None):
            plt.ioff()

        self.clustermap_ = sns.clustermap(self.df_matrix, row_linkage=self.linkage_, col_cluster=False, cbar_pos=None, figsize=(8,4), row_colors=row_colors_)
        if title is None:
            plt.title(self.uniprot_id_bait)
        else:
            plt.title(title)
        if not (fname is None):
            plt.savefig(fname, bbox_inches='tight', transparent=True)
            plt.clf()
            plt.ion()

    def to_pymol(self, fname=None):
        """
            fname is:
                None => show pymol commands to execute locally
                file name => execute pymol on the cluster & store session as .pse
        """
        df_interactors_top = self.df_interactors.sort_values('pdockq', ascending=False).groupby('labels').head(1) # Align to top interaction model
        #df_interactors_top = self.df_interactors.sort_values(['labels', 'pdockq'], ascending=False)#.head(1) # Align to top interaction model
        cmd_ = []
        cmd_.append('delete all')
        cmd_.append('bg_color white')
        ref_ = df_interactors_top.head(1).squeeze() # Reference structure to align against
        for i, r in df_interactors_top.iterrows():
            #if fname is None:
            #    #fp_ = os.path.join('~/work-euler/', r.pdb.removeprefix('/cluster/work/beltrao/jjaenes/'))
            #    #fp_ = os.path.join('~/work-euler/', r.pdb.removeprefix('/cluster/work/beltrao/jjaenes/'))
            #else:
            #    fp_ = r.pdb

            print(f'load {r.pdb_trim_bf}')
            cmd_.append(f'load {r.pdb_trim_bf}')

            if r.interaction_id != ref_['interaction_id']:
                #cmd_.append(f'align {r.interaction_id} & chain {r.bait_chain}, {ref_.interaction_id} & chain {ref_.bait_chain}')
                #chain_ = 'A' if r.interaction_id.startswith(r.bait_id) else 'B'
                #chain_ref_ = 'A' if ref_.interaction_id.startswith(ref_.bait_id) else 'B'

                cmd_.append(f'align {r.interaction_id} & chain {r.bait_chain}, {ref_.interaction_id} & chain {ref_.bait_chain}')

            col_ = '0x' + matplotlib.colors.to_hex(matplotlib.colormaps['tab10'].colors[r.labels - 1])[1:]
            interactor_chain = 'B' if r.bait_chain == 'A' else 'A'
            cmd_.append(f'color gray, {r.interaction_id} & chain {r.bait_chain}')
            cmd_.append(f'color {col_}, {r.interaction_id} & chain {interactor_chain}')

            #if r.interaction_id.startswith(ref_.bait_id):
            #    cmd_.append(f'color gray, {r.interaction_id} & chain A')
            #    cmd_.append(f'color {col_}, {r.interaction_id} & chain B')
            #else:
            #    cmd_.append(f'color {col_}, {r.interaction_id} & chain A')
            #    cmd_.append(f'color gray, {r.interaction_id} & chain B')

        for label, df_label in df_interactors_top.groupby('labels'):
            name = f'cluster{label}'
            members = ' '.join(df_label['interaction_id'])
            cmd_.append(f'group {name}, {members}')

        if fname is None:
            print('\n'.join(cmd_))
        else:
            cmd_.append(f'save {fname}')
            run_pymol(cmd_)

# Adapted from: https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/protein_contact_map/
def calc_residue_dist(residue_one, residue_two):
	"""Returns the C-alpha distance between two residues (C-alpha)"""
	diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
	return np.sqrt(np.sum(diff_vector * diff_vector))

def calc_residue_dist_all(residue_one, residue_two):
	"""Returns the C-alpha distance between two residues (all-atom)"""
	min_dist = float('inf')
	for atom_one in residue_one:
		for atom_two in residue_two:
			diff_vector = atom_one.coord - atom_two.coord
			min_dist = min(min_dist, np.sum(diff_vector * diff_vector))
	return np.sqrt(min_dist)

def calc_dist_matrix(chain_one, chain_two):
	"""Returns a matrix of C-alpha distances between two chains"""
	answer = np.zeros((len(chain_one), len(chain_two)), float)
	for row, residue_one in enumerate(chain_one) :
		for col, residue_two in enumerate(chain_two) :
			#answer[row, col] = calc_residue_dist(residue_one, residue_two)
			answer[row, col] = calc_residue_dist_all(residue_one, residue_two)
	return answer

def calc_plddt_matrix(chain_one, chain_two):
    """Returns a matrix of C-alpha distances between two chains"""
    answer = np.zeros((len(chain_one), len(chain_two)), float)
    for row, residue_one in enumerate(chain_one):
        for col, residue_two in enumerate(chain_two):
            for atom in residue_one:
                bfactor1 = atom.bfactor
            for atom in residue_two:
                bfactor2 = atom.bfactor
            answer[row, col] = min(bfactor1, bfactor2)
    return answer

def calc_interface_residues(fname, dist_threshold=5.0, plddt_threshold=70):
    """
    res_A, res_B, dist_, plddt_ = calc_interface_residues('/cluster/work/beltrao/jjaenes/23.12.06_ppi_reselect/af2-models-split/P62306/P62306_Q9Y333.pdb', plddt_threshold=70)
    """
    # all-to-all with dist_threshold 5, 8, 10 re-produces read_ppi_reselect()
    #pdb_code = 'P06730_Q9NRG4'
    #pdb_filename = 'P06730_Q9NRG4.pdb'
    structure = Bio.PDB.PDBParser(QUIET=True).get_structure(fname, fname)
    model = structure[0]

    dist_matrix = calc_dist_matrix(model['A'], model['B'])
    plddt_matrix = calc_plddt_matrix(model['A'], model['B'])

    contact_map = (dist_matrix <= dist_threshold)
    contact_map = (dist_matrix <= dist_threshold) & (plddt_matrix > plddt_threshold)

    #print("Minimum distance", numpy.min(dist_matrix))
    #print("Maximum distance", numpy.max(dist_matrix))
    res_A = np.squeeze(np.where(np.any(contact_map, axis=1))) + 1
    res_B = np.squeeze(np.where(np.any(contact_map, axis=0))) + 1
    #print('A:', res_A)
    #print('B:', res_B)
    return res_A, res_B#, dist_matrix, plddt_matrix

def get_chain_seq(fname, target_chain):
    seq = ''
    parser = Bio.PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(fname, fname)
    for chains in structure:
        for chain in chains:
            if chain.id == target_chain:
                for residue in chain:
                    seq += residue.get_resname()
    return Bio.SeqUtils.seq1(seq)

def read_biogrid(pairs_only=False, summary=False, counts=False):
    with zipfile.ZipFile(workpath('23.11.21_ppi_evidence/biogrid/BIOGRID-ALL-4.4.227.tab3.zip'), 'r') as zf:
        df_ = pd.read_csv(zf.open('BIOGRID-ALL-4.4.227.tab3.txt'), sep='\t',
            dtype={
                'Entrez Gene Interactor A': str,
                'Entrez Gene Interactor B': str,
            },
            na_values={
                'Score': '-',
            },
            #nrows=1000
        )
    printlen(df_, 'raw records')
    df_ = df_.query('(`SWISS-PROT Accessions Interactor A` in @af2_uniprot_id()) & (`SWISS-PROT Accessions Interactor B` in @af2_uniprot_id())')
    printlen(df_, 'mapped to AF2 human single-fragment structures')

    df_ = df_.query('`Experimental System Type` == "physical"').copy()
    printlen(df_, 'with `Experimental System Type` == "physical"')

    if pairs_only:
        pairs_ = set(interaction_id(r['SWISS-PROT Accessions Interactor A'], r['SWISS-PROT Accessions Interactor B']) for i, r in df_.iterrows())
        printlen(pairs_, 'unique pairs')
        return pairs_
    elif summary:
        df_['interactor_id'] = [ interaction_id(r['SWISS-PROT Accessions Interactor A'], r['SWISS-PROT Accessions Interactor B']) for i, r in df_.iterrows() ]
        df_['isin_biogrid'] = True
        df_ = df_[['interactor_id', 'isin_biogrid']].drop_duplicates()
        printlen(df_, 'after filtering for unique interactions')
        return df_.set_index('interactor_id')
    elif counts:
        df_['interaction_id'] = [ interaction_id(r['SWISS-PROT Accessions Interactor A'], r['SWISS-PROT Accessions Interactor B']) for i, r in df_.iterrows() ]
        return df_.sort_values('interaction_id')[['interaction_id', 'Publication Source']].groupby('interaction_id').size().to_frame(name='n_references')
    else:
        return df_

def read_string(combined_score_min=400, pairs_only=False, summary=False, scores=False):
    df_alias_ = pd.read_csv(workpath('23.11.21_ppi_evidence/string/9606.protein.aliases.v12.0.txt.gz'), sep='\t')\
        .query('source == "UniProt_AC" & alias in @af2_uniprot_id()')\
        .rename({'#string_protein_id': 'string_protein_id', 'alias': 'uniprot_id'}, axis=1)[['string_protein_id', 'uniprot_id']]
    
    #fp_ = workpath('23.11.21_ppi_evidence/string/9606.protein.links.v12.0.txt.gz')
    fp_ = workpath('23.11.21_ppi_evidence/string/9606.protein.physical.links.v12.0.txt.gz')
    df_ = pd.read_csv(fp_, delim_whitespace=True)\
        .merge(df_alias_, left_on='protein1', right_on='string_protein_id').rename({'uniprot_id': 'protein1_uniprot_id'}, axis=1).drop('string_protein_id', axis=1)\
        .merge(df_alias_, left_on='protein2', right_on='string_protein_id').rename({'uniprot_id': 'protein2_uniprot_id'}, axis=1).drop('string_protein_id', axis=1)\
        [['protein1_uniprot_id', 'protein2_uniprot_id', 'combined_score', 'protein1', 'protein2']]
    printlen(df_, 'raw records from', fp_)
    df_ = df_.query('combined_score >= @combined_score_min').copy()
    printlen(df_, 'after filtering for min combined score')
    df_['interaction_id'] = [ interaction_id(r['protein1_uniprot_id'], r['protein2_uniprot_id']) for i, r in df_.iterrows() ]
    printlen(df_, 'after keeping uniprot_id-s in AF2 single fragment structures and filtering for combined_score_min')
    if pairs_only:
        printlen(df_['interactor_id'], 'unique pairs')
        return set(df_['interactor_id'])
    elif summary:
        df_ = df_[['interaction_id', 'combined_score']].groupby('interaction_id').agg(combined_score=('combined_score', np.max))
        df_['isin_string'] = True
        printlen(df_, 'after aggregating directionality')
        return df_[['isin_string']]
    elif scores:
        return df_[['interaction_id', 'combined_score']].groupby('interaction_id').agg(string_combined_score=('combined_score', np.max))
    else:
        return df_

def read_af2_human_interactions(drop_negatome=True, pdockq=.5):
    fp_ = workpath('24.01.30_af2_human_interations/24.01.30_af2_human_interactions.tsv')
    df_ = pd.read_csv(fp_, sep='\t')
    printlen(df_, 'raw records from', fp_)
    if drop_negatome:
        df_ = df_.query('source != "negatome"')
        printlen(df_, 'after removing negatome-only interactions')
    if not (pdockq is None):
        df_ = df_.query('pdockq > @pdockq')
        printlen(df_, 'after filtering for pdockq')
    return df_.reset_index(drop=True)

def read_evidence():
    df_evidence = pd.read_csv(workpath('23.11.01_human_protein_map/uniprot_evidence_23.05.1.tsv'), sep='\t')
    #df_struct.query('resid_pdb != resid_pdb & resid_swiss != resid_swiss')
    #df_struct.merge(df_evidence, left_on='uniprot_id', right_on='accession'#, how='left')
    #m_ = df_evidence.duplicated(subset='accession', keep=False) 
    #df_evidence[m_] #.groupby('accession')['accession'].value_counts()
    return df_evidence.drop_duplicates(subset='accession').reset_index(drop=True)

def strip_fragment_id(af2_id):
    # A0A385XJ53-F1	=> A0A385XJ53
    fid_ = af2_id[::-1].split('-', maxsplit=1)[0][::-1]
    assert fid_[0] == 'F'
    return af2_id[:-(len(fid_) + 1)]

@functools.cache
def read_clinvar(*args, add_clinvar_label=True, drop_conflicting=True, **kwargs):
    kwargs.setdefault('dtype', {
        'CHROM': str,
        'ONC': str,
        'ONCREVSTAT': str,
        'SCI': str,
        'SCIREVSTAT': str,
        'CLNDISDBINCL': str,
        'CLNDNINCL': str,
        'CLNSIGINCL': str,
        'DBVARID': str,
        'SCIDISDB': str,
        'SCIDN': str,
        'ONCDISDB': str,
        'ONCDN': str,
        'Chromosome': str,
    })
    fp_ = workpath('23.06.02_clinvar/24.04.22_protvar_out/clinvar_mapped.tsv')
    df_ = pd.read_csv(fp_, sep='\t', *args, **kwargs)
    printlen(df_, 'rows from', fp_)
    def clinvar_label_(r):
        if r.CLNSIG in ['Benign', 'Likely_benign', 'Benign/Likely_benign']:
            return 'Benign'
        elif r.CLNSIG in ['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic']:
            return 'Pathogenic'
        elif r.CLNSIG in ['Uncertain_significance']:
            return 'VUS'
        else:
            return 'Conflicting/other'

    if add_clinvar_label:
        df_['clinvar_label'] = df_.apply(clinvar_label_, axis=1)
    
    if drop_conflicting:
        df_ = df_.query('clinvar_label != "Conflicting/other"').reset_index(drop=True)
        printlen(df_, 'after removing conflicting/other variants')
    
    return df_

def rgyr(fp_):
    # rgyr is corelated to protein size: https://doi.org/10.1134/S0026893308040195
    # https://github.com/biopython/biopython/blob/master/Bio/PDB/Entity.py#L298-L329
    # https://adkgromacstutorial.readthedocs.io/en/latest/analysis/rgyr.html
    # https://github.com/sarisabban/Rg/blob/master/Rg.py
    struct = Bio.PDB.PDBParser().get_structure('_', fp_)
    com_coord = struct.center_of_mass()

    df_ = pd.DataFrame.from_records([(atom.coord[0] - com_coord[0], atom.coord[1] - com_coord[1], atom.coord[2] - com_coord[2], atom.mass)\
        for atom in struct.get_atoms()
    ], columns=['dx', 'dy', 'dz', 'mass'])
    df_['mass'] = df_['mass'] / df_['mass'].sum()
    return math.sqrt(np.dot((df_['dx']**2 + df_['dy']**2 + df_['dz']**2), df_['mass']))

def read_autosite_cl(fp):
    # PDBParser does not grok AutoSite output due to duplicate atom names; address by reading manually using read_fwf, and constructing Atom objects adhoc
    #clust = Bio.PDB.PDBParser(QUIET=True).get_structure(fp_cl, fp_cl)
    # https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
    # https://gist.github.com/tbrittoborges/929ced78855945f3e296
    colspecs = [(0, 6), (6, 11), (12, 16), (16, 17), (17, 20), (21, 22), (22, 26), (26, 27), (30, 38),
                (38, 46), (46, 54), (54, 60), (60, 66), (76, 78), (78, 80)]
    names = ['ATOM', 'Atom serial number', 'Atom name', 'Alternate location indicator',
             'Residue name', 'Chain identifier', 'Residue sequence number', 'Insertion code',
             'X', 'Y', 'Z', 'Occupancy', 'Temperature factor', 'Segment identifier', 'Element symbol']
    df_ = pd.read_fwf(fp, names=names, colspecs=colspecs)
    return df_

def resid_center_of_mass(fp, resid):
    CA_x, CA_y, CA_z = [], [], []
    struct = Bio.PDB.PDBParser(QUIET=True).get_structure(fp, fp)
    for chains in struct:
        for chain in chains:
            for residue in chain:
                resname = residue.get_resname()
                hetflag, resseq, icode = residue.get_id() # http://biopython.org/DIST/docs/api/Bio.PDB.Residue-pysrc.html#Residue.__repr__
                if resseq in resid:
                    for atom in residue:
                        if atom.get_name() == 'CA':
                            (atom_x, atom_y, atom_z) = atom.get_vector()
                            CA_x.append(atom_x)
                            CA_y.append(atom_y)
                            CA_z.append(atom_z)
    return np.mean(CA_x), np.mean(CA_y), np.mean(CA_z)

'''
def read_protvarmap_pharmgkb():
    fp_protvarmap = workpath('23.07.03_pharmgkb/24.02.08_pharmgkb_protvar/a2178b7a-a8cd-489f-9264-bc0529287fd7.csv')
    df_protvarmap = pd.read_csv(fp_protvarmap, sep=',')
    printlen(df_protvarmap, 'unique positions in pharmgkb')
    print(df_protvarmap['Consequences'].value_counts())
    df_protvarmap = pd.read_csv(fp_protvarmap, sep=',').query('Consequences == "missense"')
    printlen(df_protvarmap, 'mapped to a missense variant by ProtVar')
    #df_protvarmap.columns
    #cols_ = ['User_input', 'Chromosome', 'Coordinate', 'ID', 'Reference_allele',
    #    'Alternative_allele', 'Notes', 'Gene', 'Codon_change', 'Strand',
    #    'CADD_phred_like_score', 'Canonical_isoform_transcripts',
    #    'Uniprot_canonical_isoform_(non_canonical)',
    #    'Alternative_isoform_mappings', 'Protein_name', 'Amino_acid_position',
    #    'Amino_acid_change', 'Consequences',
    #]
    mapper_ = {
        'Uniprot_canonical_isoform_(non_canonical)': 'uniprot_id',
        'Amino_acid_position': 'resid_pos',
        'Amino_acid_change': 'resid_change',
        'User_input': 'Variant/Haplotypes',
    }
    df_protvarmap = df_protvarmap[mapper_.keys()].rename(mapper_, axis=1).astype({'resid_pos': int}).reset_index(drop=True)
    printlen(df_protvarmap, 'after cleaning up columns')

    df_protvarmap = df_protvarmap.groupby(['uniprot_id', 'resid_pos', 'Variant/Haplotypes']).agg({'resid_change': lambda x: ';'.join(x)}).reset_index()
    printlen(df_protvarmap, 'after aggregating codon changes')
    return df_protvarmap
'''

def read_pharmgkb():
    # Drop unsupported annotations: https://www.pharmgkb.org/page/clinAnnLevels
    evidence_ = set(['1A', '1B', '2A', '2B', '3'])
    df_ = pd.read_csv(workpath('23.07.03_pharmgkb/24.04.24_protvar.tsv'), sep='\t')
    return df_.query('`Level of Evidence` in @evidence_')

@functools.cache
def get_uniprot_sequence(accession):
    url = f'https://rest.uniprot.org/uniprotkb/search?query={accession}&fields=sequence'
    return requests.get(url).json()['results'][0]['sequence']['value']

def read_am_aa_substitutions(*args, **kwargs):
    fp_ = workpath('23.10.05_VEP_scores/dm_alphamissense/AlphaMissense_aa_substitutions.tsv.gz')
    df_ = pd.read_csv(fp_, sep='\t', comment='#', *args, **kwargs)#, index_col=['uniprot_id', 'protein_variant']
    df_.insert(loc=0, column='variant_id', value=[*map(lambda uniprot_id, protein_variant: f'{uniprot_id}/{protein_variant}', df_['uniprot_id'], df_['protein_variant'])])
    df_[['aa_pos', 'aa_ref', 'aa_alt']] = df_.apply(lambda r: parse_variant_id(r['protein_variant']), axis=1, result_type='expand')
    return df_

def any_check(x):
    # Check that all elements in iterator are identical & return any/first one
    assert len(set(iter(x))) == 1
    return next(iter(x))

def all_dup(df, *args, **kwargs):
    return df[df.duplicated(keep=False, *args, **kwargs)]

def print_jaccard(a, b):
    print(len(a & b) / len(a | b), a - b, b - a)

def get_AFmultimer_ifresid(file, min_distance=8, min_pLDDT=0, include_resname=True):
    """
    Identify interface residues in AF-Multimer structures by filtering on all-atom intra-chain distance (min_distance), and residue-level model confidence (min_pLDDT)
    """
    parser = Bio.PDB.PDBParser(QUIET=True)
    struct = parser.get_structure(file, file)
    try:
        chain1, chain2 = struct[0].get_chains()
    except ValueError:
        return list(), list()
    # Chain 1/2 atoms, filtered by min_pLDDT
    chain1_atoms = list(filter(lambda a: a.get_bfactor() >= min_pLDDT, Bio.PDB.Selection.unfold_entities(struct[0][chain1.id], 'A')))
    chain2_atoms = list(filter(lambda a: a.get_bfactor() >= min_pLDDT, Bio.PDB.Selection.unfold_entities(struct[0][chain2.id], 'A')))
    if len(chain1_atoms) == 0 or len(chain2_atoms) == 0:
        return list(), list()
    # Atoms in chain1/2 within radius of the other chain
    chain1_ns = Bio.PDB.NeighborSearch(chain1_atoms)
    chain2_ns = Bio.PDB.NeighborSearch(chain2_atoms)
    chain1_ifatoms = list(filter(lambda atom: len(chain2_ns.search(center=atom.get_coord(), radius=min_distance, level='A')) > 0, chain1_atoms))
    chain2_ifatoms = list(filter(lambda atom: len(chain1_ns.search(center=atom.get_coord(), radius=min_distance, level='A')) > 0, chain2_atoms))

    def get_resseq(atom):
        return atom.get_parent().get_id()[1]

    def get_resname(atom):
        resname3 = str(atom.get_parent().get_resname()).capitalize()
        return Bio.Data.IUPACData.protein_letters_3to1[resname3]

    # Return resseq positions of the corresponding interface residues
    #def get_resseqs(atoms): return sorted(set(((atom.get_parent().get_id()[1], atom.get_parent().get_resname()) for atom in atoms)))
    def get_resseqs(atoms): return sorted(set((atom.get_parent().get_id()[1] for atom in atoms)))
    def get_resseqs_resnames(atoms): return [ f'{resname}{resseq}' for (resseq, resname) in sorted(set(((get_resseq(atom), get_resname(atom)) for atom in atoms))) ]
    if include_resname:
        return get_resseqs_resnames(chain1_ifatoms), get_resseqs_resnames(chain2_ifatoms)
    else:
        return get_resseqs(chain1_ifatoms), get_resseqs(chain2_ifatoms)

def dispall(frame, max_rows=100, max_columns=None, max_colwidth=None):
    with pd.option_context('display.max_rows', max_rows, 'display.max_columns', max_columns, 'display.max_colwidth', max_colwidth):
        IPython.display.display(frame)

@functools.cache
def read_intact_mutations():
    df_ = pd.read_csv(workpath('23.11.21_ppi_evidence/intact_mutations/mutations.tsv'), sep='\t', engine='python').rename({'#Feature AC': 'Feature AC'}, axis=1)
    printlen(df_, 'IntAct mutations')
    df_ = df_.query('`Interaction participants` == `Interaction participants`').copy()
    printlen(df_, 'with defined `Interaction participants`')
    df_['Interaction participants'] = df_['Interaction participants'].str.split('|')
    df_['n_interaction_participants'] = [* map(len, df_['Interaction participants']) ]
    df_ = df_.query('`n_interaction_participants` == 2').copy().reset_index(drop=True)
    printlen(df_, 'with exactly two participants')

    df_participants_ = pd.DataFrame(df_['Interaction participants'].to_list(), columns = ['Participant A', 'Participant B'])
    df_ = pd.concat([df_, df_participants_], axis=1)
    printlen(df_, 'after parsing participants')

    prefix = 'uniprotkb:'
    suffix = '(protein(MI:0326), 9606 - Homo sapiens)'
    df_ = df_.query('`Participant A`.str.startswith(@prefix) & `Participant A`.str.endswith(@suffix)').copy().reset_index(drop=True)
    printlen(df_, 'after filtering participant A for human')
    df_ = df_.query('`Participant B`.str.startswith(@prefix) & `Participant B`.str.endswith(@suffix)').copy().reset_index(drop=True)
    printlen(df_, 'after filtering participant B for human')

    df_['Participant A'] = df_['Participant A'].str.removeprefix(prefix).str.removesuffix(suffix)
    df_['Participant B'] = df_['Participant B'].str.removeprefix(prefix).str.removesuffix(suffix)

    df_[['Feature range(s) A', 'Feature range(s) B']] = df_['Feature range(s)'].str.split('-', expand=True)
    df_ = df_.query('(`Feature range(s) A` == `Feature range(s) B`) & (`Resulting sequence`.str.len() == 1)').copy().reset_index(drop=True)
    printlen(df_, 'single-residue substititions')
    df_['Feature range(s) A'] = pd.to_numeric(df_['Feature range(s) A'], errors='coerce')

    labels_ = { #https://www.ebi.ac.uk/intact/download/datasets#mutations
        'mutation with no effect(MI:2226)': 'neutral',
        'mutation disrupting(MI:0573)': 'disrupt/decrease',
        'mutation disrupting strength(MI:1128)': 'disrupt/decrease',
        'mutation disrupting rate(MI:1129)': 'disrupt/decrease',
        'mutation decreasing(MI:0119)': 'disrupt/decrease',
        'mutation decreasing strength(MI:1133)': 'disrupt/decrease',
        'mutation decreasing rate(MI:1130)': 'disrupt/decrease',
        'mutation causing(MI:2227)': 'cause/increase',
        'mutation increasing(MI:0382)': 'cause/increase',
        'mutation increasing strength(MI:1132)': 'cause/increase',
        'mutation increasing rate(MI:1131)': 'cause/increase',
        'mutation(MI:0118)': 'no data on effect',
    }
    df_['effect'] = df_['Feature type'].map(labels_)
    printlen(df_, 'effects -', ' '.join(f'{k}: {v}' for k, v in df_['effect'].value_counts().items()))

    df_['interaction_id'] = [ interaction_id(r['Participant A'], r['Participant B']) for i, r in df_.iterrows() ]
    df_['resid_id'] = [ f'{r["Original sequence"]}{r["Feature range(s) A"]}' for i, r in df_.iterrows() ]

    def apply_(r):
        return r['Participant A'] + '/' + r['Original sequence'] + str(r['Feature range(s) A']) + r['Resulting sequence']
    df_['variant_id'] = df_.apply(apply_, axis=1)

    cols_ = collections.OrderedDict([
        ('variant_id', 'variant_id'),
        ('Participant A', 'uniprot_id1'),
        ('Feature range(s) A', 'pos'),
        ('Original sequence', 'resid'),
        ('resid_id', 'resid_id'),
        ('Participant B', 'uniprot_id2'),
        ('interaction_id', 'interaction_id'),
        ('effect', 'effect'),
    ])
    #return df_
    df_ = df_[cols_.keys()].rename(cols_, axis=1).reset_index(drop=True)
    df_ = df_.drop_duplicates(keep='first').reset_index(drop=True)
    printlen(df_, 'after de-duplication')
    return df_

def get_AFmultimer_resid(file):
    parser = Bio.PDB.PDBParser(QUIET=True)
    struct = parser.get_structure(file, file)
    try:
        chain1, chain2 = struct[0].get_chains()
    except ValueError:
        return list(), list()

    chain1_resid = list(Bio.PDB.Selection.unfold_entities(struct[0][chain1.id], 'R'))
    chain2_resid = list(Bio.PDB.Selection.unfold_entities(struct[0][chain2.id], 'R'))

    def get_resseq(resid):
        return resid.get_id()[1]

    def get_resname(resid):
        resname3 = str(resid.get_resname()).capitalize()
        return Bio.Data.IUPACData.protein_letters_3to1[resname3]

    def get_resseqs_resnames(resids): return [ f'{resname}{resseq}' for (resseq, resname) in sorted(set(((get_resseq(resid), get_resname(resid)) for resid in resids))) ]
    return get_resseqs_resnames(chain1_resid), get_resseqs_resnames(chain2_resid)

def aa_3to1(resname3):
    return Bio.Data.IUPACData.protein_letters_3to1[resname3.capitalize()]

def read_protvar(fp_):
    df_ = pd.read_csv(fp_)
    printlen(df_, 'raw records')
    df_ = df_.query('Amino_acid_position == Amino_acid_position & Consequences == "missense"').copy()
    printlen(df_, 'missense mappings')

    def apply_(r):
        uniprot_id = r['Uniprot_canonical_isoform_(non_canonical)']
        aa_pos = int(r.Amino_acid_position)
        aa_ref, aa_alt = r.Amino_acid_change.split('/')
        return f'{uniprot_id}/{aa_3to1(aa_ref)}{aa_pos}{aa_3to1(aa_alt)}'

    df_['variant_id'] = df_.apply(apply_, axis=1)
    df_ = df_[['variant_id', 'User_input']].drop_duplicates(keep='first')
    printlen(df_, 'after dedup')
    return df_.reset_index(drop=True)

class ROC:
    def __init__(self):
        self.l_y_true = []
        self.l_y_score = []
        self.l_label = []

    def add(self, y_true, y_score, label):
        self.l_y_true.append(y_true)
        self.l_y_score.append(y_score)
        self.l_label.append(label)
        
    def barh(self, max_fpr=None, yticklabels=True, *args, **kwargs):
        l_auc = []
        for y_true, y_score, label in zip(self.l_y_true, self.l_y_score, self.l_label):
            auc = sk.metrics.roc_auc_score(y_true=y_true, y_score=y_score, max_fpr=max_fpr)
            l_auc.append(auc)
            print('barh:', label, auc)

        #plt.barh(range(len(l_auc)), l_auc, *args, **kwargs)
        df_auc = pd.DataFrame({'label': self.l_label, 'auc': l_auc})
        sns.barplot(data=df_auc, x='auc', y='label', *args, **kwargs)

        if max_fpr is None:
            plt.gca().set_xlabel('AUC')
        else:
            plt.gca().set_xlabel(f'partial AUC at FDR={max_fpr}')

        plt.gca().set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 1.0], minor=False)
        #plt.gca().set_xlim(0.5, 1)
        plt.gca().set_xlim(.48, 1.02)
        plt.gca().xaxis.grid(True)

        plt.gca().set_yticks(range(len(l_auc)))
        if yticklabels:
            plt.gca().set_yticklabels(self.l_label)
        else:
            plt.gca().set_yticklabels(['' for label in self.l_label])
        #plt.gca().invert_yaxis()

    def roc(self, *args, **kwargs):
        l_auc = []
        for y_true, y_score, label in zip(self.l_y_true, self.l_y_score, self.l_label):
            fpr, tpr, thresholds = sk.metrics.roc_curve(y_true=y_true, y_score=y_score)
            auc = sk.metrics.roc_auc_score(y_true=y_true, y_score=y_score)
            plt.plot(fpr, tpr, label=f'{label}, AUC={auc:.2f}')

        plt.gca().set_xlim(-.02, 1.02)
        plt.gca().set_ylim(-.02, 1.02)
        plt.gca().set_xlabel('FPR')
        plt.gca().set_ylabel('TPR')
        plt.gca().set_aspect('equal')
        plt.plot([0, 1], [0, 1], color='gray')
        plt.legend(*args, **kwargs)

    def aucs(self, max_fpr=None):
        l_auc = []
        for y_true, y_score, label in zip(self.l_y_true, self.l_y_score, self.l_label):
            auc = sk.metrics.roc_auc_score(y_true=y_true, y_score=y_score, max_fpr=max_fpr)
            l_auc.append(auc)
            #print(label, auc)

        s_ = pd.Series(data=l_auc, index=self.l_label, name='auc')
        s_.index.name='metric'
        return s_

def edc(pdb_file, disease_resseq, min_pLDDT=70):
    resseq_pLDDT = collections.OrderedDict()
    parser = Bio.PDB.PDBParser(QUIET=True)
    struct = parser.get_structure(pdb_file, pdb_file)

    def is_CA(a):
        return a.get_id() == 'CA'
    def is_disease(a):
        return a.get_parent().get_id()[1] in disease_resseq
    def is_pLDDT(a):
        return a.get_bfactor() >= min_pLDDT

    filtered_atoms = list(filter(lambda a: is_CA(a), Bio.PDB.Selection.unfold_entities(struct[0], 'A')))
    disease_atoms = list(filter(lambda a: is_disease(a), filtered_atoms))
    healthy_atoms = list(filter(lambda a: not is_disease(a), filtered_atoms))
    #print(len(filtered_atoms), len(disease_atoms), len(healthy_atoms), 'filtered/disease/healthy atoms')

    def kth(k, l):
        return l[k], l[:k] + l[(k + 1):]

    def kth_disease_mindist(k):
        atm, etc = kth(k, disease_atoms)
        return min(atm - atm_etc for atm_etc in etc)

    disease_mindist = [* map(kth_disease_mindist, range(len(disease_atoms))) ]
    healthy_mindist = [ min(healthy_atom - disease_atom for disease_atom in disease_atoms) for healthy_atom in healthy_atoms ]
    return np.mean(np.log(healthy_mindist)) / np.mean(np.log(disease_mindist))

__all__.extend([name for (name, thing) in locals().items() if callable(thing)]) #https://stackoverflow.com/questions/18451541/getting-a-list-of-locally-defined-functions-in-python
