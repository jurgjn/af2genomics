
import ast, collections, datetime, functools, inspect, itertools, math, os, pandas as pd, requests, sqlite3, subprocess, sys, zipfile
import numpy as np, scipy as sp, scipy.stats, scipy.stats.contingency, matplotlib, matplotlib.pyplot as plt, seaborn as sns

try:
    import Bio, Bio.PDB, Bio.SVDSuperimposer, Bio.SeqUtils
except ImportError:
    print('biopython not found; if needed, install with: conda install conda-forge::biopython')

try:
    import prody
except ImportError:
    print('prody not found; if needed, install with: conda install conda-forge::prody')

__all__ = ['RANDOM_SEED']
# Fix `RANDOM_SEED` for (partial) reproducibility
RANDOM_SEED = 4 # https://xkcd.com/221

def projectpath(path):
    #dir_ = os.path.dirname(__file__)
    dir_ = '/cluster/project/beltrao/jjaenes'
    return os.path.join(dir_, path)

def workpath(path):
    #dir_ = os.path.dirname(__file__)
    dir_ = '/cluster/work/beltrao/jjaenes'
    return os.path.join(dir_, path)

def uf(x):
    return '{:,}'.format(x)

def ul(x):
    return uf(len(x))

def printsrc(*args, **kwargs):
    """
        https://stackoverflow.com/questions/3056048/filename-and-line-number-of-python-script
        https://stackoverflow.com/questions/3711184/how-to-use-inspect-to-get-the-callers-info-from-callee-in-python
        https://github.com/snakemake/snakemake/blob/main/snakemake/exceptions.py#L17
    """
    #pprint(dir(inspect.currentframe().f_back))
    #pprint(dir(inspect.getframeinfo(inspect.currentframe().f_back)))
    frameinfo_ = inspect.getframeinfo(inspect.currentframe().f_back)
    #pprint(frameinfo_)
    #pprint(dir(frameinfo_))
    filename = frameinfo_.filename
    lineno = frameinfo_.lineno
    #lineno = workflow.linemaps[filename][ frameinfo_.lineno ]
    print(f'{os.path.basename(filename)}:{lineno}', *args, **kwargs)

def printlen(x, *args, **kwargs):
    name_ = inspect.stack()[1][3] #https://stackoverflow.com/questions/5067604/determine-function-name-from-within-that-function-without-using-traceback
    if name_ != '<module>':
        print(f'{name_}:', uf(len(x)), *args, **kwargs)
    else:
        print(uf(len(x)), *args, **kwargs)

def penc(s):
    """Encode exotic characters within identifiers using percent-encoding, e.g.:
        BF2.649 => BF2%2e649
        UH-232-(+) => UH-232-%28%2b%29
        bis(maltolato)oxovanadium(IV) => bis%28maltolato%29oxovanadium%28IV%29
    """
    #return urllib.parse.quote_plus(s) # Does not decode dots, e.g. BF2.649
    def penc_char(c):
        if c.isalnum() or c in ['_', '-']:
            return c
        else:
            return "%{0:0>2}".format(format(ord(c), "x"))

    return "".join(map(penc_char, s))

def pdec(s):
    """Opposite of pdec(). Sanity check:
    for id_ in ['BF2.649', 'UH-232-(+)', 'bis(maltolato)oxovanadium(IV)']:
        print(id_, '=>', penc(id_))
        assert id_ == pdec(penc(id_))
    """
    # unquote_plus seems to accept characters that aren't encoded by urllib.parse.quote_plus
    return urllib.parse.unquote_plus(s)

def pfile(asset_id=None, compound_id=None, struct_id=None, screen_id=None, step='{prev_steps}', base='{base}', suffix='', v=False):
    """
        Possible alternative to pf() above "asset id"
        For a small, finite set of "entity identifiers" (e.g. model; drug, model, cavity)
        define adhoc manual prefixes (using wildcard restrictions to ease parsing)
        - swiss_model_id:
        - swiss_protein_id:
    """
    if asset_id is None:
        l_asset_id = []
        if not(compound_id is None):
            if compound_id != '{}':
                (compound_ns, compound_name) = compound_id.split('_', maxsplit=1)
                if compound_ns == 'prism' and compound_name[0] != '%':
                    compound_pref = f'{compound_ns}_{compound_name[0]}'
                elif compound_ns == 'prism' and compound_name[0] == '%':
                    compound_pref = f'{compound_ns}_'
                else:
                    compound_pref = compound_ns
            else:
                compound_pref = '{compound_pref}'
                compound_id = '{compound_id}'
            l_asset_id.append(compound_pref)
            l_asset_id.append(compound_id)

        if not(struct_id is None):
            if struct_id != '{}':
                struct_pref = os.path.join(struct_id[:2], struct_id[2:4], struct_id[4:6])
                #struct_pref = struct_id[:2]
            else:
                struct_pref = '{struct_pref}'
                struct_id = '{struct_id}'
            l_asset_id.append(struct_pref)
            l_asset_id.append(struct_id)

        if not(screen_id is None):
            l_asset_id.append(screen_id)

        asset_id = os.path.join(*l_asset_id)

    # add remaining elements: base, step, suffix
    filepath = os.path.join(base, step, f'{asset_id}{suffix}')
    if v: print('pfile: ', filepath)
    return filepath

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
    else:
        return set(map(int, x))

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
    return pd.read_csv(workpath('23.11.01_human_protein_map/structures_23.11.1.tsv'), sep='\t')

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

def read_pockets_resid():
    # Pockets, one per 
    #Q96EK7-F1      157
    cols_ = ['uniprot_id', 'struct_id', 'pocket_id', 'pocket_resid', 'pocket_score_combined_scaled']
    q_ = 'pocket_score_combined_scaled >= 900'
    return read_pockets().query(q_).explode('pocket_resid')[cols_].sort_values('pocket_score_combined_scaled', ascending=False).groupby(['uniprot_id', 'pocket_resid'])\
        .head(1)#.query('struct_id == "Q96EK7-F1" & pocket_resid == 157')

@functools.cache
def read_ppi_resid(pdockq=.5):
    df_models = read_af2_human_interactions(pdockq=pdockq)
    cols_ = ['uniprot_id', 'ifresid']
    q_ne_ = 'protein1 != protein2'
    q_eq_ = 'protein1 == protein2'
    df_interfaces = pd.concat([
        df_models.query(q_ne_).rename({'protein1': 'uniprot_id', 'residues1': 'ifresid',}, axis=1)[cols_],
        df_models.query(q_eq_).rename({'protein1': 'uniprot_id', 'residues1': 'ifresid',}, axis=1)[cols_],
        df_models.query(q_ne_).rename({'protein2': 'uniprot_id', 'residues2': 'ifresid',}, axis=1)[cols_],
    ], axis=0)
    df_interfaces['ifresid'] = df_interfaces['ifresid'].map(parse_resid)
    df_ifresid = df_interfaces.explode('ifresid').drop_duplicates(keep='first')
    return df_ifresid

def merge_missense(frame, variant_col, *args, **kwargs):
    cols_out = frame.columns.tolist()
    printlen(frame, 'raw records')
    frame[['_uniprot_id', '_aa_pos', '_aa_ref', '_aa_alt']] = frame.apply(lambda r: parse_varstr(r[variant_col]), axis=1, result_type='expand')

    variant_frame = query_missense(frame[variant_col]).set_index('variant_id')[['am_pathogenicity', 'am_class', 'pred_ddg']]
    frame = frame.merge(variant_frame, left_on=variant_col, right_index=True, *args, **kwargs)
    printlen(frame, 'records matched to predictions')

    frame['am_label'] = frame['am_class'] == 'pathogenic'
    frame['pred_ddg_label'] = frame['pred_ddg'].map(lambda pred_ddg: pred_ddg > 2)
    printlen(frame.query('pred_ddg_label'), 'annotated as destabilizing')
    cols_out += ['am_pathogenicity', 'am_class', 'am_label', 'pred_ddg', 'pred_ddg_label']

    frame = frame.merge(read_pockets_resid(), left_on=['_uniprot_id', '_aa_pos'], right_on=['uniprot_id', 'pocket_resid'], how='left')
    frame['pocket_label'] = frame['pocket_score_combined_scaled'].map(lambda pocket_score_combined_scaled: pocket_score_combined_scaled == pocket_score_combined_scaled)
    printlen(frame.query('pocket_label'), 'annotated with pockets')
    cols_out += ['pocket_label']

    frame = frame.merge(read_ppi_resid(), left_on=['_uniprot_id', '_aa_pos'], right_on=['uniprot_id', 'ifresid'], suffixes=('', '_ppi'), how='left')
    frame['interface_label'] = frame['ifresid'].map(lambda ifresid: ifresid == ifresid)
    printlen(frame.query('interface_label'), 'annotated with interfaces')
    cols_out += ['interface_label']

    return frame[cols_out]

def pvalue_label(pvalue):
    if pvalue <= .0001:
        return '****'
    elif pvalue <= .001:
        return '***'
    elif pvalue <= .01:
        return '**'
    elif pvalue <= .05:
        return '*'
    else:
        return 'ns'

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
def read_pockets(insert_uniprot_id=True):
    fp_ = workpath('23.10.16_af2_human_pockets/23.10.16_af2_human_pocket_summary.tsv.gz')
    df_ = pd.read_csv(fp_, sep='\t')
    df_['pocket_resid'] = df_['pocket_resid'].map(parse_resid)
    if insert_uniprot_id:
        df_.insert(loc=0, column='uniprot_id', value=df_['struct_id'].map(strip_fragment_id))
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
        print(self.uniprot_id_bait, 'bait uniprot id')
        self.df_interactors = df_interactors.copy()
        printlen(self.df_interactors, 'interaction models sharing the bait')
        self.nresid_bait = read_chain_len(self.df_interactors.head(1).pdb.squeeze(), self.df_interactors.head(1).bait_chain.squeeze())
        print(self.nresid_bait, 'number of bait residues')

    def build_matrix(self):
        def apply_(r):
            return [resid in parse_resid(r['bait_ifresid']) for resid in range(1, self.nresid_bait + 1)]

        self.df_matrix = pd.DataFrame(self.df_interactors.apply(apply_, axis=1).to_list(), index=self.df_interactors['interactor_id'], columns=range(1, self.nresid_bait + 1))
        print(self.df_matrix.shape, 'shape matrix dimensions')
        print(sum(self.df_matrix), 'non-zero entries in residue matrix')

    def linkage(self, criterion='distance', t=.9):
        self.linkage_ = scipy.cluster.hierarchy.linkage(self.df_matrix, method='average', metric='jaccard')
        #self.df_interactors['labels'] = scipy.cluster.hierarchy.fcluster(Z=self.linkage_, criterion='distance', t=.95)
        #self.df_interactors['labels'] = scipy.cluster.hierarchy.fcluster(Z=self.linkage_, criterion='maxclust', t=2)
        self.df_interactors['labels'] = scipy.cluster.hierarchy.fcluster(Z=self.linkage_, criterion=criterion, t=t)
        print(self.df_interactors['labels'].value_counts())

    #def plot_clustermap_default(self):
    #    sns.clustermap(self.df_matrix, metric='jaccard', row_cluster=True, col_cluster=False, cbar_pos=None, figsize=(8,4))

    def clustermap(self, fname=None):
        row_colors_ = [* self.df_interactors['labels'].map(lambda label: matplotlib.colormaps['tab10'].colors[label - 1]) ]
        if not (fname is None):
            plt.ioff()

        self.clustermap_ = sns.clustermap(self.df_matrix, row_linkage=self.linkage_, col_cluster=False, cbar_pos=None, figsize=(8,4), row_colors=row_colors_)
        plt.title(self.uniprot_id_bait)
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
        #df_interactors_top = self.df_interactors.sort_values('pdockq', ascending=False).groupby('labels').head(1) # Align to top interaction model
        df_interactors_top = self.df_interactors.sort_values(['labels', 'pdockq'], ascending=False)#.head(1) # Align to top interaction model
        cmd_ = []
        cmd_.append('delete all')
        cmd_.append('bg_color white')
        ref_ = df_interactors_top.head(1).squeeze() # Reference structure to align against
        for i, r in df_interactors_top.iterrows():
            if fname is None:
                fp_ = os.path.join('~/work-euler/', r.pdb.removeprefix('/cluster/work/beltrao/jjaenes/'))
            else:
                fp_ = r.pdb
            cmd_.append(f'load {fp_}')

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
        return df_.sort_values('interaction_id')[['interaction_id', 'Publication Source']].groupby('interaction_id').size().to_frame(name='biogrid_nrecords')
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

def read_clinvar():
    df_clinvar = pd.read_csv(workpath('23.06.02_clinvar/clinvar_mapped.tsv'), sep='\t',
        dtype={
            'CHROM': str,
            'RS': str,
            'CLNDISDBINCL': str,
            'CLNDNINCL': str,
            'CLNSIGINCL': str,
            'DBVARID': str,
            'Chromosome': str,
            'CLNSIG': str,
            'CLNVC': str,
            #'MC': str,
            'Amino_acid_position': int,
        },
        #nrows=10,
    )
    printlen(df_clinvar, 'raw records')
    df_clinvar = df_clinvar.query('CLNVC == "single_nucleotide_variant"').copy()
    printlen(df_clinvar, 'after selecting for single nucleotide variants')
    missense_ = df_clinvar.MC.str.endswith("missense_variant").fillna(False)
    df_clinvar = df_clinvar[ missense_ ]
    printlen(df_clinvar, 'after selecting for missense variants')
    return df_clinvar

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

def get_uniprot_sequence(accession):
    url = f'https://rest.uniprot.org/uniprotkb/search?query={accession}&fields=sequence'
    return requests.get(url).json()['results'][0]['sequence']['value']

__all__.extend([name for (name, thing) in locals().items() if callable(thing)]) #https://stackoverflow.com/questions/18451541/getting-a-list-of-locally-defined-functions-in-python
