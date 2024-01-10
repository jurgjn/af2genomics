
import ast, collections, datetime, functools, inspect, itertools, math, os, pandas as pd, requests, sqlite3

import scipy as sp, scipy.stats
import Bio, Bio.PDB, Bio.SVDSuperimposer
import prody

__all__ = ['RANDOM_SEED']
# Fix `RANDOM_SEED` for (partial) reproducibility
RANDOM_SEED = 4 # https://xkcd.com/221

__all__ = ['workpath']
def workpath(path):
    #dir_ = os.path.dirname(__file__)
    dir_ = '/cluster/work/beltrao/jjaenes'
    return os.path.join(dir_, path)

__all__.append('uf')
def uf(x):
    return '{:,}'.format(x)

__all__.append('ul')
def ul(x):
    return uf(len(x))

__all__.append('printsrc')
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

__all__.append('printlen')
def printlen(x, *args, **kwargs):
    # printlen(af2_uniprot_id(), 'human AF2 single-fragment structures')
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

__all__.append('phead')
def phead(df, n=3):
    print(uf(len(df)), 'records')
    return df.head(n).transpose()

__all__.append('pmerge')
def pmerge(left, right, *args, **kwargs):
    left_merge_right = left.merge(right, *args, **kwargs)
    print('merge:', uf(len(left_merge_right)), 'left:', uf(len(left)), 'right:', uf(len(right)))
    return left_merge_right

__all__.append('read_af2_pos')
def read_af2_pos(add_residue_id=False, **kwargs):
    kwargs.setdefault('filepath_or_buffer', workpath('23.10.05_VEP_scores/23.10.18_af2_pos.tsv'))
    kwargs.setdefault('sep', '\t')
    df_ = pd.read_csv(**kwargs)
    if add_residue_id:
        df_.insert(loc=1, column='residue_id', value=[ f'{r.uniprot_id}/{r.pos}' for i, r in df_.iterrows() ])
    return df_

__all__.append('read_human_missense_scores')
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

__all__.append('read_interface_full')
def read_interface_full():
    fp_ = workpath('23.09.07_dburke_af2interactions_mutfunc/interface_full.list')
    df_ppi = pd.read_csv(os.popen(f'grep -v "Error" {fp_}'), sep=' ')
    df_ppi['pair_nmodels'] = df_ppi['pair'].str.count('_') + 1
    df_ppi = df_ppi.query('pair_nmodels == 2').copy()
    df_ppi[['pairA', 'pairB']] = df_ppi['pair'].str.split('_', expand=True)
    #print(df_ppi.info())
    return df_ppi

__all__.append('read_interface_full_strict')
def read_interface_full_strict():
	fp_ = workpath('23.09.07_dburke_af2interactions_mutfunc/interface_full_strict.list')
	df_ppi = pd.read_csv(os.popen(f'grep -v "Error" {fp_}'), sep=' ')
	df_ppi['pair_nmodels'] = df_ppi['pair'].str.count('_') + 1
	df_ppi = df_ppi.query('pair_nmodels == 2').copy()
	df_ppi[['pairA', 'pairB']] = df_ppi['pair'].str.split('_', expand=True)
	return df_ppi

__all__.append('parse_resid')
def parse_resid(s):
    """ Return various list-of-residues representations as set-of-ints, e.g.:
    parse_resid('')
    parse_resid('1')
    parse_resid('1,2,3,4')
    parse_resid('{1,2,3,4}')
    parse_resid('[1,2,3,4]')
    """
    if s == '' or s != s:
        return set()
    x = ast.literal_eval(s)
    if isinstance(x, int):
        return set([x])
    else:
        return set(map(int, x))

__all__.append('fillna_set')
def fillna_set(s):
    # Fill NA-s with an empty set
    # https://stackoverflow.com/questions/33199193/how-to-fill-dataframe-nan-values-with-empty-list-in-pandas
    return s.apply(lambda d: d if isinstance(d, set) else set())

__all__.append('g_convert')
def g_convert(query=["CASQ2", "CASQ1", "GSTO1", "DMD", "GSTM2"], target='UNIPROTSWISSPROT_ACC', organism='hsapiens', numeric_namespace='ENTREZGENE_ACC'):
    """
    Query for HGNC gene names using g:convert (https://biit.cs.ut.ee/gprofiler/convert)
    """
    r = requests.post(url='https://biit.cs.ut.ee/gprofiler/api/convert/convert/', json=locals())
    df_ = pd.DataFrame(r.json()['result'])
    return df_

__all__.append('pquery')
def pquery(df_, *args, **kwargs):
    df_q_ = df_.query(*args, **kwargs)
    printsrc('pquery:', uf(len(df_q_)), 'of', uf(len(df_)), 'rows selected with', args[0])
    return df_q_

__all__.append('calver')
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

__all__.append('read_af2_pLDDT')
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

__all__.append('read_swiss_resid')
def read_swiss_resid():
    df_swiss = pd.read_csv(workpath('23.10.31_SWISS_human_2023-09-14/SWISS-MODEL_Repository/INDEX'), comment='#', sep='\t')
    df_swiss['pos'] = [* map(lambda from_, to_: set(range(from_, to_ + 1)), df_swiss['from'], df_swiss['to']) ]
    df_resid = df_swiss.groupby(['UniProtKB_ac', 'provider']).agg(pos=('pos', lambda x: set.union(*(x_i for x_i in x)))).unstack()
    df_resid.columns = df_resid.columns.droplevel()
    df_resid.rename({'PDB': 'resid_pdb', 'SWISSMODEL': 'resid_swiss'}, axis=1)
    df_resid = df_resid.reset_index().rename({'UniProtKB_ac': 'uniprot_id', 'PDB': 'resid_pdb', 'SWISSMODEL': 'resid_swiss'}, axis=1)
    return df_resid

__all__.append('nan_to_set')
def nan_to_set(x):
    return x if x == x else set()

__all__.append('read_Cheng2023_s5')
def read_Cheng2023_s5():
    # Supplementary Data S5: List of variants and AlphaMissense predictions for the ClinVar benchmark.
    df_ = pd.read_csv(workpath('23.10.05_VEP_scores/science.adg7492_data_s1_to_s9/science.adg7492_data_s5.csv'), sep=',', dtype={'label': int})
    df_['protein_variant'] = df_['protein_variant'].str.replace(':','/')
    return df_

__all__.append('read_Cheng2023_s6')
def read_Cheng2023_s6():
    # Supplementary Data S6: List of variants and AlphaMissense predictions for the Cancer hotspot mutations benchmark.
    df_ = pd.read_csv(workpath('23.10.05_VEP_scores/science.adg7492_data_s1_to_s9/science.adg7492_data_s6.csv'), sep=',')
    df_['protein_variant'] = df_['protein_variant'].str.replace(':','/')
    return df_

__all__.append('read_Cheng2023_s7')
def read_Cheng2023_s7():
    # Supplementary Data S7: List of variants and AlphaMissense predictions for the Deciphering Developmental Disorders benchmark.
    df_ = pd.read_csv(workpath('23.10.05_VEP_scores/science.adg7492_data_s1_to_s9/science.adg7492_data_s7.csv'), sep=',', dtype={'label': int})
    df_['protein_variant'] = df_['protein_variant'].str.replace(':','/')
    return df_

__all__.append('read_Cheng2023_s8')
def read_Cheng2023_s8():
    # Supplementary Data S8: List of variants and AlphaMissense predictions for the ProteinGym benchmark.
    df_ = pd.read_csv(workpath('23.10.05_VEP_scores/science.adg7492_data_s1_to_s9/Supplementary_Data_S8_proteingym.csv'), sep=',')#, dtype={'label': int})
    #df_['protein_variant'] = df_['protein_variant'].str.replace(':','/')
    return df_

__all__.append('read_structures')
def read_structures():
    return pd.read_csv(workpath('23.11.01_human_protein_map/structures_23.11.1.tsv'), sep='\t')

__all__.append('read_missense')
def read_missense():
    with sqlite3.connect(workpath('23.11.01_human_protein_map/missense_23.11.1.sqlite')) as db:
        df_ = pd.read_sql_query(sql='SELECT * FROM missense', con=db)
    return df_

__all__.append('query_missense')
def query_missense(variants):
    #https://stackoverflow.com/questions/28735213/pandas-read-sql-with-a-list-of-values-for-where-condition
    with sqlite3.connect(workpath('23.11.01_human_protein_map/missense_23.11.1.sqlite')) as db:
        df_ = pd.read_sql_query(sql='SELECT * FROM missense WHERE variant_id in ' + str(tuple(variants)), con=db)
    return df_

__all__.append('query_missense_uniprot_id')
def query_missense_uniprot_id(uniprot_id):
    with sqlite3.connect(workpath('23.11.01_human_protein_map/missense_23.11.1.sqlite')) as db:
        df_ = pd.read_sql_query(sql=f'SELECT * FROM missense WHERE variant_id GLOB "{uniprot_id}*"', con=db)
    return df_

__all__.append('jaccard')
def jaccard(a, b):
    return len(a & b) / len(a | b)

__all__.append('interaction_id')
def interaction_id(uniprot_id_1, uniprot_id_2):
    # sorted pair of uniprot_id-s to compare interacting pairs from different sources
    return '_'.join(sorted([uniprot_id_1, uniprot_id_2]))

__all__.append('af2_uniprot_id')
@functools.cache
def af2_uniprot_id():
    # AF2 human single-fragment structures (set used for pocket detection, interface modelling, etc)
    return set(read_structures()['uniprot_id'])

def flatten(l):
    return [item for sublist in l for item in sublist]

__all__.append('read_summary_source')
def read_summary_source(summary=False):
    df_ = pd.read_csv(workpath('23.10.02_dburke_kcl_OneDrive/summary_source.out.bz2'), na_values={'pdockq_fd': 'chain'}, nrows=None)\
        .rename({'#protdir': 'protdir'}, axis=1)
    print(uf(len(df_)), 'raw records')
    def parse_(s_):
        l_id_ = s_.split('/')[1].split('_')
        if len(l_id_) == 2:
            return interaction_id(*l_id_)
        else:
            return '.'
    df_.insert(0, 'interaction_id', [* df_['protdir'].map(parse_) ])
    df_.insert(2, 'folding_method', df_['protdir'].map(lambda s: s.split('/')[-1]))
    df_ = df_.query('interaction_id != "."').copy()
    print(uf(len(df_)), 'after discarding non-dimers')
    df_[['uniprot_id_A', 'uniprot_id_B']] = df_['interaction_id'].str.split('_', expand=True)
    df_ = df_.query('(uniprot_id_A in @af2_uniprot_id()) & (uniprot_id_B in @af2_uniprot_id())').copy()
    print(uf(len(df_)), 'after keeping uniprot_id-s in AF2 single fragment structures')
    print(uf(len(set(df_['interaction_id']))), 'unique interaction_id-s')
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

__all__.append('read_ppi_reselect')
def read_ppi_reselect():
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
    printlen(df_, 'raw records')

    def parse_(s_):
        l_id_ = s_.split('_')
        if len(l_id_) == 2:
            return interaction_id(*l_id_)
        else:
            return '.'
    df_.insert(0, 'interaction_id', [* df_['pair'].map(parse_) ])
    df_ = df_.query('interaction_id != "."').copy()
    printlen(df_, 'after discarding non-dimers')
    df_[['uniprot_id_A', 'uniprot_id_B']] = df_['interaction_id'].str.split('_', expand=True)
    df_ = df_.query('(uniprot_id_A in @af2_uniprot_id()) & (uniprot_id_B in @af2_uniprot_id())').copy()
    printlen(df_, 'after keeping uniprot_id-s in AF2 single fragment structures')
    df_ = df_.drop_duplicates(subset=['interaction_id'], keep='first')
    printlen(df_, 'after naive de-duplication of interaction_id')
    printlen(df_.query('pdockq > .23'), 'with pdockq > .23')
    printlen(df_.query('pdockq > .5'), 'with pdockq > .5')
    df_['pdb'] = df_.apply(lambda r: workpath(f'23.12.06_ppi_reselect/af2-models-split/{r.pair.split("_")[0]}/{r.pair}.pdb'), axis=1)
    return df_.reset_index(drop=True)

__all__.append('read_pockets')
@functools.cache
def read_pockets():
    fp_ = workpath('23.10.16_af2_human_pockets/23.10.16_af2_human_pocket_summary.tsv.gz')
    df_ = pd.read_csv(fp_, sep='\t')
    df_['pocket_resid'] = df_['pocket_resid'].map(parse_resid)
    return df_

__all__.append('pLDDT_palette')
# pLDDT visualisation to match AFDB, e.g. https://alphafold.ebi.ac.uk/entry/P00533
pLDDT_palette = {
    'Very low (pLDDT < 50)': '#ff7d45ff',
    'Low (70 > pLDDT > 50)': '#ffdb13ff',
    'High (90 > pLDDT > 70)': '#64cbf3ff',
    'Very high (pLDDT > 90)': '#0053d6ff',
}
__all__.append('pLDDT_bins')
pLDDT_bins = [0, 50, 70, 90, 100]

__all__.append('to_pymol_init')
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

__all__.append('to_pymol_struct')
def to_pymol_struct(struct_id, base='~/euler-home/project/22.12_pocketomes/results/23.04_bindfunc'):
    print('load', pfile(struct_id=struct_id, step='af2', suffix='.pdb', base=base))
    print('color grey50,', struct_id)

__all__.append('to_pymol_pocket')
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
__all__.append('bait_rmsd')
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

__all__.append('read_chain_len')
def read_chain_len(fp_, chain_):
    # number of residues in a chain
    structure_ = Bio.PDB.PDBParser(QUIET=True).get_structure(fp_, fp_)
    return len([residue['CA'].coord for residue in structure_[0][chain_]]) #https://biopython.org/DIST/docs/tutorial/Tutorial.html#sec203

__all__.append('parse_varstr')
def parse_varstr(s):
    # df_var[['uniprot_id', 'aa_pos', 'aa_ref', 'aa_alt']] = df_var.apply(lambda r: parse_varstr(r['protein_variant']), axis=1, result_type='expand')
    uniprot_id, variant_id = s.split('/')
    aa_pos = variant_id[1:-1]
    aa_ref = variant_id[0]
    aa_alt = variant_id[-1]
    #print(uniprot_id, aa_pos, aa_ref, aa_alt)
    return uniprot_id, aa_pos, aa_ref, aa_alt

__all__.append('VariantEnrichment')
class VariantEnrichment:
    # struct_data, struct_uniprot_id, struct_resid
    # var_data, var_uniprot_id, var_resid
    def __init__(self, struct_data, struct_uniprot_id, struct_resid, struct_bin, var_data, var_uniprot_id, var_resid):
        self.struct_data = struct_data[[struct_uniprot_id, struct_resid, struct_bin]].copy().rename({struct_uniprot_id: 'struct_uniprot_id', struct_resid: 'struct_resid', struct_bin: 'struct_bin'}, axis=1)
        self.struct_data['struct_resid'] = self.struct_data['struct_resid'].map(parse_resid)

        # Attach set of all residues
        df_af2 = pd.read_csv('results/23.04_bindfunc/af2.tsv', sep='\t').query('n_frags == 1')[['uniprot_id', 'n_resid']]
        df_af2['resid'] = df_af2['n_resid'].map(lambda n_resid: set(range(1, n_resid + 1)))
        self.struct_data = self.struct_data.merge(df_af2[['uniprot_id', 'resid']], left_on='struct_uniprot_id', right_on='uniprot_id').drop(['uniprot_id'], axis=1)

        self.var_data = var_data[[var_uniprot_id, var_resid]].copy().rename({var_uniprot_id: 'var_uniprot_id', var_resid: 'var_resid'}, axis=1)
        #self.var_data['var_resid'] = self.var_data['var_resid'].map(parse_resid)

        self.merge_data = self.struct_data.merge(self.var_data, left_on='struct_uniprot_id', right_on='var_uniprot_id', how='left')
        def fillna_(x): return x if x == x else set()
        self.merge_data['var_uniprot_id'] = self.merge_data['var_uniprot_id'].map(fillna_)
        self.merge_data['var_resid'] = self.merge_data['var_resid'].map(fillna_)

        self.merge_data['n_struct_var'] = self.merge_data.apply(lambda r: len(r['struct_resid'] & r['var_resid']), axis=1)
        self.merge_data['n_struct_novar'] = self.merge_data.apply(lambda r: len(r['struct_resid'] - r['var_resid']), axis=1)
        self.merge_data['n_nostruct_var'] = self.merge_data.apply(lambda r: len((r['resid'] - r['struct_resid']) & r['var_resid']), axis=1)
        self.merge_data['n_nostruct_novar'] = self.merge_data.apply(lambda r: len((r['resid'] - r['struct_resid']) - r['var_resid']), axis=1)

    def fisher_exact_(k, l, m, n):
        return sp.stats.fisher_exact([[k, l], [m, n]])

    def fisher_exact(self):
        self.merge_data_sums_ = self.merge_data[['n_struct_var', 'n_struct_novar', 'n_nostruct_var', 'n_nostruct_novar']].sum(axis=0)
        return VariantEnrichment.fisher_exact_(
            k=self.merge_data_sums_.n_struct_var,
            l=self.merge_data_sums_.n_struct_novar,
            m=self.merge_data_sums_.n_nostruct_var,
            n=self.merge_data_sums_.n_nostruct_novar,
        )

    def fisher_exact_groupby(self):
        sums_ = self.merge_data.groupby('struct_bin')[['n_struct_var', 'n_struct_novar', 'n_nostruct_var', 'n_nostruct_novar']].sum()
        sums_[['statistic', 'pvalue']] = sums_.apply(lambda r: VariantEnrichment.fisher_exact_(r.n_struct_var, r.n_struct_novar, r.n_nostruct_var, r.n_nostruct_novar), axis=1, result_type='expand')
        return sums_