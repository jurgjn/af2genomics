
import collections, io, subprocess

import numpy as np,  pandas as pd

from af2genomics.common import *

try:
    import Bio, Bio.PDB, Bio.SVDSuperimposer, Bio.SeqUtils
except ImportError:
    print('biopython not found; if needed, install with: conda install conda-forge::biopython')

""" Exponential consensus ranking, constant sigma
Ref: https://doi.org/10.1038/s41598-019-41594-3
"""
def ecr_score(rank, sigma=5):
    return math.exp(-rank / sigma) / sigma

def ecr(ranks, sigma):
    return sum(map(ecr_score, ranks))

def ecr_var(ranks, sigmas):
    return sum(map(ecr_score, ranks, sigmas))

#print(ecr([1,2,3], 10))
#print(ecr_var([1,2,3], [10,10,10]))

def bimodality_coefficient(x):
    # https://doi.org/10.1038/s43018-019-0018-6
    g = scipy.stats.skew(x.dropna())
    k = scipy.stats.kurtosis(x.dropna())
    n = len(x.dropna())
    return (g**2 + 1) / (k + (3*(n-1)**2) / ((n-2)*(n-3)))

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

def edc(pdb_file, disease_resseq, feature_resseq=None, healthy_resseq=None, restrict_chain=True, min_atoms=5, min_pLDDT=0, agg_func=np.nanmedian):
    """
    https://doi.org/10.1371/journal.pone.0307312.t001
    > EDC is calculated as previously described [5], however, alpha carbon atoms with a pLDDT < 70 are excluded from the calculation
    > and only proteins with at least 5 pathogenic variants after this procedure are used for the analysis.
    """

    def iter_chain_resseq(chain_resseq):
        for chain, resseq in chain_resseq.items():
            for resseq_i in resseq:
                yield chain, resseq_i

    disease_resseq_pairs = list(iter_chain_resseq(disease_resseq))
    feature_resseq_pairs = disease_resseq_pairs if feature_resseq is None else list(iter_chain_resseq(feature_resseq))
    if not(healthy_resseq is None):
        healthy_resseq_pairs = list(iter_chain_resseq(healthy_resseq))

    parser = Bio.PDB.PDBParser(QUIET=True)
    struct = parser.get_structure(pdb_file, pdb_file)

    def get_resseq(a): return a.get_parent().get_id()[1]
    def get_chain(a): return a.get_parent().get_parent().get_id()
    def is_CA(a): return a.get_id() == 'CA'
    def is_pLDDT(a): return a.get_bfactor() >= min_pLDDT
    def check_restrict_chain(a):
        if restrict_chain: # flag set (default), restrict only to atoms from chains with at least one disease mutation (as in original implementation)
            return get_chain(a) in disease_resseq.keys()
        else:
            return True

    def is_disease(a): 
        return (get_chain(a), get_resseq(a)) in disease_resseq_pairs
    def is_feature(a): 
        return (get_chain(a), get_resseq(a)) in feature_resseq_pairs
    def is_healthy(a):
        if not(healthy_resseq is None):
            return (get_chain(a), get_resseq(a)) in healthy_resseq_pairs
        else:
            return not is_disease(a)

    all_atoms = list(Bio.PDB.Selection.unfold_entities(entity_list=struct[0], target_level='A'))
    filtered_atoms = list(filter(lambda a: is_CA(a) and is_pLDDT(a) and check_restrict_chain(a), all_atoms))
    disease_atoms = list(filter(lambda a: is_disease(a), filtered_atoms))
    feature_atoms = list(filter(lambda a: is_feature(a), filtered_atoms))
    healthy_atoms = list(filter(lambda a: is_healthy(a), filtered_atoms))

    if not ((len(disease_atoms) >= min_atoms) and (len(feature_atoms) >= min_atoms) and (len(healthy_atoms) >= min_atoms)):
        return float('nan')

    def neq_dist(atm1, atm2): return float('nan') if atm1 == atm2 else atm1 - atm2
    def feature_mindist(atm): return np.nanmin([neq_dist(atm, feature_atm) for feature_atm in feature_atoms])

    disease_mindist = [* map(feature_mindist, disease_atoms) ]
    healthy_mindist = [* map(feature_mindist, healthy_atoms) ]
    return agg_func(np.log(healthy_mindist)) / agg_func(np.log(disease_mindist))

'''
# Expected: 1.152377
af2genomics.metrics.edc(
    pdb_file=workpath('24.10.23_goflof/Gerasimavicius2022/marshlab-mutation_clustering-5655fdd/1grnA.pdb'),
    disease_resseq={'A': [23, 66, 64, 42, 21, 81, 68, 116, 83, 159, 171, 46, 23, 34, 170,]},
    min_pLDDT=0,
    agg_func=np.mean,
)
'''
'''
# KCTD1/ENSG00000134504/Q719H9 - paper: 1.90
print(edc('/cluster/work/beltrao/jjaenes/24.11.25_atlas_of_homo_oligomerization/AF_dimer_models_full_length_relaxed/Q719H9_V1_4_relaxed.pdb', [20, 30, 31, 33, 62, 69, 74], min_pLDDT=0))
print(edc(pfile(struct_id=f'Q719H9-F1', step='af2', suffix='.pdb', base='/cluster/work/beltrao/jjaenes/24.06.10_af2genomics/results/human'), [20, 30, 31, 33, 62, 69, 74], min_pLDDT=0))
pos_ = [20, 30, 31, 33, 51, 62, 69, 74, 75, 78, 113, 140, 145, 146, 197, 208, 237]
print(edc('/cluster/work/beltrao/jjaenes/24.11.25_atlas_of_homo_oligomerization/AF_dimer_models_full_length_relaxed/Q719H9_V1_4_relaxed.pdb', pos_, min_pLDDT=0))
print(edc(pfile(struct_id=f'Q719H9-F1', step='af2', suffix='.pdb', base='/cluster/work/beltrao/jjaenes/24.06.10_af2genomics/results/human'), pos_, min_pLDDT=0))

# ADSL/ENSG00000239900/P30566 - paper: 0.89
pos_ = [242, 450, 472, 259, 332, 452, 300, 90, 225, 426, 438, 26, 141, 100, 114, 204, 318, 190, 396, 141, 204, 303, 246, 396, 426, 114]
print(edc('/cluster/work/beltrao/jjaenes/24.11.25_atlas_of_homo_oligomerization/AF_dimer_models_full_length_relaxed/P30566_V1_3_relaxed.pdb', pos_, min_pLDDT=0))
print(edc(pfile(struct_id=f'P30566-F1', step='af2', suffix='.pdb', base='/cluster/work/beltrao/jjaenes/24.06.10_af2genomics/results/human'), pos_, min_pLDDT=0))
'''