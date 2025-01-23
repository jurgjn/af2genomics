
import pandas as pd

from af2genomics.common import *

def read_struct():
    return pd.read_csv(workpath('24.06.10_af2genomics/results/human/af2.tsv'), sep='\t').query('n_frags == 1').reset_index(drop=True)
