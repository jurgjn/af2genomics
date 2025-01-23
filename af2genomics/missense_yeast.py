
import pandas as pd

from af2genomics.common import *

def read_struct():
    return pd.read_csv(workpath('24.06.10_af2genomics/results/yeast/af2.trim_bf.tsv'), sep='\t')
