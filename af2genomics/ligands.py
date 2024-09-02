
import pandas as pd, openbabel, openbabel.pybel

def read_gnina(fp):
    def get_(mol_, key_):
        return float(mol_.data[key_]) if key_ in mol_.data else float('NaN')
    cols_ = ['minimizedAffinity', 'CNNscore', 'CNNaffinity', 'CNN_VS']
    df_ = pd.DataFrame.from_records([[
        get_(mol, 'minimizedAffinity'),
        get_(mol, 'CNNscore'),
        get_(mol, 'CNNaffinity'),
        get_(mol, 'CNN_VS'),
    ] for mol in openbabel.pybel.readfile('sdf', fp) ], columns=cols_)
    return df_
