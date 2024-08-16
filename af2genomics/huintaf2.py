
import io, subprocess, pandas as pd

from af2genomics.common import *

def pDockQ(pdb):
    py_ = af2genomicspath('software/huintaf2/bin/pDockQ.py')
    ps_ = subprocess.run([py_, '-p', pdb], capture_output=True)
    pDockQ = float(ps_.stdout)
    return pDockQ

def pDockQ2(pdb):
    py_ = af2genomicspath('software/huintaf2/bin/pDockQ2.py')
    ps_ = subprocess.run([py_, '-p', pdb, '--verbose'], capture_output=True)
    #print(f'{py_} -p {pdb} --verbose')
    #print(ps_.stdout)
    #names_ = ['Name', 'id', 'Chain', 'pDockQ',   'NumRes',   'IF_plDDT',    'plDDTNumDiso1+90',        'NumDiso1-70-90', 'NumDiso1-50-70', 'NumDiso1-50', 'Length'] # script header
    #names_ = ['Name', 'i',  'id',    'pdockq[i]','NumRes[i]','IF_plDDT[i]', 'plDDT[i]','NumDiso[i][0]','NumDiso[i][1]',  'NumDiso[i][2]',  'NumDiso[i][3]', 'NumResOverlap[i]', 'length[i]'] # script variables
    names_ =  ['Name', 'id', 'Chain', 'pDockQ',   'NumRes',   'IF_plDDT',    'plDDT',   'NumDiso1+90',  'NumDiso1-70-90', 'NumDiso1-50-70', 'NumDiso1-50',   'NumResOverlap', 'Length'] # merge
    df_ = pd.read_csv(io.BytesIO(ps_.stdout), sep=',', names=names_, skiprows=1)
    return df_