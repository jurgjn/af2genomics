#!/usr/bin/env python

""" 
Alphafold output the predicted quality of the amino acid (plddt) as a b_factor (0->100). 
In general, below 50 means the conformation is unreliable
Here we remove amino acids where the rolling average over a given window size is less than a threshold.
For the fragment to be kept, it needs to be a minimum length ie not lots of single amino acid fragments
Version 1.0
David F Burke
March 2022
"""

import os.path
import argparse
from Bio import PDB
from Bio.PDB import PDBIO
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--pdbfile", dest="pdbfile", default="model.pdb", help="PDB file (model.pdb)")
parser.add_argument("--outpdb", dest="outpdb", default="model_trimmed.pdb", help="Output PDB file(model_trimmed.pdb)")
parser.add_argument("--bfactor", dest="bfactor", type=int, default=50, help="Minimum Bfactor to keep(50)")
parser.add_argument("--window", dest="window", type=int, default=10, help="Window size(10)")
parser.add_argument("--loop", dest="looplength", type=int, default=10, help="Minimum length of fragment to be removed(10)")
parser.add_argument("--length", dest="fraglength", type=int, default=10, help="Minimum length of fragment to keep(10)")
options = parser.parse_args()
pdb1 =options.pdbfile

io=PDBIO()
pdbparser = PDB.PDBParser()

# read in bfactors ( plddt scores) for each residue
structure = pdbparser.get_structure("model", pdb1)

model = structure[0]
# Find list of ids to remove
for chain in structure.get_chains():
   res=[]
   residue_ids_to_remove=[]
   residue_ids_to_keep=[]
   bfactors=[]
   for residue in chain.get_residues():
       for atom in residue.get_unpacked_list():
          if atom.id == 'CA':
             bfactors.append(atom.bfactor)
             res.append(residue.id[1])

   # Make rolling average
   df = pd.DataFrame(bfactors)
   df.columns=['bfactors']
   bf_average=df.rolling(window=options.window, min_periods=1).mean()
   for residue in chain.get_residues():
       if bf_average['bfactors'][int(residue.id[1]-1)]<options.bfactor:
           print("Chain",chain,"Removing ",residue.full_id, bf_average['bfactors'][int(residue.id[1]-1)])
           residue_ids_to_remove.append(residue.id[1])
       else:
           print("Chain ",chain,"Keeping ",residue.full_id, bf_average['bfactors'][int(residue.id[1]-1)])
           residue_ids_to_keep.append(residue.id[1])

   # If fragment selected to be removed is too short, then choose to keep it. ie dont remove short loops
   fraglength=1
   for resindex in range(1,len(residue_ids_to_remove)):
       if residue_ids_to_remove[resindex]-residue_ids_to_remove[resindex-1]==1:
          fraglength=fraglength+1
       else:
          if fraglength<options.looplength:
            for res in range(resindex-fraglength,resindex):
              id=residue_ids_to_remove[res]
              print("Chain ",chain,"Loop too short - Keeping ",id,res)
              residue_ids_to_remove[res]=-1
              residue_ids_to_keep.append(id)
            fraglength=1

   residue_ids_to_keep.sort()

   # If fragment to be kept is too short, the remove it. 
   fraglength=1
   for resindex in range(1,len(residue_ids_to_keep)):
        if residue_ids_to_keep[resindex]-residue_ids_to_keep[resindex-1]==1:
            fraglength=fraglength+1
        else:
            if fraglength<options.fraglength:
               for res in range(resindex-fraglength,resindex):
                  residue_ids_to_remove.append(residue_ids_to_keep[res])
                  print("Chain ",chain,"Fragment too short - Remove ",residue_ids_to_keep[res])
            fraglength=1

   # now get rid of short loops (set to -1) from the remove list 
   residue_ids_to_remove.sort()
   residue_ids_to_remove=[i for i in residue_ids_to_remove if i>-1]
   print("Chain ",chain,"Removing",residue_ids_to_remove)

   #now remove the residues
   for id in residue_ids_to_remove:
        chain.detach_child((' ', id, ' '))

# save structure
io.set_structure(model)
io.save(options.outpdb)
