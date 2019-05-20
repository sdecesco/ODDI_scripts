#!/usr/bin/env python
from __future__ import print_function
# ----------------------------------------------------------------------#
# 			                   CovFinder v0.1	                		#
# --------------------------- Compatibility ----------------------------#
# 					Compatible with Python 2.7 & 3.5					#
# ----------------------Prerequiste installations-----------------------#
#																		#
# 		                                                            	#
# ----------------------------------------------------------------------#
# (C) Dr. De Cesco Stephane - v3.1 - 23/01/2017							#
# ----------------------------------------------------------------------#
# -------------------------- MODULES IMPORT ----------------------------#

import os
import os.path
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors as Desc
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols as fps
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
import numpy as np
import argparse
import operator
import progressbar

widgets = [' [', progressbar.Timer(), '] ', progressbar.Bar(), ' [', progressbar.Percentage(), '] ', ' (',
           progressbar.ETA(), ') ']

#library file that you wish to compore other library to
original = '/ssddata/sdecesco/data/Scripts/CovFinder/CovFinder/data/libraries/Ubiquigen_DUBtarget.sdf'

#List of other libraries
list_to_compare = [
    '/ssddata/sdecesco/data/Scripts/CovFinder/CovFinder/data/libraries/2016-11 Asinex Soft Electrophiles - 8209.sdf',
    '/ssddata/sdecesco/data/Scripts/CovFinder/CovFinder/data/libraries/Chem Div covalent_inhibitor library_5637.sdf',
    '/ssddata/sdecesco/data/Scripts/CovFinder/CovFinder/data/libraries/Enamine_Covalent_Fragments.sdf',
    '/ssddata/sdecesco/data/Scripts/CovFinder/CovFinder/data/libraries/Enamine_Cys_focused_Covalent_Fragments.sdf',
    '/ssddata/sdecesco/data/Scripts/CovFinder/CovFinder/data/libraries/Life Chemicals_Covalent_Fragment_Library.sdf',
    '/ssddata/sdecesco/data/Scripts/CovFinder/CovFinder/data/libraries/UORSY_covalent_fragments.sdf']

original_mol = Chem.SDMolSupplier(original)

compare = []
for i in list_to_compare:
    compare_mol = Chem.SDMolSupplier(i)
    library_name = i.split('/')[-1]
    counter_bar = 0
    bar = progressbar.ProgressBar(maxval=len(compare_mol), widgets=widgets)
    bar.start()
    for j in compare_mol:
        bar.update(counter_bar)
        counter_bar+=1
        compare.append((library_name, j, fps.FingerprintMol(j)))
    bar.finish()


bar = progressbar.ProgressBar(maxval=len(original_mol), widgets=widgets)
bar.start()
counter_bar = 0
with open('/ssddata/sdecesco/data/Scripts/CovFinder/CovFinder/data/libraries/ubiquigen_similar.csv','w') as output:
    out = 'Original_smiles_ubiquigen,Similar_structure_smiles,library,tanimoto\n'
    output.write(out)
    for m in original_mol:
        bar.update(counter_bar)
        counter_bar += 1
        counter = 0
        out_temp = str(Chem.MolToSmiles(m)) + ','
        mw_original = Desc.MolWt(m)
        closest = 0
        for mol_comp in compare:
            mw_compare = Desc.MolWt(mol_comp[1])
            if abs(mw_compare-mw_original) >= 50:
                continue
            sim = DataStructs.FingerprintSimilarity(mol_comp[2], fps.FingerprintMol(m))
            if sim > closest:
                closest = sim
                out_temp = str(Chem.MolToSmiles(mol_comp[1])) + ',' + str(
                    mol_comp[0]) + ',' + str(round(sim, 2)) + ','
            elif sim == 1:
                out_temp = str(Chem.MolToSmiles(mol_comp[1])) + ',' + str(
                    mol_comp[0]) + ',' + str(round(sim, 2))
                break
        to_write = str(Chem.MolToSmiles(m)) + ','+out_temp
        print(to_write)
        output.write(to_write+'\n')
    bar.finish()
