#!/usr/bin/env python
from __future__ import print_function
# ----------------------------------------------------------------------#
# 			                   PairFinder v0.1	                		#
# --------------------------- Compatibility ----------------------------#
# 					Compatible with Python 2.7 & 3.5					#
# ----------------------Prerequiste installations-----------------------#
#																		#
# 		                                                            	#
# ----------------------------------------------------------------------#
# (C) Dr. De Cesco Stephane - v0.1 - 05/05/2017							#
# ----------------------------------------------------------------------#
# -------------------------- MODULES IMPORT ----------------------------#

import argparse
import progressbar
import os
import os.path
import subprocess
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors as Desc
from rdkit.Chem import SaltRemover
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols as fps
import numpy as np
import operator
import time as t

# To remove danger of using input in python2 environments.
try:
    input = raw_input
except NameError:
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help='Name of the sdf input file [MANDATORY]', metavar='')
    parser.add_argument('-u', '--input_uncapped', help='Name of the uncapped sdf input file [MANDATORY]', metavar='')
    parser.add_argument('-o', '--output',
                        help='If you want an SDF with the annotated molecule please use this option followed by the name of the output file',
                        metavar='')
    args = parser.parse_args()
    widgets = [' [', progressbar.Timer(), '] ', progressbar.Bar(), ' [', progressbar.Percentage(), '] ', ' (',
               progressbar.ETA(), ') ']

    if args.input_file and args.input_uncapped:
        if os.path.exists(args.input_file) and os.path.exists(args.input_uncapped):
            mol_template = Chem.SDMolSupplier(args.input_file)
            mol_uncapped = Chem.SDMolSupplier(args.input_uncapped)
        else:
            print('ERROR : file ' + str(args.input_file) + ' does not exist')
            sys.exit()

    mol_uncapped = Chem.SDMolSupplier("/ssddata/sdecesco/data/Scripts/PairFinder/data/amines_uncapped.sdf")
    mol_template = Chem.SDMolSupplier("/ssddata/sdecesco/data/Scripts/PairFinder/data/amines_capped.sdf")
    # mol_template = []
    # mol_uncapped = []

    if not mol_uncapped and not mol_template:
        with open('/ssddata/sdecesco/data/Scripts/PairFinder/data/emolecules_all.txt', 'r') as smiles:
            counter = 0
            template_capped = Chem.MolFromSmarts("[H][#7](-[#6])-[#6](=O)C([H])([H])[H]")
            template_capped_2 = Chem.MolFromSmarts("[H][#7](-[#6])-[#6](=O)C([H])([H])[CH3]")
            template_amine = Chem.MolFromSmarts("N([H])[H]")
            remover = SaltRemover.SaltRemover()

            lines = smiles.readlines()

            bar = progressbar.ProgressBar(maxval=len(lines), widgets=widgets)
            bar.start()
            for line in lines:
                if counter == 0:
                    counter += 1
                    continue
                line = line.replace('\n', '')
                data = line.split(' ')
                mol = Chem.MolFromSmiles(data[0])
                if mol is None:
                    counter += 1
                    continue
                if 500 < Desc.MolWt(mol) < 100:
                    counter += 1
                    continue
                mol = Chem.AddHs(mol)
                mol = remover.StripMol(mol, dontRemoveEverything=True)
                if mol is None:
                    continue
                if mol.HasSubstructMatch(template_amine):
                    mol = Chem.RemoveHs(mol)
                    mol.SetProp('_Name', data[1])
                    mol_uncapped.append(mol)
                elif mol.HasSubstructMatch(template_capped_2) or mol.HasSubstructMatch(template_capped):
                    mol = Chem.RemoveHs(mol)
                    mol.SetProp('_Name', data[1])
                    mol_template.append(mol)
                    # print('Capped: ',data[0])

                counter += 1
                bar.update(counter)
                # print('Amine: ',data[0])
                # if counter == 20:
                #     sys.exit()
                # """For thest purpose"""
                # if len(mol_template)>500 and len(mol_uncapped) > 1000:
                #     break
        bar.finish()

        # mol_uncapped_list = []
        # for m in mol_uncapped:
        #     if m is None:
        #         continue
        #     mol_uncapped_list.append((m,fps.FingerprintMol(m)))
        # mol_uncapped = mol_uncapped_list
        out_1 = Chem.SDWriter('/ssddata/sdecesco/data/Scripts/PairFinder/data/amines_capped.sdf')
        out_2 = Chem.SDWriter('/ssddata/sdecesco/data/Scripts/PairFinder/data/amines_uncapped.sdf')
        for m in mol_template:
            out_1.write(m)
        out_1.close()
        for m in mol_uncapped:
            out_2.write(m[0])
        out_2.close()
    print(len(mol_template))
    print(len(mol_uncapped))

    mol_to_output = []
    bar = progressbar.ProgressBar(maxval=len(mol_uncapped), widgets=widgets)
    bar.start()
    counter = 0
    mol_amines = {'0-50': [], '50-100': [], '100-120': [], '120-140': [], '140-150': [], '150-160': [], '160-170': [],
                  '170-180': [],
                  '180-190': [], '190-200': [], '200-210': [], '210-220': [], '220-230': [], '230-240': [],
                  '240-250': [],
                  '250-260': [], '260-270': [], '270-280': [], '280-290': [], '290-300': [], '300-310': [],
                  '310-320': [],
                  '320-330': [], '330-340': [], '340-350': [], '350-360': [], '360-370': [], '370-380': [],
                  '380-390': [],
                  '390-400': [], '400-410': [], '410-420': [], '420-430': [], '430-440': [], '440-450': [],
                  '450-460': [],
                  '460-470': [], '470-480': [], '480-490': [], '490-500': [], 'other': []}

    for m in mol_uncapped:
        counter += 1
        bar.update(counter)
        if m is None:
            continue
        mw = Desc.MolWt(m)
        cat = 'other'
        for keys in mol_amines.keys():
            if keys == 'other':
                continue
            low_up = str(keys).split('-')
            if float(low_up[1]) >= mw > float(low_up[0]):
                cat = keys
                break
        mol_amines[cat].append((m, fps.FingerprintMol(m)))
    bar.finish()

    for keys, values in mol_amines.items():
        print(keys, '= ', len(values))

    bar = progressbar.ProgressBar(maxval=len(mol_template), widgets=widgets)
    bar.start()
    counter = 0
    for capped in mol_template:
        counter += 1
        rxn = AllChem.ReactionFromSmarts("[#6:1]-[NH:2]-[#6](-[CH3])=[O]>>[#6:1]-[NH2:2]")
        ps = rxn.RunReactants([capped])
        try:
            temp = Chem.MolToSmiles(ps[0][0])
        except IndexError:
            rxn = AllChem.ReactionFromSmarts("[#6:1]-[NH:2]-[#6](-[C][CH3])=[O]>>[#6:1]-[NH2:2]")
            ps = rxn.RunReactants([capped])
            try:
                temp = Chem.MolToSmiles(ps[0][0])
            except IndexError:
                print(Chem.MolToSmiles(capped))
                continue
        temp_mol = Chem.MolFromSmiles(temp)
        if temp_mol is None:
            continue
        fp_template = fps.FingerprintMol(temp_mol)
        mw_capped = Desc.MolWt(capped)
        cat = 'other'
        for keys in mol_amines.keys():
            if keys == 'other':
                continue
            low_up = str(keys).split('-')
            if float(low_up[1]) >= mw_capped > float(low_up[0]):
                cat = keys
                break
        # print(Chem.MolToSmiles(capped))

        for m in mol_amines[cat]:
            if DataStructs.FingerprintSimilarity(m[1], fp_template) == 1:
                # print("Match number ", match_number, ': ', Chem.MolToSmiles(mol[0]))
                capped.SetProp('Matched_amine_smiles', Chem.MolToSmiles(m[0]))
                capped.SetProp('Matched_amine_emolid', m[0].GetProp('_Name'))
                mol_to_output.append(capped)
                break
        bar.update(counter)
    bar.finish()
    out = Chem.SDWriter('/ssddata/sdecesco/data/Scripts/PairFinder/data/output.sdf')
    for m in mol_to_output:
        out.write(m)
    out.close()
