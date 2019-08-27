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
from rdkit.Chem import Draw
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
import numpy as np
import argparse
import operator
import progressbar
import pandas as pd
from sklearn.decomposition import PCA

# To remove danger of using input in python2 environments.
try:
    input = raw_input
except NameError:
    pass

widgets = [' [', progressbar.Timer(), '] ', progressbar.Bar(), ' [', progressbar.Percentage(), '] ', ' (',
           progressbar.ETA(), ') ']


def dict_of_cov_group():
    cov_smarts = {}
    cov_stats = {}
    cov_group_list = [
        '1_2_dicarbonyl_oxalyl<>[$([CD3&!R][#6,#8,#7,F,Cl,I,Br]),$([CD3&!R][#1]),$([CD2&!R])](=[OD1&!R])[$([CD3&!R][#6,#8,#7,F,Cl,Br,I]),$([CD3&!R][#1]),$([CD2&!R])]=[OD1&!R]',
        '4_subst_n_alkyltetrahydropyridines<>[a]C1=CCN([A])CC1', '4_vinyl_pyridine<>c1cnccc1[CX3&!R]=[CX3&!R]',
        'acetylene_alkyne<>[#1,#6]C#CC',
        'acrylamide<>[$([NX3H2&!R]),$([NX3H1&!R][c,C]),$([NX3&!R]([c,C])[c,C])][CX3&!R](=[OX1])[CX3&!R]=[$([CX3H2&!R]),$([CX3H1&!R][c,C]),$([CX3H0&!R]([c,C])[c,C])]',
        'acyclic_acetal<>[#6][OX2&!R][$([CX4&!RH2]),$([CX4&!RH1][#6]),$([CX4&!R]([#6])[#6])][OX2&!R][#6]',
        'acyclic_acid_halide_acyl_halide<>*[CD3&!R]([F,Cl,I,Br])=[OX1]', 'acyclic_alkyne<>[#6][CX2&!R]#[CX2&!R][#6]',
        'acyclic_ketone<>[#6&!R][CD3&!R](=[OD1&!R])[#6&!R]', 'acyl_amide<>[C,c][C;!R](=O)[N;!R][C;!R](=O)[C,c]',
        'acyl_cyanides<>[NX1]#[CX2&!R][CX3&!R](=[OX1])[#6]', 'aldehyde<>[*][CX3&!RH1]=[OD1]',
        'alkyl_halide_noF<>[CX4&!R][Br,I,Cl]', 'alkynes<>[#6][#6]#[#6][#6,#1]', 'alpha_haloamines<>[F,Cl,I,Br][#7]',
        'alphahalo_ketone_carbonyl<>[#6][CD3&!c](=[OD1])[CX4&!c][Cl,Br,I,F]',
        'anhydride<>[OX1]=[CD3]([*])[OD2][CD3](=[OX1])[*]',
        'aziridine<>N1[$([CX4H2]),$([CX4H1][#6]),$([CX4]([#6])[#6])][$([CX4H2]),$([CX4H1][#6]),$([CX4]([#6])[#6])]1',
        'azo<>[#6][NX2]=[NX2][#6]', 'azocyanamide<>[N;R0]=[N;R0]C#N',
        'Benzoxazinones<>[#6]-1=[ND2&!R]-[#6]=[C]-[#6](=[OD1&!R])-[OD2&!R]-1',
        'beta_heterosubstituted_carbonyl<>[OD1&!R]=[$([CD3&!R]);!$([CD3&!R][#8])]CC([F,Br,I,Cl])',
        'betalactams<>[OX1&!R]=C1CCN1', 'carbamate_thiocarbamate<>[NX3][CX3&!R](=[OX1,SX1])[OX2&!R,#16&!R][#6]',
        'carbamic_acid<>[$([NX4+]),$([ND1&!R]),$([NX3&!H2]),$([NX3&!RH1][#6]),$([NX3&!R]([#6])[#6])][CD3&!R](=[OD1&!R])[OX2&!RH1,OX1-]',
        'carbazide<>O=*N=[N+]=[N-]', 'carbodiimide<>[#6][NX2]=C=[NX2][#6]', 'chloramidine<>[Cl]C([C&R0])=N',
        'coumarines<>[OX1&!R]=c1ccc2ccccc2o1', 'cyanohydrins<>[ND1&!R]#[CD2&!R][CX4&!R][OD2H1]',
        'cyanophosphonate<>P(OC)(OC)(=O)C#N', 'cyclopropylamine<>C1CC1[NX3&!R]', 'diazonium<>[c,C][N+&!R]#[N&!R]',
        'enamine<>[CX3]=[CX3&!R][NX3&!R]', 'epoxide<>[OX2r3]1[#6r3][#6r3]1',
        'formic_acid_esters<>[CX3H1&!R](=[OX1])[OX2&!R][#6]',
        'halo_alkene<>[Br,F,I,Cl][$([CX3&!R][#6]),$([CX3&!RH1])]=[$([CX3&!R]([#6])[#6]),$([CX3&!RH1][#6]),$([CX3&!RH2])]',
        'halo_amine<>[CX4&!R][$([NX3&!RH1]),$([NX3&!R][#6])][Br,F,I,Cl]', 'halopyrimidine<>c1cnc([F,Cl,Br,I])nc1',
        'hemiaminal<>[OX2&!RH1][$([CX4&!RH2]),$([CX4&!RH1][#6]),$([CX4&!R]([#6])[#6])][$([NX3&!RH2]),$([NX3&!RH1][#6]),$([NX3&!R]([#6])[#6])]',
        'hemiketal<>[OX2&!RH1][$([CX4&!RH2]),$([CX4&!RH1][#6]),$([CX4&!R]([#6])[#6])][$([OX2&!RH1]),$([OX2&!R][#6])]',
        'heteroatom_heteroatom_N_N<>[$([#6][NH1;!R]-[NH1;!R][#6]),$([#6][NH1;!R]-[N;!R]([#6])[#6]),$([#6][N;!R]([#6])-[N;!R]([#6])[#6])]',
        'heteroatom_heteroatom_N_S<>[$([#6][NH1;!R]-[SD2;!R][#6]),$([#6][N;!R]([#6])-[SD2;!R][#6])]',
        'heteroatom_heteroatom_O_N<>[$([#6][#8;!R][NH1;!R][#6]),$([#6][#8;!R][N;!R]([#6])[#6])]',
        'heteroatom_heteroatom_S_O<>[#6][SD2;!R]-[#8;!R][#6]', 'heteroatom_heteroatom_S_S<>[#6][SD2;!R][SD2;!R][#6]',
        'hydantoin<>[OX1]=C1[NX3H1]C(=[OX1])CN1', 'hydrazide<>[NX3&!R][NX3&!R][#6X3](=[OX1])[#6]',
        'hydrazine<>[#6][NX3&!RH1][NX3&!RH2]', 'hydrazone<>[CX3&!R]=[NX2&!R][NX3&!R]',
        'hydroxamic_acid<>[CX3&!R;$([H0][#6])](=[OD1])[#7X3;$([H1])][OX2&!RH1,OX1-]',
        'hydroxylamine<>[c,C;!$([#6]=[#8])][$([NX3&!RH1]),$([NX3&!R][#6])][OX2&!RH1]',
        'imide<>[C,c][C;!R](=O)[N;!R][C;!R](=O)[C,c]',
        'imidoyl_halide<>[$([CX3&!RH1]),$([CX3&!RH0][#6])](=[NX2&!R])[F,Br,I,Cl]', 'imine<>[CX3]=[NX2!R]',
        'isocyanate<>[OD1&!R]=[CD2&!R]=[ND2&!R][#6]', 'isocyanide_isonitrile<>[#6][#7+]#[#6-]',
        'isothiocyanate<>[SD1&!R]=[CD2&!R]=[ND2&!R][#6]', 'lawesson_reagent_derivative<>[#6]P1(=S)SP(S1)(=S)[#6]',
        'maleimide<>C1=CC(=[OX1])NC1=[OX1]',
        'michael_acceptors_double<>[C]=[C][$([#6]=[OX1&!R]),$([#6][#7&!R]),$([#6][#16!R])]',
        'michael_acceptors_triple<>[#6]#[#6][$([#6]=[OX1&!R]),$([#6][#7&!R]),$([#6][#16!R])]',
        'nitrile<>[ND1&!R]#[CD2!R][#6,#7,#16]', 'o_thiocarbamate<>[#6][OX2&!R][CX3&!R](=[SX1])[NX3&!R][#6]',
        'orthoquinone<>[OD1&!R]=C1C(=[OD1&!R])C=CC=C1', 'oxime<>[CX3]=[NX2][OX2H1]',
        'para_hydroquinone<>[OH]c1ccc([OH])cc1', 'stilbene<>c1ccc([#6&!R]=[#6&!R]c2ccccc2)cc1',
        'pentafluorophen_ester<>C(=O)(C)Oc1c(F)c(F)c(F)c(F)c1(F)',
        'perhaloketone<>[#6][CD3&!R](=[OD1&!R])[CD4&!R]([F,Br,I,Cl])([F,Br,I,Cl])[F,Br,I,Cl]',
        'peroxide<>[#6][#8][#8][#6]', 'phosphorane<>C=P',
        'phosphoric_acid_or_ester<>[#6,#1][OX2&!R][PX4&!R](=[OX1])([OX2&!R][#6,#1])[OX2&!R][#6,#1]',
        'polyenes<>[#6]=[CX3&!R][CX3&!R]=[$([CX3&!R][#6&!R]),$([CX3&!RH2])]', 'propiolactone<>[OX1&!R]=C1OCC1',
        'propiosultone<>[OX1&!R]=[SX4]1(=[OX1&!R])OCCC1', 'quinone<>O=[#6]1[#6]:,=[#6][#6](=O)[#6]:,=[#6]1',
        's_thiocarbamate<>[#6][SX2&!R][CX3&!R](=[OX1])[NX3&!R][#6]',
        'sulfone<>[$([#16X4&!R](=[OX1])(=[OX1])([#6])[#6]),$([#16X4+2&!R]([OX1-])([OX1-])([#6])[#6])]',
        'sulfonic_acid_ester<>[#6][SX4&!R](=[OX1])(=[OX1])[$([OX2&!R][#6,#1]),$([OX1-])]', 'sulfonium<>[SX3+]',
        'sulfonyl_cyanide<>[#6][SD4&!R,SD4+2&!R](=[OD1&!R])(=[OD1&!R])[CD2&!R]#[NX1]',
        'sulfonyl_halide<>[#6,#7][#16]([F,Cl,I,Br])(=[#8])=[#8]',
        'sulfonyl_urea<>[NX3H1][CX3](=[OX1])[NX3H1][SX4](=[OX1])=[OX1]',
        'sulfoxide<>[$([#16X3](=[OX1])([#6])[#6]),$([#16X3+]([OX1-])([#6])[#6])]', 'sulphanylamino<>[#6][NX3H1][SX2H1]',
        'terminal_vinyl<>[CX3H1&!R]=[CX3H2&!R]', 'thioacetal<>[#6][CX4&!R]([#16&!R][#6])[#16&!R][#6]',
        'thioamide<>[#6][CX3&!R](=[SX1&!R])[NX3]([#6])[#6]', 'thiocarbonyl_aromatic<>[C,c]=[SX1&!RH0]',
        'thioester<>[OX1&!R]=[#6]([#6])[SX2]',
        'thioic_acid<>[$([CX3]([OX2H1])=S),$([CX3]([SX2H1])=O),$([CX3]([SX2H1])=S)]', 'thioketone<>[SD1&!R]=[#6]',
        'thiol<>[#6][SH1X2&!R]',
        'thiourea<>[$([NX3&!RH2]),$([NX3&!RH1][#6]),$([NX3&!R]([#6])[#6])][CX3&!R](=[SX1&!R])[$([NX3&!RH2]),$([NX3&!RH1][#6]),$([NX3&!R]([#6])[#6])]',
        'triflate<>OS(=O)(=O)C(F)(F)F', 'boronic_acid_ester<>[#6]-[BX3](-[O])-[O]',
        'vinyl_sulfone<>[#6,#7][SX4](=[OX1])(=[OX1])[C]=[C]',
        'oxodiazolone_1<>[#8X1&!R]=[#6]-1-[#7]-[#6]=[#7X2]-[#8X2]-1',
        'cyclic_acetal<>[$([#6r5]),$([#6r6])]-[$([#8r5]),$([#8r6])]-[$([CX4H2]),$([CX4H1][#6]),$([CX4]([#6])[#6])]-[$([#8r5]),$([#8r6])]-[$([#6r5]),$([#6r6])]',
        'sulfonyl_pyrazole<>[#6]-1=[#7]-[#7](-[#6]=[#6]-1)[S&X4&!R]([#6])(=[O&X1])=[O&X1]',
        'oxodiazolone_2<>[#8X1&!R]=[#6]-1-[#7]-[#7X2]=[#6]-[#8X2]-1',
        'sulfonyl_pyrazole_2<>[#7X2]-1=[#6!R2]-[#7](-[#6!R2]=[#6!R2]-1)[S&X4&!R]([#6])(=[O&X1])=[O&X1]',
        'oxazolone<>[#8X1&!R]=[#6]-1-[#6R1]-[#6R1]=[#7X2]-[#8X2]-1',
        'oxodiazolone_3<>[#8X1&!R]=[#6]-1-[#7X2]=[#6]-[#7]-[#8X2]-1',
        'carboxamide_1<>[$([#7H]-[#6]),$([#7](-[#6])[#6])]-[#6X3&!R](=[OX1&!R])-[#7]-1-[#6,#7]=[#6,#7]-[#7,#6]=[#6,#7]-1',
        'tertiary_halo<>[#6]C([#6])([#6])[Br,I,Cl]',
        'thiazolidinone_ene<>[#6]=[#6]-1-[#16R1X2]-[#6]-[#7]-[#6X3R1]-1=[OX1&!R]',
        'thiazolidinedione<>S1C(=[OX1])[N]C(=[OX1])C1', '1_2_thiazol_3_one<>S1[NX3H1]C(=[OX1])C=C1',
        'thiazolidinone<>[OX1&!R]=[#6]-1-[#7]-[#6]-[#6]-[#16]-1', '2-chloroheteroN<>Cl-[#6&R1]-[#7]',
        '2-chloroheteroN_2<>Cl-[#6&R1]=[#7]']
    for line in cov_group_list:
        line_splitted = str(line).split('<>')
        name = line_splitted[0]
        smart = line_splitted[1]
        cov_smarts[name] = smart.strip('\n')
        cov_stats[name] = 0
    cov_stats['unknown'] = 0
    return cov_smarts, cov_stats


def find_covgroup(mol_list):
    new_list = []
    bar = progressbar.ProgressBar(maxval=len(mol_list), widgets=widgets)
    bar.start()
    counter = 0
    for mol in mol_list:
        counter += 1
        bar.update(counter)
        if mol is None:
            continue
        if not args.only_stat:
            Chem.Kekulize(mol)
            list_of_match = []
            # if mol.GetProp('Name') !='Z1378215196':
            #     continue
            for group_name, pattern in cov_smarts.items():
                # if group_name != '':
                #     continue
                subst = Chem.MolFromSmarts(pattern)
                if mol.HasSubstructMatch(subst):
                    list_of_match.append(group_name)
            group_list = ''

            if len(list_of_match) > 1 or '2-chloroheteroN_2' in list_of_match:
                for i in reversed(range(len(list_of_match))):
                    if list_of_match[i] == 'carbamate_thiocarbamate':
                        list_of_match.pop(i)
                        continue
                    if list_of_match[i] == 'terminal_vinyl' and 'michael_acceptors_double' in list_of_match:
                        list_of_match.pop(i)
                        continue
                    if list_of_match[i] == 'terminal_vinyl' and 'vinyl_sulfone' in list_of_match:
                        list_of_match.pop(i)
                        continue
                    if list_of_match[i] == 'polyenes' and 'michael_acceptors_double' in list_of_match:
                        list_of_match.pop(i)
                        continue
                    if list_of_match[i] == 'polyenes' and 'acrylamide' in list_of_match:
                        list_of_match.pop(i)
                        continue
                    if list_of_match[i] == 'sulfone' and 'vinyl_sulfone' in list_of_match:
                        list_of_match.pop(i)
                        continue
                    if list_of_match[i] == 'nitrile' and 'michael_acceptors_double' in list_of_match:
                        subst = Chem.MolFromSmarts('[C]=[C](C#N)[$([#6]=[OX1&!R]),$([#6][#7&!R]),$([#6][#16!R])]')
                        if mol.HasSubstructMatch(subst):
                            list_of_match.pop(i)
                            continue
                    if list_of_match[i] == 'nitrile' and 'vinyl_sulfone' in list_of_match:
                        subst = Chem.MolFromSmarts('[#6,#7][SX4](=[OX1])(=[OX1])[C](C#N)=[C]')
                        if mol.HasSubstructMatch(subst):
                            list_of_match.pop(i)
                            continue
                    if list_of_match[i] == 'michael_acceptors_double' and 'acrylamide' in list_of_match:
                        list_of_match.pop(i)
                        continue
                    if list_of_match[i] == 'michael_acceptors_double' and 'vinyl_sulfone' in list_of_match:
                        list_of_match.pop(i)
                        continue
                    if list_of_match[i] == 'acyclic_ketone' and 'michael_acceptors_double' in list_of_match:
                        list_of_match.pop(i)
                        continue
                    if str(list_of_match[i]).startswith('2-chloroheteroN_'):
                        list_of_match[i] = '2-chloroheteroN'
                        continue
            list_of_match = list(set(list_of_match))
            for i in list_of_match:
                group_list += i + ', '
            group_list = group_list.rstrip(', ')
            if group_list == '':
                group_list = 'unknown'
            mol.SetProp('Covalent_group', group_list)
        mol.SetProp('MW', str(Desc.MolWt(mol)))
        mol.SetProp('LogP', str(Desc.MolLogP(mol)))
        mol.SetProp('TPSA', str(Desc.TPSA(mol)))
        mol.SetProp('HBA', str(Desc.NumHAcceptors(mol)))
        mol.SetProp('HBD', str(Desc.NumHDonors(mol)))
        new_list.append(mol)
    bar.finish()
    return new_list


def generate_stat(new_list):
    descriptors_stat = {'TPSA': [], 'MW': [], 'LogP': [], 'HBA': [], 'HBD': []}

    if not args.no_entropy:
        entropy = get_entropy_measures(new_list)
        chem_coverage_score = get_chemical_coverage_score(new_list)

    for m in new_list:
        descriptors_stat['TPSA'].append(float(m.GetProp('TPSA')))
        descriptors_stat['MW'].append(float(m.GetProp('MW')))
        descriptors_stat['LogP'].append(float(m.GetProp('LogP')))
        descriptors_stat['HBA'].append(float(m.GetProp('HBA')))
        descriptors_stat['HBD'].append(float(m.GetProp('HBD')))
        if not args.only_stat:
            for group_name in cov_stats.keys():
                if group_name in m.GetProp('Covalent_group').split(', '):
                    cov_stats[group_name] += 1

    labels = []
    values = []
    labels_other = []
    values_other = []
    labels_other_less = []
    values_other_less = []
    temp_other_less = 0
    temp_other = 0

    if not args.only_stat:
        cov_stats_sorted = sorted(cov_stats.items(), key=operator.itemgetter(1), reverse=True)
        for keys, val in cov_stats_sorted:
            if keys != 'other':
                if (val / (len(new_list))) <= 0.015:
                    if val >= 10:
                        labels_other.append(keys)
                        values_other.append(val)
                    elif 0 < val < 10:
                        labels_other_less.append(keys)
                        values_other_less.append(val)
                        temp_other_less += val
                    temp_other += val
                elif val > 0:
                    labels.append(keys)
                    values.append(val)
        labels.append('other')
        values.append(temp_other)
        labels_other.append('Others(less than 10)')
        values_other.append(temp_other_less)

    descriptors_stat['TPSA'] = np.array(descriptors_stat['TPSA'])
    descriptors_stat['MW'] = np.array(descriptors_stat['MW'])
    descriptors_stat['LogP'] = np.array(descriptors_stat['LogP'])
    descriptors_stat['HBA'] = np.array(descriptors_stat['HBA'])
    descriptors_stat['HBD'] = np.array(descriptors_stat['HBD'])

    if not args.only_pie:
        gs = gridspec.GridSpec(8, 2)
        props = dict(boxstyle='round', facecolor='wheat')

        plt.figure(1, (20, 20))
        plt.subplot(gs[0, 0])
        plt.xlabel('TPSA')
        plt.ylabel('TPSA', fontsize=20, weight='bold', bbox=props)
        # plt.title('TPSA binned distribution')
        plt.hist(descriptors_stat['TPSA'], 20, color='g')

        plt.subplot(gs[0, 1])
        plt.xlabel('TPSA')
        # plt.title('TPSA boxed')
        plt.boxplot(descriptors_stat['TPSA'], 0, 'rs', 0, 0.98)

        plt.subplot(gs[1, 0])
        plt.xlabel('MW')
        plt.ylabel('MW', fontsize=20, weight='bold', bbox=props)

        # plt.title('MW binned distribution')
        plt.hist(descriptors_stat['MW'], 30, color='g')

        plt.subplot(gs[1, 1])
        plt.xlabel('MW')
        # plt.title('MW boxed')
        plt.boxplot(descriptors_stat['MW'], 0, 'rs', 0, 0.98)

        plt.subplot(gs[2, 0])
        plt.xlabel('LogP')
        plt.ylabel('LogP', fontsize=20, weight='bold', bbox=props)

        # plt.title('LogP binned distribution')
        plt.hist(descriptors_stat['LogP'], 20, color='g')

        plt.subplot(gs[2, 1])
        plt.xlabel('LogP')
        # plt.title('LogP boxed')
        plt.boxplot(descriptors_stat['LogP'], 0, 'rs', 0, 0.98)

        plt.subplot(gs[3, 0])
        plt.xlabel('HBA')
        plt.ylabel('HBA', fontsize=20, weight='bold', bbox=props)

        # plt.title('HBA binned distribution')
        plt.hist(descriptors_stat['HBA'], 6, color='g')

        plt.subplot(gs[3, 1])
        plt.xlabel('HBA')
        # plt.title('HBA boxed')
        plt.boxplot(descriptors_stat['HBA'], 0, 'rs', 0, 0.98)

        plt.subplot(gs[4, 0])
        plt.xlabel('HBD')
        plt.ylabel('HBD', fontsize=20, weight='bold', bbox=props)
        # plt.title('HBD binned distribution')
        plt.hist(descriptors_stat['HBD'], 3, color='g')

        plt.subplot(gs[4, 1])
        plt.xlabel('HBD')
        # plt.title('HBD boxed')
        plt.boxplot(descriptors_stat['HBD'], 0, 'rs', 0, 0.98)

        if not args.no_entropy:
            ax = plt.subplot(gs[5, :])

            plt.title('Similarity-based (Tanimoto>0.75) entropy score', fontsize=20, weight='bold', bbox=props)
            # ax1 = ax.add_axes([0.05, 0.8, 0.9, 0.15], frameon=True)
            rect = mpl.patches.Rectangle((entropy['SSE_similarity'] + 0.005, 0), 1, 0.01, color='b', angle=90.0)
            ax.add_patch(rect)
            ax.annotate(round(entropy['SSE_similarity'], 3), ((entropy['SSE_similarity'], 0.5)), color='b',
                        weight='bold',
                        fontsize=20, ha='center', va='center', bbox=props)
            ax.annotate('No diversity', (0, 0.9), weight='bold',
                        fontsize=12, ha='center', va='center', bbox=props)
            ax.annotate('high diversity', (1, 0.9), weight='bold',
                        fontsize=12, ha='center', va='center', bbox=props)
            cmap = mpl.cm.RdYlGn
            norm = mpl.colors.Normalize(vmin=0, vmax=1)

            cb1 = mpl.colorbar.ColorbarBase(ax, cmap=cmap,
                                            norm=norm,
                                            orientation='horizontal')

            ax2 = plt.subplot(gs[6, :])

            plt.title('Scaffold-based (Murcko decomposition) entropy score', fontsize=20, weight='bold', bbox=props)
            rect = mpl.patches.Rectangle((entropy['SSE_scaffold'] + 0.005, 0), 1, 0.01, color='b', angle=90.0)
            ax2.add_patch(rect)
            ax2.annotate(round(entropy['SSE_scaffold'], 3), ((entropy['SSE_scaffold'], 0.5)), color='b', weight='bold',
                         fontsize=20, ha='center', va='center', bbox=props)

            ax2.annotate('No diversity', (0, 0.9), weight='bold',
                         fontsize=12, ha='center', va='center', bbox=props)
            ax2.annotate('High diversity', (1, 0.9), weight='bold',
                         fontsize=12, ha='center', va='center', bbox=props)
            cmap = mpl.cm.RdYlGn
            norm = mpl.colors.Normalize(vmin=0, vmax=1)

            cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cmap,
                                            norm=norm,
                                            orientation='horizontal')
            ax3 = plt.subplot(gs[7, :])

            plt.title('Chemical space coverage score', fontsize=20, weight='bold', bbox=props)
            rect = mpl.patches.Rectangle((chem_coverage_score + 0.005, 0), 1, 0.01, color='b', angle=90.0)
            ax3.add_patch(rect)
            ax3.annotate(round(chem_coverage_score, 3), ((chem_coverage_score, 0.5)), color='b', weight='bold',
                         fontsize=20, ha='center', va='center', bbox=props)

            ax3.annotate('Low coverage', (0, 0.9), weight='bold',
                         fontsize=12, ha='center', va='center', bbox=props)
            ax3.annotate('High coverage', (1, 0.9), weight='bold',
                         fontsize=12, ha='center', va='center', bbox=props)
            cmap = mpl.cm.RdYlGn
            norm = mpl.colors.Normalize(vmin=0, vmax=1)

            cb1 = mpl.colorbar.ColorbarBase(ax3, cmap=cmap,
                                            norm=norm,
                                            orientation='horizontal')

        plt.subplots_adjust(top=0.92, bottom=0.08, left=0.10, right=0.95, hspace=0.95,
                            wspace=0.15)
    if not args.only_stat:
        plt.figure(2, (20, 20))
        cmap = plt.cm.Set3
        colors = cmap(np.linspace(0, 1, len(values)))
        p, tx, autotexts = plt.pie(values, shadow=True, startangle=90, colors=colors, autopct='%1.1f%%',
                                   counterclock=False)
        plt.setp(autotexts, size=8, weight="bold")
        for i, a in enumerate(autotexts):
            a.set_text("{}".format(str(a.get_text()) + '\n(' + str(values[i]) + ')'))
        plt.legend(labels, loc="best",fontsize='medium')
        plt.axis('equal')
        plt.title('Covalent groups repartition')

        plt.figure(3, (20, 20))
        plt.title('Expanded view for "others"')
        cmap = plt.cm.Set3
        colors = cmap(np.linspace(0, 1, len(values_other)))
        p, tx, autotexts = plt.pie(values_other, shadow=True, startangle=90, colors=colors,
                                   autopct='%1.1f%%',
                                   counterclock=False)
        plt.setp(autotexts, size=8, weight="bold")
        for i, a in enumerate(autotexts):
            a.set_text("{}".format('(' + str(values_other[i]) + ')'))
        plt.legend(labels_other, loc="best",fontsize='medium')
        plt.axis('equal')
        text = ''
        for i in range(len(labels_other_less)):
            text += str(labels_other_less[i]) + ': ' + str(values_other_less[i]) + '\n'
        text = text.rstrip('\n')
        props = dict(boxstyle='round', facecolor='wheat')
        plt.text(-1.0, 1.0, text, verticalalignment='top', bbox=props)
    plt.show()


def bin_PCA(df, num_bins=5):
    pca1 = pd.DataFrame(data=pd.cut(df['PCA 1'], num_bins, labels=False).values, columns=['PCA1_bin'])
    pca2 = pd.DataFrame(data=pd.cut(df['PCA 2'], num_bins, labels=False).values, columns=['PCA2_bin'])
    pca3 = pd.DataFrame(data=pd.cut(df['PCA 3'], num_bins, labels=False).values, columns=['PCA3_bin'])
    pca4 = pd.DataFrame(data=pd.cut(df['PCA 4'], num_bins, labels=False).values, columns=['PCA4_bin'])
    pca5 = pd.DataFrame(data=pd.cut(df['PCA 5'], num_bins, labels=False).values, columns=['PCA5_bin'])
    df = df.join(pca1).join(pca2).join(pca3).join(pca4).join(pca5)
    df['Bin_ID'] = df[['PCA1_bin', 'PCA2_bin', 'PCA3_bin', 'PCA4_bin', 'PCA5_bin']].apply(
        lambda x: ''.join(map(str, x)), axis=1)
    return df


def get_chemical_coverage_score(mol_list):
    fps = []
    for m in mol_list:
        arr = np.zeros((1,))
        f = AllChem.GetMorganFingerprintAsBitVect(m, 3)
        DataStructs.ConvertToNumpyArray(f, arr)
        fps.append(arr)
    fps_array = np.asarray(fps)
    pca = PCA(n_components=5)
    pca.fit(fps_array)
    Xred = pca.transform(fps_array)
    PCA_df = pd.DataFrame(Xred, columns=['PCA 1', 'PCA 2', 'PCA 3', 'PCA 4', 'PCA 5'])
    frame_PCA = bin_PCA(PCA_df, num_bins=5)
    score = frame_PCA.Bin_ID.nunique() / (6 ** 5)
    return score


def get_entropy_measures(mol_list):
    list_cores = []
    fps_list = []
    for m in mol_list:
        core_smi = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
        list_cores.append(core_smi)
        fps_list.append(AllChem.GetMorganFingerprintAsBitVect(m, 3))

    unique_cores = list(set(list_cores))
    core_count = []
    for i in unique_cores:
        counter = 0
        for j in list_cores:
            if i == j:
                counter += 1
        core_count.append((i, counter))
        core_count = sorted(core_count, key=operator.itemgetter(1), reverse=True)

    top_12_core = core_count[:12]
    mol_core = []
    legend = []
    for core in top_12_core:
        m = Chem.MolFromSmiles(core[0])
        AllChem.Compute2DCoords(m)
        mol_core.append(m)
        ratio = str(round(core[1] / len(mol_list) * 100, 1))
        legend.append('Count = ' + str(core[1]) + ' (' + ratio + '%)')

    img = Draw.MolsToGridImage(mol_core, molsPerRow=4, subImgSize=(200, 200), legends=legend)
    img.save(args.output + '_top_cores.png')

    s_entropy = 0
    for i in core_count:
        prob_core = i[1] / len(list_cores)
        sum_term = prob_core * np.log(prob_core)
        s_entropy += sum_term
    s_entropy = -s_entropy
    ss_entropy = s_entropy / (np.log(len(unique_cores)))

    # time_start = t.time()
    list_of_dist = np.zeros((len(fps_list), len(fps_list)))

    for i in range(len(fps_list)):
        for j in range(i, len(fps_list)):
            if list_of_dist[i][j] == 0 and list_of_dist[j][i] == 0 and i != j:
                sim = DataStructs.FingerprintSimilarity(fps_list[i], fps_list[j])
                list_of_dist[i][j] = sim
                list_of_dist[j][i] = sim
            elif i == j:
                list_of_dist[i][j] = 1

    # time_stop = t.time()
    # print('It took : ', str(time_stop - time_start), ' seconds to calculate the matrix')

    s_entropy_fps = 0
    already_clustered = []
    for i in range(len(list_of_dist)):
        counter = 0
        for j in range(i, len(list_of_dist)):
            if list_of_dist[i][j] >= 0.75 and j not in already_clustered:
                counter += 1
                already_clustered.append(j)
        if counter != 0:
            prob = counter / len(fps_list)
            sum_term = prob * np.log(prob)
            s_entropy_fps += sum_term
    s_entropy_fps = -s_entropy_fps
    ss_entropy_fps = s_entropy_fps / (np.log(len(fps_list)))

    return {'SE_scaffold': s_entropy, 'SSE_scaffold': ss_entropy, 'SE_similarity': s_entropy_fps,
            'SSE_similarity': ss_entropy_fps}


if __name__ == "__main__":
    cov_smarts, cov_stats = dict_of_cov_group()

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_file', help='Name of the input file [MANDATORY]', metavar='')
    parser.add_argument('-o', '--output',
                        help='If you want an SDF with the annotated molecule please use this option followed by the name of the output file',
                        metavar='')
    parser.add_argument('-s', '--no_stat',
                        help='Do not display the stats',
                        action='store_true', default=False)
    parser.add_argument('-p', '--only_pie',
                        help='Only display the pie chart',
                        action='store_true', default=False)
    parser.add_argument('-S', '--only_stat',
                        help='Only display the library stat no covalent group search',
                        action='store_true', default=False)
    parser.add_argument('-e', '--no_entropy',
                        help='do not calculate entropy',
                        action='store_true', default=False)
    args = parser.parse_args()

    if args.input_file:
        if os.path.exists(args.input_file):
            mol_list = Chem.SDMolSupplier(args.input_file)
        else:
            print('ERROR : file ' + str(args.input_file) + ' does not exist')
            sys.exit()
    else:

        parser.print_help()
        print('\n', '/!\============= ERROR : Please use the [-i] option to input your file ===============/!\\')
        sys.exit()

    # mol_list = Chem.SDMolSupplier('data/Enamine_Covalent_Fragments_3038cmpds.sdf')

    # mol_list = Chem.SDMolSupplier('data/Enamine_Cys_focused_Covalent_Fragments_995cmpds.sdf')

    new_list = find_covgroup(mol_list)

    if args.no_stat:
        pass
    else:
        generate_stat(new_list)

    if args.output:
        w = Chem.SDWriter(args.output + '.sdf')
        for m in new_list:
            w.write(m)
