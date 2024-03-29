#!/usr/bin/env python
# ----------------------------------------------------------------------#
# 			CNS MPO Score and Solubility forecaster index (3.2)			#
# --------------------------- Compatibility ----------------------------#
# 					Compatible with Python 2.7 & 3.5					#
# ----------------------Prerequiste installations-----------------------#
# 																		#
# 		- Chemaxon Marvin suite WITH license							#
# 				cxcalc module required to work 	(used for pKa pred)		#
# 					(make sure to add cxcalc to $PATH)					#
# by default : C:\Program Files (x86)\ChemAxon\MarvinBeans\bin			#
# How to use :															#
# 		- input the sdf file name when requested (.sdf included)		#
# 		- the script output an sdf with #name_out.sdf					#
# 			This sdf contains the fields with : 						#
# 						- CNS MPO Score									#
# 						- Solubility forecaster index (SFI)[optional]	#
# 						- bpKa,logD(7.4),logP,MW,HBD,#Ar,TPSA			#
# 						- All individual components of MPO Score		#
# ----------------------------------------------------------------------#
# (C) Dr. De Cesco Stephane - v3.2 - 22/03/2017							#
# ----------------------------------------------------------------------#
# -------------------------- MODULES IMPORT ----------------------------#
from __future__ import print_function

import argparse
import os
import os.path
import subprocess
import sys
import math
import rdkit.Chem as Chem

# To remove danger of using input in python2 environments.
try:
	input = raw_input
except NameError:
	pass


# ----------------------------- Desc object -----------------------------#

def monotonic_score(value, lower, upper):
	# =======================================================================================
	#         upper
	# 1|---------
	#  |		  \
	#  |           \@value
	#  |			\
	# 0|_____________\______________________
	# 	             lower
	# Function to return a score between 0 and 1 depending of the value of the parameters.
	# 		 | 1 if value < upper
	# Score ={ f(upper and lower) if upper < value < lower
	# 		 | 0 if value > lower
	# =======================================================================================
	try:
		value = float(value)
	except ValueError as message:
		print(message)
		print(value)
	upper = float(upper)
	lower = float(lower)
	if value <= lower:
		score = 1
	elif value > upper:
		score = 0
	else:
		score = 1 - ((value - lower) * (1 / (upper - lower)))
	return score


def hump_score(value, low1, up1, up2, low2):
	# =======================================================================================
	#         	up1		 up2
	# 1|		  --------
	#  |		 / 		  \
	#  |        /      	   \
	#  |	   /			\
	# 0|______/______________\_______________
	# 	     low1			  low2
	# Function to return a score between 0 and 1 depending of the value of the parameters.
	# 		 | 0 if value < low1
	# 		 | f (low1 and up1) if low1 < value < up1
	# Score ={ 1 if up1 < value < up2
	# 		 | f (up2 and low2) if up2 < value < low2
	# 		 | 0 if value > lower
	# =======================================================================================

	value, low1, up1, up2, low2 = float(value), float(low1), float(up1), float(up2), float(low2)
	score = 0
	if value <= low1:
		score = 0
	elif up1 < value <= up2:
		score = 1
	elif value > low2:
		score = 0
	elif low1 < value <= up1:
		score = ((value - low1) * (1 / (up1 - low1)))
	elif up2 < value <= low2:
		score = 1 - ((value - up2) * (1 / (low2 - up2)))
	return score


def calc_BBBScore(prop_dict):
	weight_dict = {'Aro_R': 1, 'HA': 1, 'MWHBN': 1.5, 'TPSA': 2, 'pKa': 0.5}
	scores_dict = {}

	scores_dict['Aro_R'] = calc_AroScore(int(prop_dict['ArRings']))
	scores_dict['HA'] = calc_HAScore(int(prop_dict['HAC']))
	scores_dict['MWHBN'] = calc_MWHBNScore(float(prop_dict['MW']), float(prop_dict['HBD']),float(prop_dict['HBA']))
	scores_dict['TPSA'] = calc_TPSAScore(float(prop_dict['TPSA']))
	scores_dict['pKa'] = calc_PKAScore(float(prop_dict['bpKa']))

	total_score = 0

	for key,value in scores_dict.items():
		total_score+=value*weight_dict[key]

	return total_score


def calc_AroScore(nAr):
	if nAr > 4:
		return 0
	elif nAr == 4:
		return 0.199399
	elif nAr == 3:
		return 0.691115
	elif nAr == 2:
		return 1
	elif nAr == 1:
		return 0.816016
	elif nAr == 0:
		return 0.336376


def calc_HAScore(HA):
	if 5 < HA <= 45:
		return (1 / 0.624231) * (0.0000443 * (HA ** 3) - 0.004556 * (HA ** 2) + 0.12775 * HA - 0.463)
	else:
		return 0


def calc_MWHBNScore(mw, hbd,hba):
	mwhbn = (hbd + hba) / math.sqrt(mw)
	if 0.05 < mwhbn <= 0.45:
		return (1 / 0.72258) * (26.733 * (mwhbn ** 3) - 31.495 * (mwhbn ** 2) + 9.5202 * mwhbn - 0.1358)
	else:
		return 0


def calc_TPSAScore(tpsa):
	if 0 < tpsa <= 120:
		return (1 / 0.9598) * (-0.0067 * tpsa + 0.9598)
	else:
		return 0


def calc_PKAScore(pka):
	if 3 < pka <= 11:
		return (1 / 0.597488) * (0.00045068 * (pka**4) - 0.016331 * (pka ** 3) + 0.18618 * (pka ** 2) - 0.71043 * pka + 0.8579)
	else:
		return 0

class Desc:
	def __init__(self, mol):
		self.mol = mol  # stores the mol coordinate into .mol
		self.Prop = {}  # Create an empty dictionnary to store all the individual Properties
		self.get_param()
		if args.sfi:
			self.Prop['SFI'] = float(self.Prop['logD']) + float(self.Prop['ArRings'])  # Calc SFI
		self.calc_mpo_score()  # Calculate the MPO score

	# self.calc_mpo_area()

	def get_param(self):
		record_data = False
		coordinates = True
		new_mol = ''
		transcription = {'LOGD': 'logD', 'LOGP': 'logP', 'PSA': 'TPSA', 'MASS': 'MW', 'AROMATIC_RINGCOUNT': 'ArRings',
		                 'pkacalculator': 'bpKa', 'DONOR_COUNT': 'HBD','ACCEPTOR_COUNT':'HBA', 'CHARGE_DISTRIBUTION': 'charge_dist(7.4)',
		                 'ATOMCOUNT': 'HAC'}
		lines = self.mol.splitlines()

		for i in range(len(lines)):
			if lines[i].startswith('>  <'):
				c_lines = lines[i].lstrip('>  <')
				chemprop = c_lines.rstrip('>')
				if chemprop in transcription.keys() or chemprop in transcription.values():
					coordinates = False
			if coordinates:
				new_mol += lines[i] + '\n'
			if record_data:
				value = str(lines[i])
				if chemprop == 'LOGD':
					value = lines[i].split('\t')[1].strip()
				if chemprop == 'pkacalculator':
					value = lines[i].split('\t')[0]
					if value.strip() == "":
						value = 0.0
				if chemprop == 'ATOMCOUNT':
					value = int(lines[i]) - int(lines[i + 1])
				if chemprop == 'CHARGE_DISTRIBUTION':
					value = lines[i].split('\t')[1].lstrip()
				if chemprop == 'Name':
					chemprop = chemprop.lower()
				if chemprop in transcription.keys():
					self.Prop[transcription[chemprop]] = value
				elif chemprop in transcription.values():
					self.Prop[chemprop] = value
				record_data = False
			if lines[i].startswith('>  <'):
				lines[i] = lines[i].lstrip('>  <')
				chemprop = lines[i].rstrip('>')
				record_data = True
		self.mol = new_mol

	"""
	Calc_Param()
	Function deprecated
	"""

	# def Calc_Param(self):
	#     with open('tmp.sdf', 'w') as tempfile:
	#         tempfile.write(self.mol)  # create and write the molecule in a temp file
	#     if os.name == 'nt':
	#
	# command_1 = 'cxcalc -i name -S /ssddata/sdecesco/data/Scripts/CNS_MPO/data/nuak1.sdf pkacalculator -t basic -b
	# 1 logD -H 7.4 logP donorcount polarsurfacearea mass aromaticringcount'
	#
	#         res = subprocess.check_output(["cmd","/C cxcalc aromaticringcount mass polarsurfacearea donorcount pKa
	# logD -H 7.4 logP tmp.sdf"])  # using ChemAxon cxcalc command line tool to measure the pKa and logD7.4
	#         res = res.replace('\r', '')  # results parsing
	#         res_lines = res.split('\n')  # results parsing
	#         title_line = res_lines[0].split('\t')  # first line is the parameters name
	#         result_line = res_lines[1].split('\t')  # second line is with the calculated parameters
	#     elif os.name == 'posix':
	#         command = 'cxcalc aromaticringcount mass polarsurfacearea donorcount pKa logD -H 7.4 logP '
	#         command += os.getcwd() + '/tmp.sdf'
	#         res = subprocess.check_output([command], shell=True)
	#         res = res.decode("utf-8")
	#         res_lines = res.split('\n')  # results parsing
	#         title_line = res_lines[0].split('\t')  # first line is the parameters name
	#         result_line = res_lines[1].split('\t')  # second line is with the calculated parameters
	#     else:
	#         print('Your operating system is not supported... now closing')
	#         sys.exit()
	#     results = dict(zip(title_line, result_line))  # zip them and store them in a dict to retrieve them easily
	#     self.Prop['bpKa'] = results['bpKa1']  # retrieving the pKa value of the basic atom
	#     if self.Prop['bpKa'] == '':
	#         self.Prop['bpKa'] = 0
	#     self.Prop['apKa'] = results['apKa1']
	#     self.Prop['logD'] = results['logD[pH=7.4]']  # retrieving the logD7.4
	#     self.Prop['logP'] = results['logP']  # retrieving the logP
	#     self.Prop['HBD'] = results['donorcount']  # retrieving HBD count
	#     self.Prop['TPSA'] = results['Polar surface area']  # retrieving TPSA
	#     self.Prop['MW'] = results['Molecular weight']  # retrieving MW
	#     self.Prop['ArRings'] = results['Aromatic ring count']
	#     os.remove('tmp.sdf')  # deleting the temporary file

	def calc_mpo_score(
			self):  # call the monotonic or hump score function for each term with the boundaries and sum them at the
		#  end to populate the CNS MPO Score
		try:
			self.Prop['bpKaScore'] = float(monotonic_score(self.Prop['bpKa'], 8,
			                                               10))  # Todo : See with hump score for pKa (https://doi.org/10.1021/acs.jmedchem.6b01469)
			self.Prop['logPScore'] = float(monotonic_score(self.Prop['logP'], 3, 5))
			self.Prop['logDScore'] = float(monotonic_score(self.Prop['logD'], 2, 4))
			self.Prop['MWScore'] = float(monotonic_score(self.Prop['MW'], 360, 500))
			self.Prop['HBDScore'] = float(monotonic_score(self.Prop['HBD'], 0.5, 3.5))
			self.Prop['TPSAScore'] = float(hump_score(self.Prop['TPSA'], 20, 40, 90, 120))
			self.Prop['MPOScore'] = self.Prop['bpKaScore'] + self.Prop['logPScore'] + self.Prop['logDScore'] \
			                        + self.Prop['MWScore'] + self.Prop['HBDScore'] + self.Prop['TPSAScore']
			self.Prop['MPOScore_v2'] = self.Prop['bpKaScore'] + self.Prop['logPScore'] + self.Prop['MWScore'] + (
					2 * self.Prop['HBDScore']) + self.Prop['TPSAScore']
			self.Prop['BBBScore'] = calc_BBBScore(self.Prop)
			self.Prop['BBB_Aro_R'] = calc_AroScore(int(self.Prop['ArRings']))
			self.Prop['BBB_HA'] = calc_HAScore(int(self.Prop['HAC']))
			self.Prop['BBB_MWHBN'] = calc_MWHBNScore(float(self.Prop['MW']), float(self.Prop['HBD']),float(self.Prop['HBA']))
			self.Prop['BBB_TPSA'] = calc_TPSAScore(float(self.Prop['TPSA']))
			self.Prop['BBB_pKa'] = calc_PKAScore(float(self.Prop['bpKa']))

		except KeyError as missing:

			print(
				"The following parameter field [", missing,
				"] in the sdf file is missing or is not spelled as required "
				"(case sensitive):\n- \"bpKa\" (basic pKa)\n- \"logD\" ("
				"logD at pH=7.4)\n- \"logP\"\n- \"HBD\" (Hydrogen bond "
				"donor count)\n- \"TPSA\" (total polar surface area)\n- \"MW\" (molecular weight)\nMolecule skipped")
			self.Prop['bpKaScore'] = "Error"
			self.Prop['logPScore'] = "Error"
			self.Prop['logDScore'] = "Error"
			self.Prop['MWScore'] = "Error"
			self.Prop['HBDScore'] = "Error"
			self.Prop['TPSAScore'] = "Error"
			self.Prop['MPOScore'] = "Error"
			self.Prop['MPOScore_v2'] = "Error"

	def print_details(self, level):  # For visualisation in the terminal only
		if level == 0:
			print("======================")
			# noinspection PyBroadException
			try:
				print(self.Prop['name'])
			except:
				print("No name inputed")
			print('CNS MPO Score = ' + str(self.Prop['MPOScore']))
			print('BBB Score = ' + str(self.Prop['BBBScore']))
			# print('MPO percentage area = ', self.Prop['MPO_area'])
			print("======================")
		if level == 1:
			print("======================")
			# noinspection PyBroadException
			try:
				print(self.Prop['name'])
			except:
				print("No name inputed")
			for key, value in self.Prop.items():
				print(key + ": " + str(value))
			print("======================")

	def sdf_writer(self):
		mol = self.mol.replace("$$$$", '')
		mol = mol.rstrip('\n') + '\n'
		for key, value in self.Prop.items():
			mol = mol + "\n>  <" + key + ">\n"
			mol = mol + str(value) + "\n"
		mol += "\n$$$$\n"
		return mol

	def calc_mpo_area(self):
		values = [self.Prop['bpKaScore'], self.Prop['logPScore'], self.Prop['logDScore'], self.Prop['MWScore'],
		          self.Prop['HBDScore'], self.Prop['TPSAScore']]
		values.sort(reverse=True)
		areas = []
		for i in range(len(values)):
			if i == 5:
				area = ((values[i] * values[0]) / 2) * 0.866  # 0.866 = sin(60deg)
			else:
				area = ((values[i] * values[i + 1]) / 2) * 0.866  # 0.866 = sin(60deg)
			areas.append(area)
		total_area = (0.866 / 2) * 6
		fraction_area = sum(areas) / total_area
		self.Prop['MPO_area'] = round(fraction_area * 100, ndigits=2)


def calc_sdf_parser(sdf_param):
	molecule_list = []
	molecule = ''
	lines = sdf_param.splitlines(keepends=True)
	for mol_line in lines:
		molecule += mol_line
		if mol_line == "$$$$\n" or mol_line == '$$$$\r\n':
			molecule_list.append(molecule)
			molecule = ''
	return molecule_list


def sdf_parser(file):
	with open(file, 'r') as sdf_file:
		molecule_list = []
		molecule = ''
		lines = sdf_file.readlines()
		for mol_line in lines:
			molecule += mol_line
			if mol_line == "$$$$\n" or mol_line == '$$$$\r\n':
				molecule_list.append(molecule)
				molecule = ''
	return molecule_list


def calc_param(file):
	command_1 = 'cxcalc -i name -S pkacalculator -t basic -b 1 logD -H 7.4 logP donorcount acceptorcount polarsurfacearea mass aromaticringcount rotatablebondcount atomcount atomcount -z 1 '
	if str(file).startswith('/'):
		command_1 += str(file)
	else:
		command_1 += os.getcwd() + '/' + str(file)
	sdf_file = subprocess.check_output([command_1], shell=True)
	sdf_file = sdf_file.decode()
	return sdf_file


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-i', '--in_file', help='Name of the input file [MENDATORY]', metavar='')
	parser.add_argument('-S', '--sfi',
	                    help='Use this flag if you want to report SFI index (solubility Forecast index) the field '
	                         '"ArRings" (count of '
	                         'aromatics rings) becomes required in the sdf fields if the -C option is used (by '
	                         'default option is off) ', action='store_true', default=False)
	parser.add_argument('-C', '--no_calc',
	                    help='Use this flag to disable automatic calculation of chemical properties (need to be '
	                         'provided in the sdf file then (beware the sdf fields need to have the same name as '
	                         'following: bpKa (basic pKa), logD (logD at pH=7.4), logP, HBD (Hydrogen bond donor '
	                         'count), TPSA (total polar surface area)'
	                         ', MW (molecular weight)', action='store_true', default=False)
	parser.add_argument('-f', '--output_folder', type=str, help='create a folder to output the file', metavar='')
	parser.add_argument('-v', '--verbose', help='more details output in the terminal', action='store_true',
	                    default=False)
	args = parser.parse_args()
	if not args.in_file:
		parser.print_help()
		sys.exit()
	while True:
		if args.in_file:
			if os.path.exists(args.in_file):
				break
			else:
				print('ERROR : file inputed as argument [-i] does not exist')
				sys.exit()
	if not args.no_calc:
		sdf_withparam = calc_param(args.in_file)
		mol_list = calc_sdf_parser(sdf_withparam)
	else:
		mol_list = sdf_parser(args.in_file)
	fpathname = args.in_file.split('/')  # parsing the file name to use it for the output
	fpath = '/'.join(fpathname[:-1]) + '/'
	if fpath == '/':
		fpath = ''
	fname = fpathname[-1].split('.')[0]
	mol_to_out = ''
	for m in mol_list:  # scanning through the molecules
		m = Desc(m)  # initializing the Desc object with all the parameters
		mol_to_out = mol_to_out + m.sdf_writer()
		if args.verbose:
			m.print_details(
				1)  # Print the details in the terminal. By default level = 0, enter Print_details(1) for more details
		else:
			m.print_details(0)
	if args.output_folder:
		if str(args.output_folder).startswith('/'):
			path_to_check = args.output_folder
		else:
			path_to_check = fpath + str(args.output_folder)
		if not os.path.exists(path_to_check):
			os.makedirs(path_to_check)
		with open(path_to_check + '/' + fname + "_out.sdf", 'w') as output:
			output.write(mol_to_out)
	else:
		with open(fpath + fname + "_out.sdf", 'w') as output:
			output.write(mol_to_out)
