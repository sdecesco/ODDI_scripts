#!/usr/bin/env python
import argparse
import os,sys


def count_mol(sdf_file):
	with open(sdf_file, 'r') as sdf:
		count = 0
		for line in sdf:
			if line.startswith('$$$$'):
				count += 1
	return count


def get_split_counts(sdf_file, n_split):
	number_of_mol = count_mol(sdf_file)
	n_mol, rest = divmod(number_of_mol, n_split)
	n_per_file = [n_mol] * n_split
	n_per_file[-1] = n_per_file[-1] + rest
	return n_per_file


def split_file(file_name, n_files=2):
	filename, file_extension = os.path.splitext(file_name)
	n_split = get_split_counts(file_name, n_split=n_files)

	with open(file_name, 'r') as sdf_in:
		f_count = 0
		mol_count = 0
		output_name = filename + '_split_' + str(f_count + 1) + file_extension
		out = open(output_name, 'w')
		for line in sdf_in:
			if line.startswith('$$$$'):
				mol_count += 1
			if mol_count == n_split[f_count]:
				mol_count = 0
				f_count += 1
				out.write(line)
				out.close()
				if f_count < len(n_split):
					output_name = filename + '_split_' + str(f_count + 1) + file_extension
					out = open(output_name, 'w')
				else:
					break
			else:
				out.write(line)


if __name__ == "__main__":

	parser = argparse.ArgumentParser(description="Split your SDF file into N files")
	parser.add_argument('-i', '--in_file', help='Name of the input file [MENDATORY]', metavar='')
	parser.add_argument('-n', '--number_of_files',type=int,default=2,help='Number of files to split your SDF file', metavar='')
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

	split_file(args.in_file,n_files=args.number_of_files)