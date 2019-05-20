#!/usr/bin/env python
import sys
file_count = 1
filename = 'chembl_22_1.sdf'
file_out = 'chembl_22_1_'+str(file_count)+'.sdf'
max_files = 2
number_of_cpds_per_file = 100

with open(filename,'r') as file:
	count=0
	out = open(file_out,'w')
	for line in file:
		if line == '$$$$\n':
			count+=1
			print 'new mol: '+str(count)
		if count == number_of_cpds_per_file:	
			count = 0
			file_count+=1
			out.write(line)
			out.close()
			file_out = 'MIDAS_SDF_EXPORT2_'+str(file_count)+'.sdf'
			out = open(file_out,'w')
			print 'new file created'
			continue
		else:
			out.write(line)
		if file_count==max_files:
			out.close()
			sys.exit()
out.close()