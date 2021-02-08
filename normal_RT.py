import sys
import collections
from scipy import stats
import math
import os
#factor = 1

def cal_RT_value(mutations, factor):

	dict_file = collections.defaultdict(lambda :0)

	with open(mutations) as f:
		for line in f:
			line = line.strip('\n').split('\t')
			if len(line) == 1:
				continue
			else:
				chrID = line[0]
				position = line[1]
				ref_base = line[2]
				alt_base = line[3]
				coverage = line[4]
				mutation = int(line[5])

				RT_value = math.log(mutation+1,2) / factor

				k = str(chrID)+'\t'+str(position)
				content = str(chrID)+'\t'+str(position)+'\t'+ str(ref_base)+'\t'+ str(alt_base)+ '\t'+ str(coverage) + '\t' + str(RT_value)

				dict_file[k] = content

	return dict_file

def compare_RT(rt1, rt2):

	value1 = []
	value2 = []
	base = []
	row = []

	rt1_tx_ls = [i for i in rt1]
	rt2_tx_ls = [i for i in rt2]

	rt1_rt2_common_tx_ls = set(rt1_tx_ls) & set(rt2_tx_ls)

	for i in rt1_rt2_common_tx_ls:
		
		value1.append(float(rt1[i].split('\t')[5]))
		value2.append(float(rt2[i].split('\t')[5]))
		base.append(rt1[i].split('\t')[2])
		new_row = str(i) + '\t' + str(rt1[i].split('\t')[5])+'\t'+str(rt2[i].split('\t')[5])+'\t'+ str(rt1[i].split('\t')[2])
		row.append(new_row)
	r,p = stats.pearsonr(x=value1,y=value2)

	return row, r, p

def cal_factor(bam_file_1,bam_file_2):
	size1 = 'bam_file_1.size'
	size2 = 'bam_file_2.size'

	samtools_command = 'samtools view -c '+bam_file_1 +' > '+ size1
	samtools_command2 = 'samtools view -c '+bam_file_2 +' > '+ size2
	os.system(samtools_command)
	os.system(samtools_command2)
	
	with open(size1) as s1, \
		open(size2) as s2:
			factor1 = float(s1.readline().strip('\n')) / 1000000
			factor2 = float(s2.readline().strip('\n')) / 1000000
	os.system('rm '+size1)
	os.system('rm '+size2)
	return factor1, factor2

bam_file_1=sys.argv[3]
bam_file_2=sys.argv[4]

factor1, factor2 = cal_factor(bam_file_1,bam_file_2)

rt1 = cal_RT_value(sys.argv[1], factor1)
rt2 = cal_RT_value(sys.argv[2], factor2)

row, r, p = compare_RT(rt1,rt2)
for i in row:
	new_row = ''
	i = i.split('\t')
	new_row = str(i[0])+'\t'+str(i[1])+'\t'+str(i[2])+'\t'+str(i[3])+'\t'+str(i[4])
print(new_row)
print(r, p)
