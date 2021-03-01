import sys
import collections
from collections import Counter

def get_AG(site_file):
	c = collections.defaultdict(lambda :0)
	with open(site_file) as f:
		for line in f:
			line = line.strip().split('\t')
			if len(line) == 1:
				continue
			else:
				chrID = line[0]
				position = line[2]
				ref_base = line[3]
				alt_base = line[4]
				coverage = line[5]
				relative_pos = line[7]
				read_length = line[8]
				base_qual = line[9]
				read_name = line[11]
				strand = line[13]

				if int(line[8]) - int(line[7]) == 1 and float(line[9]) > 26 and line[3] == 'A' and line[4] == 'G':
					content = str(chrID)+'\t'+str(position)+'\t'+ str(ref_base)+'\t'+ str(alt_base)+ '\t'+ str(coverage)+'\t'+strand
					c[content] +=1
				else:
					continue
	
	for k,v in c.items():
		if v >=4 and float(v)/float(k.split('\t')[4]) > 0.1:
			print(k,"\t",v)
for i in sys.argv:
	get_AG(i)
