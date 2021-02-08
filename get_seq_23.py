#!/usr/bin/env python
import sys
import pysam
import numpy as np
import collections
from collections import Counter
values = 12

def get_site_arround(sites_path, fasta_path, values):

	fasta = pysam.FastaFile(fasta_path)
	sequences = []
	with open(sites_path) as f:
		for line in f:
			line = line.strip().split('\t')
			chr_id = line[0]
			start = int(line[1]) - values
			end = int(line[1]) + values -1
			seq = fasta.fetch(chr_id, start, end).upper()
			sequences.append(seq)

	return sequences

sites_path = sys.argv[1]
fasta_path = sys.argv[2]

sequences = get_site_arround(sites_path, fasta_path, values)

BASES = "A T C G N".split()
d = collections.defaultdict(list)
c = collections.defaultdict(lambda :0)
for i in sequences:
	for j in range(len(i)):
		base = i[j]
		c[base] +=1
		d[j].append(base)


#for k, v in d.items():
#	Counter(v)
#	print(k,Counter(v))
count_A,count_T,count_C,count_G = 0,0,0,0

for k,v in c.items():
	if k == 'A':
		count_A = v/23
	elif k == 'T':
		count_T = v/23
	elif k == 'C':
		count_C = v/23
	elif k == 'G':
		count_G = v/23
print("A\tT\tC\tG")

for k, v in d.items():
	num = {}
	num = Counter(v)
	for i, j in num.items():
		if i =='A':
			base_A = j/count_A
		elif i =='T':
			base_T = j/count_T
		elif i =='C':
			base_C = j/count_C
		elif i =='G':
			base_G = j/count_G
	row = str(base_A) + '\t' +str(base_T) + '\t' +str(base_C) + '\t' +str(base_G)
	print(row)
