#!/usr/bin/env python
import sys
import pysam
import numpy as np
import collections
from collections import Counter
values = 5

def complement(s):
	basecomplement = {
         "A":"T",
          "T":"A",
          "G":"C",
          "C":"G",
          "a":"t",
          "t":"a",
          "g":"c",
          "c":"g",
			"N":"N",}
	letters = list(s)
	letters = [basecomplement[base] for base in letters]
	letters = letters[::-1]

	return ''.join(letters)

def get_site_arround(sites_path, fasta_path, values):

	fasta = pysam.FastaFile(fasta_path)
	sequences = []
	with open(sites_path) as f:
		for line in f:
			line = line.strip().split('\t')
			chr_id = line[0]
			strand = line[5]
			start = int(line[1]) - values -1
			end = int(line[1]) + values
			seq = fasta.fetch(chr_id, start, end).upper()
			if 'neg' in strand:
				seq = complement(seq)
			#	print(seq)
			else:
				seq = seq
			sequences.append(seq)

	return sequences

sites_path = sys.argv[1]
fasta_path = sys.argv[2]

sequences = get_site_arround(sites_path, fasta_path, values)
s=0
for i in sequences:
	s+=1
	print(">"+str(s))
	print(i)
#for k,v in Counter(sequences).items():
#	print(k,v)
#d = collections.defaultdict(list)
#c = collections.defaultdict(lambda :0)
#for i in sequences:
#	c = collections.defaultdict(lambda :0)
#	for j in range(len(i)):
#		base = i[j]
#		c[base] +=1
#		d[j].append(base)
#	print(c)

