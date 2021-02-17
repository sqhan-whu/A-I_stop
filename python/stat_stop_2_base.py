## usage: python stat_stop_2_base.py 1.rough.txt 1neg.rough.txt 
## get stop read position -2 base (qual > 26)


import sys
from collections import Counter

def get_count(pos_file,neg_file):
	base = []
	with open(pos_file) as p:
		for line in p:
			line = line.strip().split('\t')
			if len(line) == 1:
				continue
			else:
				alt_base = line[4]
				abs_pos = int(line[7])
				read_len = int(line[8])
				qual = int(line[9])
				if read_len - abs_pos == 1 and qual > 26:
					base.append(alt_base)
				else:
					continue

	with open(neg_file) as n:
		for line in n:
			line = line.strip().split('\t')
			if len(line) == 1:
				continue
			else:
				alt_base = line[4]
				abs_pos = int(line[7])
				read_len = int(line[8])
				qual = int(line[9])
				pos_alt_base = ''
				if abs_pos == 2 and qual > 26:
					if alt_base == 'A':
						pos_alt_base = 'T'
					elif alt_base == 'T':
						pos_alt_base = 'A'
					elif alt_base == 'C':
						pos_alt_base = 'G'
					elif alt_base == 'G':
						pos_alt_base = 'C'
						
					base.append(pos_alt_base)
				else:
					continue

	return base

pos_file = sys.argv[1]
neg_file = sys.argv[2]

base = get_count(pos_file,neg_file)
base_number = Counter(base)

print('A:'+'\t'+str(base_number['A']))
print('T:'+'\t'+str(base_number['T']))
print('C:'+'\t'+str(base_number['C']))
print('G:'+'\t'+str(base_number['G']))
