import sys

with open(sys.argv[1]) as f:
	for line in f:
		line = line.strip().split('\t')
		if line[2] == 'T':
			ref_base = 'A'
		elif line[2] == 'A':
			ref_base = 'T'
		elif line[2] == 'C':
			ref_base = 'G'
		elif line[2] == 'G':
			ref_base = 'C'
		if line[3] == 'T':
			alt_base = 'A'
		elif line[3] == 'A':
			alt_base = 'T'
		elif line[3] == 'C':
			alt_base = 'G'
		elif line[3] == 'G':
			alt_base = 'C'
		row = line[0] +'\t'+ line[1]+ '\t'+ref_base +'\t'+ alt_base +'\t'+line[4] +'\t'+line[5]
		print(row)
