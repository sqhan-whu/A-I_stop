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
			elif float(line[9]) > 26 and line[3] == 'A':
				row = line[3]+line[4]
				c[row] +=1
			else:
				continue
	
	for k,v in c.items():
		print(k,v,site_file)
for i in sys.argv:
	get_AG(i)
