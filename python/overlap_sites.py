#### over lap for sanple sites
#### usage : python overlap_sites.py <sample 1> <sample 2> ...<sample n> > overlap.txt
### modify parameter : value  to get minimal overlap times
import sys
import collections
from collections import Counter

value = 2
dict_c = collections.defaultdict(list)
c = collections.defaultdict(lambda:0)

def get_AG(site_file):
	c = collections.defaultdict(list)
	with open(site_file) as f:
		for line in f:
			line = line.strip().split('\t')
			if len(line) ==1:
				continue
			chrID = str(line[0])
			pos = str(line[1])
			ref_base = str(line[2])
			alt_base = str(line[3])
			coverage = str(line[4])
			strand = str(line[5])
			var_num = str(line[6])

			name = chrID+' '+pos+' '+ref_base+' '+alt_base+' '+strand
			content = coverage + ' '+var_num
			c[name].append(content)

	return c

for i in sys.argv:
	c = get_AG(i)
	for k,v in c.items():
		dict_c[k].append(v)

for k,v in dict_c.items():
	if len(v) >= value:
		temp_li = []
		for res in v:
			temp_li.extend(res)
		print(str(k)+"\t"+"\t".join(temp_li))
