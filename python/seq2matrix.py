import sys
import collections
from collections import Counter
d = collections.defaultdict(list)

def get_matrix(seq_fa):
	with open(seq_fa) as f:
		for line in f:
			line = line.strip()
			if line.startswith('>'):
				continue
			else:
				for i in range(len(line)):
					if i ==5:
						base = 'I'
					else:
						base = line[i]
					d[i].append(base)
#					print(line[i])
	return d

d = get_matrix(sys.argv[1])
#print("A U C G I")
for k,v in d.items():
	#print('A'+' '+str(Counter(v)['A'])+' '+'U'+' '+str(Counter(v)['T'])+' '+'C'+' '+str(Counter(v)['C'])+' '+'G'+' '+str(Counter(v)['G'])+' '+'I'+' '+str(Counter(v)['I']))
	print(str(Counter(v)['A'])+' '+str(Counter(v)['T'])+' '+str(Counter(v)['C'])+' '+str(Counter(v)['G'])+' '+str(Counter(v)['I']))
