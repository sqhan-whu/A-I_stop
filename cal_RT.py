from os.path import basename
import sys
import collections

def cal_RT(mutation_file):

	d = collections.defaultdict(lambda :0)
	#x = []
	with open(mutation_file) as RT:
		for line in RT:
			line = line.strip('\n').split('\t')
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

				if (int(read_length) - int(relative_pos)) == 1 and int(base_qual) > 26:
					content = str(chrID)+'\t'+str(position)+'\t'+ str(ref_base)+'\t'+ str(alt_base)+ '\t'+ str(coverage)
			#	x.append(content)

					d[content] +=1
				else:
					continue
			#x = set(x)

		for k, v in d.items():
			print(str(k)+"\t"+str(v))

mutation_file  = sys.argv[1]

cal_RT(mutation_file)
