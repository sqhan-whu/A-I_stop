# Date: 20200301
# author: Han Shaoqing
# usages: analysis for weiqi A-I stop site( Ni Nt)
#####################################################################

from os.path import basename
from sys import argv

def get_RT_dic(samf):
	with open(samf) as f:
		forward_reads = {}
		reverse_reads = {}
		f_pos_coverage = {}
		r_pos_coverage = {}
		freads_num = 0
		rreads_num = 0
		for line in f:
			if not line.strip().startswith("@"):
				flag = line.split('\t')[1]
				ch = line.split('\t')[2]
				pos = line.split('\t')[3]
				seq = line.split('\t')[9]
				tlen= len(seq)
				if (str(flag) == '0'):
					position = "+\t"+str(ch)+"\t"+str(int(pos))
					forward_reads[position] = forward_reads.get(position,0) + 1
					freads_num +=1
					for l in range(tlen+1):
						position = "+\t"+str(ch)+"\t"+str(int(pos)+l)
						f_pos_coverage[position] = f_pos_coverage.get(position,0) + 1
				elif (str(flag) == '16'):
					position = "-\t"+str(ch)+"\t"+str(int(pos))
					reverse_reads[position] = reverse_reads.get(position,0) + 1
					rreads_num +=1
					for l in range(tlen+1):
						position = "-\t"+str(ch)+"\t"+str(int(pos)-l)
						r_pos_coverage[position] = r_pos_coverage.get(position,0) + 1
				else:
					next
		for k,v in reverse_reads.items():
			if k in forward_reads:
				forward_reads[k] += v
			else:
				forward_reads[k] = v
	
	f.close()

	return forward_reads, freads_num, rreads_num, freads_num+rreads_num, f_pos_coverage, r_pos_coverage

forward_reads, freads_num, rreads_num, totalreads_num, f_pos_coverage, r_pos_coverage = get_RT_dic(argv[1])

print("dirction\tchrID\tposition\tNi\tNt")

for k,v in forward_reads.items():
	dirction, chrID, pos = k.split('\t')
	if dirction == '+':
		cov = f_pos_coverage[k]
		print(str(k)+"\t"+str(v)+"\t"+str(cov))
	else:
		cov = r_pos_coverage[k]
		print(str(k)+"\t"+str(v)+"\t"+str(cov))
#if __name__ == '__main__':
#	main()
