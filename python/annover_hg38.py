####usage: python script sample_sites.txt 

import sys, os
import collections

def get_annvar_input(sites_file):
	c = collections.defaultdict(list)
	with open(sites_file) as f, \
	open('tmp_'+sites_file,'w') as out:
		for line in f:
			line = line.strip().split('\t')
			chr_id = line[0]
			#chr_id.replace("chr", "")
			pos = line[1]
			ref_base = line[2]
			alt_base = line[3]
			row = chr_id + '\t' + pos+'\t'+pos+'\t'+ref_base+'\t'+alt_base+'\n'
			out.write(row)


def annvar_process(sites_file):
	input_file = 'tmp_'+sites_file

	command1 = 'perl ~/project/software/opt/annovar/annotate_variation.pl -filter -out filter_snp_'+sites_file+' -build hg38 -dbtype snp151Common '+'tmp_'+sites_file+' ~/project/software/opt/annovar/hg38db'
	command2 = 'perl ~/project/software/opt/annovar/annotate_variation.pl -out function_'+sites_file+' -dbtype ensGene -build hg38 '+'filter_snp_'+sites_file+'.hg38_snp151Common_filtered ~/project/software/opt/annovar/hg38db'
	
	remove = 'rm '+ input_file
	os.system(command1)
	os.system(command2)
	os.system(remove)


def get_raw_information(sites_file):
	input_file = 'filter_snp_'+sites_file+'.hg38_snp151Common_filtered'

	command3 = 'cut -f 1-2 '+input_file+' > tmp'+sites_file
#	command4 = 'cat '+sites_file+' |grep -w -f tmp'+sites_file+' > final_'+sites_file
	command4 ='awk \'NR==FNR{s[$1"\t"$2]=$0;next}NR>FNR{print s[$1"\t"$2]}\' '+sites_file+' tmp'+sites_file+' >final_'+sites_file
	command5 = 'rm tmp'+sites_file
	os.system(command3)
	os.system(command4)
	os.system(command5)

sites_file = sys.argv[1]
get_annvar_input(sites_file)
annvar_process(sites_file)
get_raw_information(sites_file)
