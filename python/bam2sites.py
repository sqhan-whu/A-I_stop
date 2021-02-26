import sys, os

hg38 = "/home/weiqi/project/00.DATABASE/hg38/hg38.fa"
hg38_size = "/home/weiqi/project/00.DATABASE/hg38/hg38.fa.size"
mm10 = "/home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa"
mm10_size = "/home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa.size"

def process_bam2sites(pos_bam, neg_bam):
	name = pos_bam.split('.')[0]
	call_pos_sites = 'python ~/bin_python/I_script/AG.py ' + pos_bam + ' ' + mm10 + ' '+mm10_size+' '+pos_bam+'.AG.bam > '+pos_bam+'.pos.rough.txt'
	call_neg_sites = 'python ~/bin_python/I_script/TC.py ' + neg_bam + ' ' + mm10 + ' '+mm10_size+' '+neg_bam+'.TC.bam > '+neg_bam+'.neg.rough.txt'
	
	trans_neg2pos = 'python ~/bin_python/I_script/trans_neg2pos.py '+ pos_bam+'.pos.rough.txt '+neg_bam+'.neg.rough.txt > '+ name+'.process.site.txt'
	stat_pos_site = 'python ~/bin_python/I_script/stat_pos_site.py '+name+'.process.site.txt > '+ name+'.site.txt'

#	print(call_pos_sites)
	os.system(call_pos_sites)
	os.system(call_neg_sites)
	os.system(trans_neg2pos)
	os.system(stat_pos_site)

pos_bam = sys.argv[1]
neg_bam = sys.argv[2]

process_bam2sites(pos_bam, neg_bam)
