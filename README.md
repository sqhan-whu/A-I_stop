
# Slic-seq
----------------------------------------
## Scripts for analysing Slic-seq data ##
----------------------------------------
Strand-specific Slic-seq data analysis process includes the basic processing tools for quality control, clean reads mapping, and reads strand-specific resolution statistics of truncated reads and so on.

Slic-seq is a novel and effective biochemical method for transcriptome-wide identification of inosine based on Endonuclease V cleavage activity, which achieved the specific ligation of inosine-cleaved sequencing. This robust and straightforward approach substantially enhances the detection and scope of A-to-I editing sites in cellular RNA while achieving the enrichment of inosine-containing RNAs.

----------------------------------------
### The link address:
Github: https://github.com/sqhan-whu/A-I_stop/
-----------------------------------------

### Data analysis process
------------------------------------	
### Quick start

**Download test data**

Users can download example Scli-seq data in Hela cells from Gene Expression Omnibus (GEO):
```

wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14082811/SRR14082811 ./ &
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14082812/SRR14082812 ./ &
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14082813/SRR14082813 ./ &
wget https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14082814/SRR14082814 ./ &

mv SRR14082811 Hela_Scli-seq.rep1.sra
mv SRR14082812 Hela_Scli-seq.rep2.sra
mv SRR14082813 Hela_control.rep1.sra
mv SRR14082814 Hela_control.rep2.sra
        
```	

### In the process of mapping steps, special attention should be paid to the parameter setting of hisat2 software, "-- RNA-strandness RF --no-softclip". Since data mutations will exist at ends of the reads, all mutations should be retained during comparison without softclip.
 
**Adapter seqeucnes were trimmed by using trim_galore. fastuniq is used to remove PCR duplication. Then use seqtk trimfq '-e 10' or '-b 6 -e 10' remove random sequences added in library construction. Clean reads were mapped to genome by using HISAT2; with parameter "--mp 4,2 --rna-strandness RF --no-softclip --no-mixed --no-discordant". A custom python script was used to parse the pileup format into a tabular format summarizing the mutation at each position. Genomes: hg38 and mm10 (https://www.gencodegenes.org/). Example: 00_mapping.pl **


** Filter poor quality sequence**

```
#!/bin/bash
bamtools filter -in Hela_Scli-seq.rep1.bam -out Hela_Scli-seq.rep1.filter.bam -isPrimaryAlignment true -isProperPair true -isMapped true -mapQuality '>1'
bamtools filter -in Hela_Scli-seq.rep2.bam -out Hela_Scli-seq.rep2.filter.bam -isPrimaryAlignment true -isProperPair true -isMapped true -mapQuality '>1'
bamtools filter -in Hela_control.rep1.bam -out Hela_control.rep1.filter.bam -isPrimaryAlignment true -isProperPair true -isMapped true -mapQuality '>1'
bamtools filter -in Hela_control.rep2.bam -out Hela_control.rep2.filter.bam -isPrimaryAlignment true -isProperPair true -isMapped true -mapQuality '>1'
```
	

**slic-seq peaks calling (optional)**	
         
Call slic-seq peaks:

```
#!/bin/bash	 
python call.rpmf.peak.py -i Hela_Scli-seq.rep1 -c Hela_control.rep1 -w 3
python call.rpmf.peak.py -i Hela_Scli-seq.rep2 -c Hela_control.rep2 -w 3
python peak_finalize.py -p Hela_Scli-seq.rep1_RPMF_out.txt -a gencode.v31.gene.slidingwindows.gtf
python peak_finalize.py -p Hela_Scli-seq.rep2_RPMF_out.txt -a gencode.v31.gene.slidingwindows.gtf

bedtools  intersect -a Hela_Scli-seq.rep1.filter.bam  -b Hela_Scli-seq.rep1_RPMF_out_final.bed  > Hela_Scli-seq.rep1.filter.region.bam
bedtools  intersect -a Hela_Scli-seq.rep2.filter.bam  -b Hela_Scli-seq.rep2_RPMF_out_final.bed  > Hela_Scli-seq.rep2.filter.region.bam
```


**Reads strand-specific resolution of slic-seq data**

```
perl $0 Hela_Scli-seq.rep1.bam 

#!/usr/bin/env perl
use Cwd qw(abs_path);
use File::Basename;
$bam = shift;
$BIN='I_script';
$CUTOFF = 50;
my ($fname,$path,$ext)=fileparse($bam,'.filter.bam');
$path =~ s/\/\s*$//;
open F, ">$fname\_cal_snp.sh";
print F "samtools view $bam |grep 'XS:A:+' > $path/$fname.pos.sam\n";
print F "samtools view $bam |grep 'XS:A:-' > $path/$fname.neg.sam\n";
print F "samtools view -H $bam > $path/$fname.header.h\n";
print F "cat $path/$fname.header.h $path/$fname.pos.sam > $path/$fname.pos.0.sam\n";
print F "cat $path/$fname.header.h $path/$fname.neg.sam > $path/$fname.neg.0.sam\n";

print F "samtools view -b -S $path/$fname.pos.0.sam > $path/$fname.pos.bam\n";
print F "samtools view -b -S $path/$fname.neg.0.sam > $path/$fname.neg.bam\n";
print F "bamtools sort -in $path/$fname.pos.bam -out $path/$fname.pos.sort.bam\n";
print F "bamtools sort -in $path/$fname.neg.bam -out $path/$fname.neg.sort.bam\n";
print F "samtools index $path/$fname.pos.sort.bam\n";
print F "samtools index $path/$fname.neg.sort.bam\n";

```

Using strand-specific bam files potential search A-I sites (RPKM):

```
python bam2sites_hg38.py Hela_Scli-seq.rep1.pos.sort.bam Hela_Scli-seq.rep1.neg.sort.bam
python bam2sites_hg38.py Hela_Scli-seq.rep2.pos.sort.bam Hela_Scli-seq.rep2.neg.sort.bam


Resluts : Hela_Scli-seq.rep1_process_sites_files.txt
```


# Run from "process_sites_files".

```
example: python $0 Hela_Scli-seq.rep1_process_sites_files.txt

import sys
import collections
from collections import Counter

def Merge(dict1, dict2): 
	res = dict()
	for k,v in dict1.items():
		if k in dict2:
			v = dict2[k] + v
			res[k] = v 
		else:
			v =  v 
			res[k] = v
	return res


def reverse_base(base):
	if base == 'A':
		rev_base = 'T'
	elif base == 'T':
		rev_base = 'A'
	elif base == 'C':
		rev_base = 'G'
	elif base == 'G':
		rev_base = 'C'

	elif base == 'N':
		rev_base = 'N'

	return rev_base

def get_AG(site_file):
	dict_I_low = dict()
	cal_3bp_sites = dict()
	filter_sites = dict()
	dict_rawmutation = dict()
	dict_I = dict()
	c = collections.defaultdict(lambda :0)
	d = collections.defaultdict(lambda :0)
	with open(site_file) as f:
		for line in f:
			line = line.strip().split('\t')
			if len(line) == 1:
				continue
			else:
				chrID = line[0]
				position = line[2]
				ref_base = line[3]
				alt_base = line[4]
				coverage = float(line[5])
				relative_pos = line[7]
				read_length = line[8]
				base_qual = line[9]
				read_name = line[11]
				strand = line[13]

				if int(line[8]) - int(line[7]) == 1 and coverage > 0 and float(line[9]) > 26:# and float(line[9]) > 10:# and line[3] == 'A' and line[4] == 'G':

					content = str(chrID)+'\t'+str(position)+'\t'+ str(ref_base)+'\t'+ str(alt_base)+ '\t'+ str(coverage)+'\t'+strand

					c[content] +=1
				if coverage > 0 and float(line[9]) > 26:

					content = str(chrID)+'\t'+str(position)+'\t'+ str(ref_base)+'\t'+ str(alt_base)+ '\t'+ str(coverage)+'\t'+strand
					d[content] +=1
						
#	dict3 = Merge(c, d)
	
#	print(dict3)
#	print(c)
#	print(d)
	for k,v in c.items():
		ref_base = k.split('\t')[2]
		alt_base = k.split('\t')[3]
		strand = k.split('\t')[5]
		position = k.split('\t')[1]
		chrID = k.split('\t')[0]
		coverage = k.split('\t')[4]
		new_content= str(chrID)+'\t'+str(position)
		new_value = str(ref_base)+'\t'+ str(alt_base)+ '\t' +str(coverage)+'\t'+strand+'\t'+ str(v)
		if v >= 6 and float(v)/float(k.split('\t')[4]) >= 0.1 and float(coverage) < 5000  and 'M' not in chrID:
			dict_I[new_content] = new_value
		if v <6 and v >=3 and float(v)/float(k.split('\t')[4]) >= 0.05 and float(coverage) < 5000 and 'M' not in chrID:
			dict_I_low[new_content] = new_value

	for k,v in d.items():
		ref_base = k.split('\t')[2]
		alt_base = k.split('\t')[3]
		strand = k.split('\t')[5]
		position = k.split('\t')[1]
		chrID = k.split('\t')[0]
		coverage = k.split('\t')[4]
		new_content= str(chrID)+'\t'+str(position)
		new_value = str(ref_base)+'\t'+ str(alt_base)+ '\t' +str(coverage)+'\t'+strand+'\t'+ str(v)
		if v >= 6 and float(v)/float(k.split('\t')[4]) >= 0.1:
			dict_rawmutation[new_content] = new_value
	
	for m,n in dict_I.items():
		chrID = m.split('\t')[0]
		position = m.split('\t')[1]
		p = 0
		for i in range(1,3):
			l1_position = int(position) - i
			r1_position = int(position) + i
			wrong_l1 = chrID + "\t" + str(l1_position)
			wrong_r1 = chrID + "\t" + str(r1_position)
			right_partten = n.split('\t')[0] + "\t" + n.split('\t')[1]
		
			if wrong_l1 in dict_rawmutation and float(dict_rawmutation[wrong_l1].split("\t")[4]) / float(dict_rawmutation[wrong_l1].split("\t")[2]) > 0.2:
				wrong_l1_partten =  dict_rawmutation[wrong_l1].split("\t")[0]+"\t"+dict_rawmutation[wrong_l1].split("\t")[1]
				if  wrong_l1_partten != right_partten or wrong_l1_partten == "A\tA":
					p = 1
			#	else:
			#		filter_sites[m]=n

			if wrong_r1 in dict_rawmutation and float(dict_rawmutation[wrong_r1].split("\t")[4]) / float(dict_rawmutation[wrong_r1].split("\t")[2]) > 0.2:
				wrong_r1_partten =  dict_rawmutation[wrong_r1].split("\t")[0]+"\t"+dict_rawmutation[wrong_r1].split("\t")[1]
				if	wrong_r1_partten != right_partten or wrong_r1_partten == "A\tA":
					p = 1
			
		if p == 0:
			filter_sites[m]=n
	
	for k,v in filter_sites.items():
		chrID = k.split('\t')[0]
		position = k.split('\t')[1]
	#	print(chrID+"\t"+position)
		for i in range(1,3):
			l1_position = int(position) - i
			r1_position = int(position) + i
			wrong_l1 = chrID + "\t" + str(l1_position)
			wrong_r1 = chrID + "\t" + str(r1_position)
			if wrong_l1 in dict_I_low:
				cal_3bp_sites[wrong_l1] = dict_I_low[wrong_l1]
			if wrong_r1 in dict_I_low:
				cal_3bp_sites[wrong_r1] = dict_I_low[wrong_r1]
	with open(site_file+"_high_sites.txt",'w') as o1:
		for k,v in filter_sites.items():
			o1.write(k+"\t"+v+"\n")
	with open(site_file+"_low_sites.txt",'w') as o2:
		for k,v in cal_3bp_sites.items():
			o2.write(k+"\t"+v+"\n")

	with open(site_file+"_total.sites.txt",'w') as o3:
		for k,v in filter_sites.items():
			o3.write(k+"\t"+v+"\n")
		for k,v in cal_3bp_sites.items():
			o3.write(k+"\t"+v+"\n")

	o1.close()
	o2.close()
	o3.close()
#	print(dict_I_low)
#print("chrom\tposition\tref_base\talt_base\tcoverage\tstrand\ttrunc_reads")
for i in sys.argv[1:]:
	get_AG(i)
```

**Results: $sample.sites_out.txt**
```
chrom	position	ref_base	alt_base	coverage	strand	trunc_reads
chr1	1007403	A	G	16.0	pos	3
chr1	1009425	A	G	14.0	neg	4
chr1	1279589	A	G	40.0	pos	4
chr1	1318496	A	G	39.0	pos	4
chr1	1318492	A	G	45.0	pos	4
chr1	1358961	A	G	26.0	pos	3
chr1	1403078	A	G	22.0	pos	4
chr1	1409620	A	G	78.0	pos	4
chr1	1445880	A	G	50.0	pos	4
chr1	1446279	A	G	74.0	pos	5
```


# Computing site motif.
```
## Get the sequence near potential A-I sites.
Usage: python $0 $sample.sites_out.txt genome.fa

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
#		f.readline()
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
```

### Run for "$sample.sites.fa".
```
python $0 $sample.sites.fa &


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
	return d

d = get_matrix(sys.argv[1])
#print("A U C G I")
for k,v in d.items():
	print(str(Counter(v)['A'])+' '+str(Counter(v)['T'])+' '+str(Counter(v)['C'])+' '+str(Counter(v)['G'])+' '+str(Counter(v)['I']))

```

# Check the proportion of mutations in the data.
The proportion of mutations at truncated sites is an important indicator of data reliability.
```
Usage: python $0 $sample.process_sites_files.txt

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
#			elif int(line[8]) - int(line[7]) == 1 and float(line[9]) > 26 and line[3] == 'A':
			elif int(line[8]) - int(line[7]) == 1 and float(line[9]) > 26:
				row = line[3]+line[4]
				c[row] +=1
			else:
				continue
	
	for k,v in c.items():
		print(k,v,site_file)
for i in sys.argv:
	get_AG(i)
```


# Annotation for A-I sites by annovar software.

```
Usages: python $0 $sample.sites_out.txt
	
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
	command4 ='awk \'NR==FNR{s[$1"\t"$2]=$0;next}NR>FNR{print s[$1"\t"$2]}\' '+sites_file+' tmp'+sites_file+' |sort -k1,1 -k2,2n >final_'+sites_file
	command5 = 'rm tmp'+sites_file
	os.system(command3)
	os.system(command4)
	os.system(command5)

sites_file = sys.argv[1]
get_annvar_input(sites_file)
annvar_process(sites_file)
get_raw_information(sites_file)
```


# Differential analysis of inosine sites.
Statistical analysis of the A-to-I editing sites in each mouse disease model was done. The calculation of differential A-to-I editing was based on the difference in the number of truncated reads enriched in different mouse disease model samples and control samples. Difference statistics are done by using the edgeR package (FDR < 0.05 and log2FC > 1).

```
Preparing Input Data as follows :

id	Alzheimer_rep1 Alzheimer_rep2	Alzheimer_control1	Alzheimer_control2
chr1_11601083_UTR3_pos	2	6	5	12
chr1_13625015_UTR3_neg	22	19	31	31
chr1_21961561_exonic_neg	14	29	6	1
chr1_22704184_intronic_neg	4	2	16	10
chr1_22705259_intronic_neg	4	2	16	8
chr1_22705274_intronic_neg	2	4	36	10
chr1_24051529_intronic_neg	2	14	18	22
chr1_52727291_UTR5_neg	12	14	4	2
chr1_59179316_intronic_neg	6	6	4	10
chr1_60478588_UTR3_pos	4	18	4	6
chr1_66468370_UTR5_pos	4	21	6	6
chr1_66672714_exonic_pos	123	364	49	43
...
```
The calculation of differential A-to-I editing was based on the difference in the number of truncated reads, differences can be calculated using R package, such as edgeR.

```
#Get raw data Alzheimer
rawdata = read.table("data.txt",sep='\t',header=TRUE)
group = factor(c(rep(1,2), rep(2,2))) 
y = DGEList(counts=rawdata[,2:ncol(rawdata)], genes=rawdata[,1],group=group)

y$samples$lib.size = colSums(y$counts)

delivery = group
data.frame(Sample=colnames(y),delivery)

design = model.matrix(~0+delivery)
rownames(design) = colnames(y)

y = estimateGLMCommonDisp(y, design, verbose=TRUE)
y = estimateGLMTrendedDisp(y,design)
y = estimateGLMTagwiseDisp(y,design)

#GLM
fit = glmFit(y, design)
contrast.design= c(1,-1)
lrt = glmLRT(fit, contrast=contrast.design)
glm = topTags(lrt, n = Inf)$table
glm_tags = topTags(lrt, n = Inf)$table %>% filter(FDR < 0.05 & abs(logFC) >1)

rawdata_norm = as.data.frame(t(t(rawdata[,-1])/colSums(rawdata[,-1])*1000000))
data= cbind(rawdata$id, rawdata_norm)
#Filter data and save only genes of interest
idx = match(glm_tags$genes, rawdata$id)
data_filtered= cbind(rawdata$id[idx], rawdata_norm[idx,])
write.csv(glm_tags,'sample.edgeR.csv', row.names = FALSE)
write.csv(data_filtered, 'sample.csv', row.names = FALSE)
```
