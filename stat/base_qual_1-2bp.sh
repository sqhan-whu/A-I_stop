seqtk sample -s100 HepG2-1_combined_R1.fastq.gz 1000000 > sub1.fq
seqtk sample -s100 H293-1_combined_R1.fastq.gz 1000000 > sub293.fq

python ~/bin_python/I_script/base_qual_12.py sub1.fq 6 > base_HepG2.txt
python ~/bin_python/I_script/base_qual_12.py sub293.fq 0 > base_H293T.txt
