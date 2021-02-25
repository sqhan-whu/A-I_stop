#!/bin/bash
python trans_neg2pos.py H293-1.sorted.bam.filter.bam.pos.sort.bam.pos.rough.txt H293-1.sorted.bam.filter.bam.neg.sort.bam.neg.rough.txt > H293-1.site.txt
python trans_neg2pos.py H293-2.sorted.bam.filter.bam.pos.sort.bam.pos.rough.txt H293-2.sorted.bam.filter.bam.neg.sort.bam.neg.rough.txt > H293-2.site.txt
python trans_neg2pos.py H293-EnV-1.sorted.bam.filter.bam.pos.sort.bam.pos.rough.txt H293-EnV-1.sorted.bam.filter.bam.neg.sort.bam.neg.rough.txt > H293-EnV-1.site.txt
python trans_neg2pos.py H293-EnV-2.sorted.bam.filter.bam.pos.sort.bam.pos.rough.txt H293-EnV-2.sorted.bam.filter.bam.neg.sort.bam.neg.rough.txt > H293-EnV-2.site.txt
python trans_neg2pos.py H293-Oxd-1.sorted.bam.filter.bam.pos.sort.bam.pos.rough.txt H293-Oxd-1.sorted.bam.filter.bam.neg.sort.bam.neg.rough.txt > H293-Oxd-1.site.txt
python trans_neg2pos.py H293-Oxd-2.sorted.bam.filter.bam.pos.sort.bam.pos.rough.txt H293-Oxd-2.sorted.bam.filter.bam.neg.sort.bam.neg.rough.txt > H293-Oxd-2.site.txt
python trans_neg2pos.py Hela-1.sorted.bam.filter.bam.pos.sort.bam.pos.rough.txt Hela-1.sorted.bam.filter.bam.neg.sort.bam.neg.rough.txt > Hela-1.site.txt
python trans_neg2pos.py Hela-2.sorted.bam.filter.bam.pos.sort.bam.pos.rough.txt Hela-2.sorted.bam.filter.bam.neg.sort.bam.neg.rough.txt > Hela-2.site.txt
python trans_neg2pos.py Hela-EnV-1.sorted.bam.filter.bam.pos.sort.bam.pos.rough.txt Hela-EnV-1.sorted.bam.filter.bam.neg.sort.bam.neg.rough.txt > Hela-EnV-1.site.txt
python trans_neg2pos.py Hela-EnV-2.sorted.bam.filter.bam.pos.sort.bam.pos.rough.txt Hela-EnV-2.sorted.bam.filter.bam.neg.sort.bam.neg.rough.txt > Hela-EnV-2.site.txt
python trans_neg2pos.py Hela-Oxd-1.sorted.bam.filter.bam.pos.sort.bam.pos.rough.txt Hela-Oxd-1.sorted.bam.filter.bam.neg.sort.bam.neg.rough.txt > Hela-Oxd-1.site.txt
python trans_neg2pos.py Hela-Oxd-2.sorted.bam.filter.bam.pos.sort.bam.pos.rough.txt Hela-Oxd-2.sorted.bam.filter.bam.neg.sort.bam.neg.rough.txt > Hela-Oxd-2.site.txt

##################### A  ------>  AG AT AC
python stat.site.mutation.py H293-1.site.txt H293-2.site.txt H293-EnV-1.site.txt H293-EnV-2.site.txt H293-Oxd-1.site.txt H293-Oxd-2.site.txt Hela-1.site.txt Hela-2.site.txt Hela-EnV-1.site.txt Hela-EnV-2.site.txt Hela-Oxd-1.site.txt Hela-Oxd-2.site.txt > site_mutaion_total.A.txt
python p2.stat.site.mutation.py H293-1.site.txt H293-2.site.txt H293-EnV-1.site.txt H293-EnV-2.site.txt H293-Oxd-1.site.txt H293-Oxd-2.site.txt Hela-1.site.txt Hela-2.site.txt Hela-EnV-1.site.txt Hela-EnV-2.site.txt Hela-Oxd-1.site.txt Hela-Oxd-2.site.txt > p2.site_mutaion_total.A.txt
