python trans_neg2pos.py 3neg.site.txt > 3neg.site.rev
python trans_neg2pos.py 4neg.site.txt > 4neg.site.rev
cat 3pos.site.rev 3neg.site.rev | sort -k1,1 -k2,2n > final.3.sites
cat 4pos.site.rev 4neg.site.rev | sort -k1,1 -k2,2n > final.4.sites


python normal_RT.py final.3.sites final.4.sites F-2-3.sorted.bam.filter.bam.region.bam F-2-4.sorted.bam.filter.bam.region.bam > RT_value.txt

