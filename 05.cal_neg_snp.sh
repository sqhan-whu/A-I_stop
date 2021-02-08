#!/bin/bash
python TC.py F-2-1.sorted.bam.filter.bam.region.bam.neg.sort.bam /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa.size F-2-1.pos.TC.bam > 1neg.rough.txt
python TC.py F-2-2.sorted.bam.filter.bam.region.bam.neg.sort.bam /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa.size F-2-2.pos.TC.bam > 2neg.rough.txt
python TC.py F-2-3.sorted.bam.filter.bam.region.bam.neg.sort.bam /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa.size F-2-3.pos.TC.bam > 3neg.rough.txt
python TC.py F-2-4.sorted.bam.filter.bam.region.bam.neg.sort.bam /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa.size F-2-4.pos.TC.bam > 4neg.rough.txt

python cal_RT_neg.py 3neg.rough.txt|awk '{if($6 >=4 && $6/$5 >=0.1)print $0}' > 3neg.site.txt
python cal_RT_neg.py 4neg.rough.txt|awk '{if($6 >=4 && $6/$5 >=0.1)print $0}' > 4neg.site.txt

cat 3neg.site.txt 4neg.site.txt |awk '{if($3 ~/T/ && $4 ~/C/)print $0}'|cut -f 1,2 |sort |uniq > neg.sites
