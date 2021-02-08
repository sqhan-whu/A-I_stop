#!/bin/bash
python AG.py F-2-1.sorted.bam.filter.bam.region.bam.pos.sort.bam /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa.size F-2-1.pos.AG.bam > 1.rough.txt
python AG.py F-2-2.sorted.bam.filter.bam.region.bam.pos.sort.bam /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa.size F-2-2.pos.AG.bam > 2.rough.txt
python AG.py F-2-3.sorted.bam.filter.bam.region.bam.pos.sort.bam /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa.size F-2-3.pos.AG.bam > 3.rough.txt
python AG.py F-2-4.sorted.bam.filter.bam.region.bam.pos.sort.bam /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa /home/weiqi/project/01.weiqi/20200929_inosine/process/index/GRCm38.p6.fa.size F-2-4.pos.AG.bam > 4.rough.txt

python cal_RT.py 3.rough.txt|awk '{if($6 >=4 && $6/$5 >=0.1)print $0}' > 3.site.txt
python cal_RT.py 4.rough.txt|awk '{if($6 >=4 && $6/$5 >=0.1)print $0}' > 4.site.txt

cat 3.site.txt 4.site.txt |awk '{if($3 ~/A/ && $4 ~/G/)print $0}'|cut -f 1,2 |sort |uniq > pos.sites
