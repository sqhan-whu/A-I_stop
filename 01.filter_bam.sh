#!/bin/bash
bamtools filter -in F-2-1.sorted.bam -out F-2-1.sorted.bam.filter.bam -isPrimaryAlignment true -isProperPair true -isMapped true -mapQuality '>1'
bamtools filter -in F-2-2.sorted.bam -out F-2-2.sorted.bam.filter.bam -isPrimaryAlignment true -isProperPair true -isMapped true -mapQuality '>1'
bamtools filter -in F-2-3.sorted.bam -out F-2-3.sorted.bam.filter.bam -isPrimaryAlignment true -isProperPair true -isMapped true -mapQuality '>1'
bamtools filter -in F-2-4.sorted.bam -out F-2-4.sorted.bam.filter.bam -isPrimaryAlignment true -isProperPair true -isMapped true -mapQuality '>1'
samtools index F-2-1.sorted.bam.filter.bam
samtools index F-2-2.sorted.bam.filter.bam
samtools index F-2-3.sorted.bam.filter.bam
samtools index F-2-4.sorted.bam.filter.bam
samtools view -c F-2-1.sorted.bam.filter.bam >> library.size
samtools view -c F-2-2.sorted.bam.filter.bam >> library.size
samtools view -c F-2-3.sorted.bam.filter.bam >> library.size
samtools view -c F-2-4.sorted.bam.filter.bam >> library.size
htseq-count -i ID -t window -s no -m union --nonunique all -f 'bam' F-2-1.sorted.bam.filter.bam /home/weiqi/project/00.DATABASE/GRmm38/Annotation/mRNA/gencode.vM24.gene.slidingwindows.gtf > HTSeq_Count_1_windows.txt
htseq-count -i ID -t window -s no -m union --nonunique all -f 'bam' F-2-2.sorted.bam.filter.bam /home/weiqi/project/00.DATABASE/GRmm38/Annotation/mRNA/gencode.vM24.gene.slidingwindows.gtf > HTSeq_Count_2_windows.txt
htseq-count -i ID -t window -s no -m union --nonunique all -f 'bam' F-2-3.sorted.bam.filter.bam /home/weiqi/project/00.DATABASE/GRmm38/Annotation/mRNA/gencode.vM24.gene.slidingwindows.gtf > HTSeq_Count_3_windows.txt
htseq-count -i ID -t window -s no -m union --nonunique all -f 'bam' F-2-4.sorted.bam.filter.bam /home/weiqi/project/00.DATABASE/GRmm38/Annotation/mRNA/gencode.vM24.gene.slidingwindows.gtf > HTSeq_Count_4_windows.txt
