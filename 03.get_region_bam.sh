#!/bin/bash
bedtools  intersect -a F-2-1.sorted.bam.filter.bam -b 31RPMF_out.txt_final.bed  > F-2-1.sorted.bam.filter.bam.region.bam
bedtools  intersect -a F-2-2.sorted.bam.filter.bam -b 42RPMF_out.txt_final.bed  > F-2-2.sorted.bam.filter.bam.region.bam
bedtools  intersect -a F-2-3.sorted.bam.filter.bam -b 31RPMF_out.txt_final.bed  > F-2-3.sorted.bam.filter.bam.region.bam
bedtools  intersect -a F-2-4.sorted.bam.filter.bam -b 42RPMF_out.txt_final.bed  > F-2-4.sorted.bam.filter.bam.region.bam
samtools index F-2-1.sorted.bam.filter.bam.region.bam
samtools index F-2-2.sorted.bam.filter.bam.region.bam
samtools index F-2-3.sorted.bam.filter.bam.region.bam
samtools index F-2-4.sorted.bam.filter.bam.region.bam
