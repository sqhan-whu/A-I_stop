#!/usr/bin/env perl

use Cwd qw(abs_path);
use File::Basename;

$bam = shift;
#$bam =abs_path($bam);
$ANNO='/home/weiqi/project/00.DATABASE/GRmm38/Annotation/gencode.vM24.annotation.bed';
$GENOME = '/home/weiqi/project/00.DATABASE/GRmm38/GRCm38.p6.fa';
$SIZE = '/home/weiqi/project/00.DATABASE/GRmm38/GRCm38.p6.fa.size';
$BIN='/home/weiqi/workfs/20200808_4sU/scripts';
$SNP = '~/workfs/20200808_4sU/mapping/snp.con_processed.txt';
$CUTOFF = 50;
my ($fname,$path,$ext)=fileparse($bam,'.coincide.sort.bam');
$path =~ s/\/\s*$//;
open F, ">$fname\_cal_snp.sh";

print F "#!/bin/bash\n#SBATCH -N 1 -c 4\n";
#print F "bamtools filter -in $bam -out $path/$fname.sorted.filter.bam -isPrimaryAlignment true -isProperPair true  -isMapped true -mapQuality '>1'\n";

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
