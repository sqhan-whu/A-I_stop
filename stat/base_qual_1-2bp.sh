seqtk sample -s100 HepG2-1_combined_R1.fastq.gz 1000000 > sub1.fq
seqtk sample -s100 H293-1_combined_R1.fastq.gz 1000000 > sub293.fq

python ~/bin_python/I_script/base_qual_12.py sub1.fq 6 > base_HepG2.txt
python ~/bin_python/I_script/base_qual_12.py sub293.fq 0 > base_H293T.txt

awk '{print $0"\tHepG2"}' base_HepG2.txt > HepG2.txt
awk '{print $0"\tHepG2"}' base_H293T.txt > base_H293T.txt

cat base_H293T.txt HepG2.txt > qual.txt

#######################

########C:\Users\50246\Desktop\scratch\03.weiqi\figure\quality_2base

library(ggplot2) 
library(ggpubr)

data <- read.table("qual.txt",head=T)
pdf('qual.pdf', w= 12/2.54, h= 12/2.54, pointsize= 10)
ggbarplot(data, x = "base", y = "qual", color="group",
          add = c("mean_se"),
           palette = "jco",
          position = position_dodge(0.8))
dev.off()
