#!/bin/bash
python call.rpmf.peak.py -i 3 -c 1 -w 3
python call.rpmf.peak.py -i 4 -c 2 -w 3
python peak_finalize.py -p 31RPMF_out.txt -a ~/project/00.DATABASE/GRmm38/Annotation/mRNA/gencode.vM24.gene.slidingwindows.gtf
python peak_finalize.py -p 42RPMF_out.txt -a ~/project/00.DATABASE/GRmm38/Annotation/mRNA/gencode.vM24.gene.slidingwindows.gtf
