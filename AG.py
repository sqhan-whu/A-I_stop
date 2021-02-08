#!/usr/bin/env python

import sys
import pysam
from collections import Counter
from operator import itemgetter
import numpy as np

# python script.py test.bam hg38.fasta hg38.sizes
script, sam_path, fasta_path, size_path, new_bam = sys.argv

def get_pileup_database(sam_path, size_path):
	sam = pysam.AlignmentFile(sam_path)
	sizes = dict()
	with open(size_path) as f:
		for line in f:
				row = line.strip("\n").split("\t")
				sizes[row[0]] = int(row[1])
	values = dict()
	for ref in sam.references:
		print(ref)
		length = sizes[ref]
		dat = np.zeros(length)
		for x in sam.pileup(ref, 0, length):
			dat[x.pos] = x.n
		values[ref] = dat
	return values

def get_mutation_site(sam_path, fasta_path, values):
	sam = pysam.AlignmentFile(sam_path)
	fasta = pysam.FastaFile(fasta_path)
	new_sam = pysam.AlignmentFile(new_bam, "wb", template=sam)
	for align in sam:
		read_mate = 'read1' if align.is_read1 else 'read2'
		if align.is_unmapped:
			continue
		qseq = align.query_alignment_sequence
		#qual = align.query_alignment_qualities

		qual = align.query_qualities
		name = align.query_name
		i = j = l = 0
		for mark, num in align.cigartuples:
			if mark == 0:
				seq1 = qseq[i:i+num].upper()
				seq2 = fasta.fetch(sam.references[align.reference_id], 
								   align.reference_start + j, 
								   align.reference_start + j + num).upper()
				qual1 = qual[i:i+num]
				k = 0
				for base1, base2, quality in zip(seq1, seq2, qual1):
					if base1 != base2:
						chrom = sam.references[align.reference_id]
						start = align.reference_start + j + k
						end = start + 1
						rbase = base2
						qbase = base1
						coverage = values[chrom][start]
						n_stop = values[chrom][start+i+num-k-1]
						#n_through = values[chrom][start+2]
						# OUTPUT
						print("\t".join(map(str, [chrom, start, end, rbase, qbase, coverage,n_stop ,k+1, i+num, quality, start+i+num-k, name, read_mate])))
						
						if rbase == "A" and qbase == "G" and l == 0:
							#pass
							new_sam.write(align)
							l = 1
						
					k += 1
				i += num
				j += num
			elif mark == 1:
				i += num
			elif mark == 2 or mark == 3:
				j += num
			elif mark == 4 or mark == 5:
				i += num
				j += num
			else:
				raise Exception()
	
get_mutation_site(sam_path, fasta_path, get_pileup_database(sam_path, size_path))
