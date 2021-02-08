#!/usr/bin/env python
import csv
import sys, getopt
import os
import HTSeq
from collections import defaultdict

###  Setup all cut-off

# RPMF cut-off
RPMF_CUTOFF = 4
# window occurence
#OCCURENCE_CUTOFF = 2

def get_peak_position(peak_file, annotation_file):
	dict_peaks = dict()
	PATH_PEAKS = peak_file
	dict_peak_window = defaultdict(list)
	with open(peak_file) as PATH_PEAKS:
		PATH_PEAKS.readline()
		for row in PATH_PEAKS:
			if row.startswith('WindowId'):
				continue
			else:
	#		print(row.split('\t')[6].strip())
				RPMF = row.split('\t')[6].strip()
	#	print(row.split('\t')[6])
				if float(RPMF) > RPMF_CUTOFF:
					dict_peaks[row.split('\t')[0]] = 0
	# annotate peak_window
	dict_peak_window = defaultdict(list)
	PATH_WINDOWS = annotation_file

	gtf_peaks_file = HTSeq.GFF_Reader(PATH_WINDOWS)
	k = 0
	for feature in gtf_peaks_file:
		window_id = feature.attr['ID']
		if window_id in dict_peaks:
			if window_id not in dict_peak_window:
				k += 1
					#if k%5000 == 0:
				#	print(k,len(dict_peaks),feature.iv)
				begin_window = feature.iv.start
				end_window = feature.iv.end
				chromo_window = feature.iv.chrom
				strand_window = feature.iv.strand
				dict_peak_window[window_id].append(chromo_window)
				dict_peak_window[window_id].append(str(begin_window))
				dict_peak_window[window_id].append(str(end_window))
				type_transcript = window_id.split('_window_')[1].rsplit('_', 1)[0]
				dict_peak_window[window_id].append(type_transcript)
				transcript_name = window_id.split('_window_')[0].rsplit('-', 1)[0]
				dict_peak_window[window_id].append(transcript_name)
				dict_peak_window[window_id].append(strand_window)

	with open(peak_file +  '_windows.bed', 'w') as annot_file:
		for key, value in dict_peak_window.items():
			#print(key,value)
			if key == 'None':
				print('none',value)
			annot_file.write('\t'.join(value) + '\n')

def merge_peaks(PATH_PEAKS):
	bedfile = PATH_PEAKS +  '_windows.bed'
	sort_bedfile = PATH_PEAKS + '_windows_sort.bed'
	merge_bedfile = PATH_PEAKS + '_windows_merge.bed'
	final_bedfile = PATH_PEAKS + '_final.bed'
	sort_command = 'sort -k1,1 -k2,2n '+bedfile+' > '+sort_bedfile
	os.system(sort_command)
	merge_command = 'bedtools merge -c 4,5 -o distinct -i '+sort_bedfile+' > '+merge_bedfile
	os.system(merge_command)
	awk_command = 'awk \'{OFS = "\\t"; print $1,$2,$3,$4,"","+"}\' '+merge_bedfile + ' > ' + final_bedfile
	os.system(awk_command)

def main(argv):
	'''
	Main function of peak_calling_fisher.py
	:param argv:
	-p --path Path of the working folder
	-a --annotation Annotation file in GTF format
	:return:
	'''
	try:
		opts, args = getopt.getopt(argv,"hp:a:", ["path=","annotation="])
	except getopt.GetoptError:
		print('Cannot run command - Help: peak_finalize.py -p <path> -e <expdesign> -t <peak_technique> -a <annotation> -b <bed_name>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('peak_finalize.py -p <path> -a <annotation>')
			sys.exit()
		elif opt in ("-p", "--path"):
			path = arg
		elif opt in ("-a", "--annotation"):
			annotation_file = arg
	get_peak_position(path, annotation_file)
	merge_peaks(path)
#	print(path,annotation_file)


if __name__ == "__main__":
	main(sys.argv[1:])
