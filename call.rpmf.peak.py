#!/usr/bin/env python
import HTSeq
import csv
import sys, getopt
import numpy
import os
import re

def calculate_rpm(PATH_LIBRARY_SIZE,data):
	with open(PATH_LIBRARY_SIZE, "r") as file_library_size:
		library_size = 0.0
		for row in file_library_size:
			bam_name = row.split('\t')[0]
			if bam_name == data:
				library_size = float(list(filter(None, re.split('\t', row)))[1].strip())/1000000
				print('RPM FACTOR :', bam_name, library_size)
	with open('HTSeq_Count_' + data + '_RPMF_RPM.txt', "w") as listPeaks:
		header = 'WindowId\tWindowcov\tRPM'
		listPeaks.write(header+'\n')
		fileName = 'HTSeq_Count_' + data + '_windows.txt'
		print('read '+ fileName)
		with open(fileName) as windowCovFile:
			for row in windowCovFile:
				window_id = row.split('\t')[0]
				if len(re.split('\t', row)) == 1:
					window_cov = 0
				else:
					window_cov = float(re.split('\t', row)[1])
				norm_cov = (window_cov / library_size)
				listPeaks.write(window_id+'\t'+str(window_cov)+'\t'+str(norm_cov)+'\n')
	print("RPM calculated")

def calculate_rpm_fold( ip_data, input_data, window_cutoff):
	with open('HTSeq_Count_' + ip_data + '_RPMF_RPM.txt', "r") as ipfile, \
		open('HTSeq_Count_' + input_data + '_RPMF_RPM.txt', "r") as inputfile, \
		open(ip_data+ input_data+'RPMF_out.txt', "w") as bed_file:
		header = ['WindowId','Windowcov','RPM','Windowcov_Input','RPM_Input','Ratio_windowcov','RPMF']
		bed_file.write('\t'.join(header)+ '\n')
		ipfile.readline()
		inputfile.readline()
		window_name_to_row = dict()
		for rowInput in inputfile:
			window_id_input =  list(filter(None, re.split('\t', rowInput)))[0]
			window_name_to_row[window_id_input] = rowInput
		# RPM 2
		for row_ip in ipfile:
			window_id_ip = list(filter(None, re.split('\t', row_ip)))[0]
			row_input = window_name_to_row[window_id_ip]
			row_input = row_input.replace(window_id_ip,'').strip()
			new_row = row_ip.strip() + '\t' + row_input.strip()

			# Calc ratio window	
			window_cov = float(list(filter(None, re.split('\t', row_ip)))[1].strip())
			window_input_cov = float(list(filter(None, re.split('\t', row_input)))[0].strip())

			if window_input_cov == 0:
				new_row += '\t'+str(window_cov)
			else:
				ratio_windows = window_cov / window_input_cov
				new_row += '\t'+str(ratio_windows)

			# Calc ratio RPM
			if window_cov > window_cutoff:
				rpm_ip = re.split('\t', row_ip)[-1].strip()
				rpm_input = re.split('\t', row_input)[-1].strip()
				rpm_ip = float(rpm_ip)
				rpm_input = float(rpm_input)
				if rpm_input > 0:
					rpm = float(rpm_ip) / float(rpm_input)
				else:
					rpm = float(rpm_ip) / 0.1
				new_row += '\t'+str(rpm)

			else:
				new_row += '\t0'
			bed_file.write(new_row + '\n')
	print("RPMF Score calculated")

def main(argv):
	try:
		opts, args = getopt.getopt(argv,"i:c:w:", ["ip_data=", "input_data=", "window_cutoff="])
	except getopt.GetoptError:
		print('Cannot run command - Help: peak_calling_rpmf.py -p <path> -b <biocond> -i <ip_data> -c <input_data> -w <window_cutoff>')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print('peak_calling_rpmf.py -p <path> -b <biocond> -i <ip_data> -c <input_data> -w <window_cutoff>')
			sys.exit()
		elif opt in ("-i", "--ipdata"):
			ip_data = arg
		elif opt in ("-c", "--inputdata"):
			input_data = arg
		elif opt in ("-w", "--windcutoff"):
			window_cutoff = int(arg)

	PATH_LIBRARY_SIZE = 'library_size.txt'
	calculate_rpm(PATH_LIBRARY_SIZE, ip_data)
	calculate_rpm(PATH_LIBRARY_SIZE, input_data)
	calculate_rpm_fold(ip_data, input_data, window_cutoff)

if __name__ == "__main__":
	main(sys.argv[1:])
