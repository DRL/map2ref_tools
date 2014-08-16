#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File   	: cas_to_barchart.py
Version : 0.1
Author 	: Dominik R. Laetsch, dominik.laetsch at gmail dot com 
Bugs 	: ...
To do 	: ...
"""
from __future__ import division
import os, re, sys, argparse, subprocess, itertools, commands
import numpy as np

def parse_assembly(assembly_file):
	'''Parses assembly file.'''

	assembly_contigs = []
	
	contig_id, contig_seq = '', ''

	with open(assembly_file) as fh:
		for line in fh:
			if line.startswith(">"):
				if contig_id == '':
					pass
				else:
					contig_n = contig_seq.count('N')
					contig_len = len(contig_seq)
					contig_atgc = contig_len - contig_n
					contig_gc = (contig_seq.count('G') + contig_seq.count('C')) / (contig_len - contig_n)
					assembly_contigs.append([(contig_id), (contig_len), (contig_atgc), (contig_n), (contig_gc)])
				contig_id = line.lstrip(">").rstrip("\n")
				contig_seq = ''
			else:
				contig_seq += line.rstrip("\n")
	contig_n = contig_seq.count('N')
	contig_len = len(contig_seq)
	contig_atgc = contig_len - contig_n
	contig_gc = (contig_seq.count('G') + contig_seq.count('C')) / (len(contig_seq) - contig_n)
	assembly_contigs.append([(contig_id), (contig_len), (contig_atgc), (contig_n), (contig_gc)])
	return assembly_contigs

def parse_cas(cas_file):
	assembly_re = re.compile(r"-d (\S+)")
	cov_start_re = re.compile(r"Contig info")
	contig_info_re = re.compile(r"\s+Contig (\d+) info")
	total_sites_re = re.compile(r"\s+Sites\s+(\d+)")
	cov_0_sites_re = re.compile(r"\s+Covered\s+0\s+times\s+(\d+)")
	cov_1_sites_re = re.compile(r"\s+Covered\s+1\s+time\s+(\d+)")
	cov_2_sites_re = re.compile(r"\s+Covered\s+2\s+times\s+(\d+)")
	cov_3_sites_re = re.compile(r"\s+Covered\s+3\+\s+times\s+(\d+)")
	avg_cov_re = re.compile(r"\s+Average coverage\s+(\d+\.\d+)")
	#error, message = commands.getstatusoutput("clc_mapping_info " + cas_file)
	#if (error):
	#	sys.exit("ERROR: Please load the CLC module ('module load clc')") 
	#mapping_info = subprocess.check_output("clc_mapping_info " + cas_file, stderr=subprocess.STDOUT, shell=True)
	read_switch = 0
	data_list = []
	contig_id, sites, cov_0_sites, cov_1_sites, cov_2_sites, cov_3_sites = '', 0, 0, 0, 0, 0
	with open(cas_file) as fh:
		for line in fh:
	#for line in mapping_info:
			line = line.rstrip("\n")
			if cov_start_re.match(line):
				read_switch = 1
			if read_switch:
				if contig_info_re.match(line):
					if not contig_id == '':
						data_list.append([(contig_id), (sites), (cov_0_sites), (cov_1_sites), (cov_2_sites), (cov_3_sites)])
					contig_id = int(contig_info_re.search(line).group(1))
				elif total_sites_re.match(line):
					sites = int(total_sites_re.search(line).group(1))
				elif cov_0_sites_re.match(line):
					cov_0_sites = int(cov_0_sites_re.search(line).group(1))
				elif cov_1_sites_re.match(line):
					cov_1_sites = int(cov_1_sites_re.search(line).group(1))
				elif cov_2_sites_re.match(line):
					cov_2_sites = int(cov_2_sites_re.search(line).group(1))
				elif cov_3_sites_re.match(line):
					cov_3_sites = int(cov_3_sites_re.search(line).group(1))
				#elif avg_cov_re.match(line):
				#	avg_cov = float(avg_cov_re.search(line).group(1))
				else:
					pass
	data_list.append([(contig_id), (sites), (cov_0_sites), (cov_1_sites), (cov_2_sites), (cov_3_sites)])
	return data_list
	#print array

def merge_data(assembly_data, mapping_data):
	data = []
	for assembly_info, mapping_info in zip(assembly_data, mapping_data):
		if assembly_info[1] == mapping_info[1]:
			#print assembly_info 
			#print mapping_info 
			percent_atgc_0_cov = float((mapping_info[2]-assembly_info[3])/assembly_info[2])
			data.append([(assembly_info[0]),int(assembly_info[1]),int(assembly_info[2]),int(assembly_info[3]),float(assembly_info[4]),int(mapping_info[2]),int(mapping_info[3]),int(mapping_info[4]),int(mapping_info[5]),float(percent_atgc_0_cov)])
		else:
			sys.exit("ERROR : " + assembly_info + "\n" + mapping_info)
	return data

def get_input():

	parser = argparse.ArgumentParser(
		prog='cas_to_barchart.py',
		usage = '%(prog)s infile [-p] [-f] [-t] [-e] [-n] [-o] [-l] [-c] [-s] [-h]',
		add_help=True)
	parser.add_argument('i', metavar = 'infile', help='Input file (blobplot.txt)')
	parser.add_argument('-p', metavar = 'max_taxa_plot', default=7, type = int, help='Maximum number of phyla to plot (Default = 7)')
	args = parser.parse_args()

def print_barcharts(array):
	from stackedBarGraph import StackedBarGrapher 
	from matplotlib import pyplot as plt

	d = np.array([[101.,0.,0.,0.,0.,0.,0.],
              [92.,3.,0.,4.,5.,6.,0.],
              [56.,7.,8.,9.,23.,4.,5.],
              [81.,2.,4.,5.,32.,33.,4.],
              [0.,45.,2.,3.,45.,67.,8.],
              [99.,5.,0.,0.,0.,43.,56.]])

	d_widths = [.5,1.,3.,2.,1.,2.]
	d_labels = ["fred","julie","sam","peter","rob","baz"]
	d_colors = ['#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777']
	fig = plt.figure()
	SBG = StackedBarGrapher()
	ax = fig.add_subplot(111)
	SBG.stackedBarPlot(ax,
                   d,
                   d_colors,
                   xLabels=d_labels,
                   yTicks=7,
                   scale=True
                  )
	plt.title("Scaled bars with set widths")

	fig.subplots_adjust(bottom=0.4)
	plt.tight_layout()
	plt.show()
	plt.close(fig)
	del fig

if __name__ == '__main__':
	assembly_file = sys.argv[1]
	assembly_data = parse_assembly(assembly_file)
	cas_file = sys.argv[2]
	mapping_data = parse_cas(cas_file)
	data = merge_data(assembly_data, mapping_data)
	array = np.array(data)
	new_array = array[array[:,1].argsort()]
	print new_array
	#print_barcharts(array)
	#SBG.demo()
