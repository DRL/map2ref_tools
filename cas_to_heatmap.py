#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
File   	: cas_to_heatmap.py
Version : 0.1
Author 	: Dominik R. Laetsch, dominik.laetsch at gmail dot com 
Bugs 	: ...
To do 	: ...
"""
from __future__ import division
import os, re, sys, argparse, subprocess, itertools, commands
import numpy as np
import matplotlib as mat
import matplotlib.pyplot as plt
from matplotlib import cm
from pylab import imshow, show, get_cmap

def get_input():

	parser = argparse.ArgumentParser(
		prog='cas_to_heatmap.py',
		usage = '%(prog)s -a ASSEMBLY -c CAS ... [-h]',
		add_help=True)
	parser.add_argument('-a', metavar = 'assembly', help='Assembly file')
	parser.add_argument('-c', metavar = 'cas', nargs='+', help='Cas file(s)') 
	args = parser.parse_args()
	assembly_file, cas_files = args.a, args.c
	if not assembly_file:
		sys.exit("ERROR : Please specify an assembly file.")
	if len(cas_files) < 1:
		sys.exit("ERROR : Please specify a cas file.")
	return assembly_file, cas_files

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
	mapping_data = []
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
						mapping_data.append([(contig_id), (sites), (cov_0_sites), (cov_1_sites), (cov_2_sites), (cov_3_sites)])
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
	mapping_data.append([(contig_id), (sites), (cov_0_sites), (cov_1_sites), (cov_2_sites), (cov_3_sites)])
	return mapping_data
	#print array

def merge_data(assembly_data, mapping_data):
	data = []
	for assembly_info, mapping_info in zip(assembly_data, mapping_data):
		if assembly_info[1] == mapping_info[1]:
			percent_atgc_0_cov = float((mapping_info[2]-assembly_info[3])/assembly_info[2])
			# id, sites, ATGC, N, GC, cov0, cov, cov2, cov3, 0cov/ATGC
			data.append([(assembly_info[0]),(int(assembly_info[1])),(int(assembly_info[2])),(int(assembly_info[3])),(float(assembly_info[4])),(int(mapping_info[2])),(int(mapping_info[3])),(int(mapping_info[4])),(int(mapping_info[5])),(float(percent_atgc_0_cov))])
		else:
			sys.exit("ERROR : " + assembly_info + "\n" + mapping_info)
	return data

def print_to_file(data, field):
	bucket = {}
	for dataset in data:
		bucket[dataset]=[]
		for line in data[dataset]:
			bucket[dataset].append("line")




def heatmap(data):
	plt.rcParams['xtick.major.pad']='50'
	plt.rcParams['ytick.major.pad']='50'
	black, grey, background_grey, white = '#262626', '#d3d3d3', '#F0F0F5', '#ffffff'
	#plt.figure(1, figsize=(250,10), dpi=1000)

	#colors = cm.get_cmap(name='Set2')
	#color_index = 1
	#color_dict={}
	#for dataset in data:
	#	color_dict[dataset] = mat.colors.rgb2hex(colors(1.0 * (color_index/len(data))))
	#	color_index += 1

	columns = 0
	rows = 0
	names = []
	for dataset in sorted(data.keys()):
		rows = len(data[dataset])
		columns += 1
		names.append(dataset)
	plt.figure(1, figsize=(100,50), dpi=100)

	plot_array = ''

	for dataset in sorted(data.keys()):
		print dataset
		array = np.array(data[dataset])
		#value = np.array(array[:,0]).astype(str)
		value = np.array(array[:,1]).astype(int) # length
		#value = np.array(array[:,2]).astype(int) # ATGCs
		#value = np.array(array[:,3]).astype(int) # N's
		#value = np.array(array[:,4]).astype(float) # GC
		#value = np.array(array[:,5]).astype(int) # cov0's
		#value = np.array(array[:,9]).astype(float) # 0cov/ATGC
		sort_idx = np.argsort(value)
		sort_idx = sort_idx[::-1]

		#outfile = "sorted_by_uncov_sites"
		
		#x_label = " Contigs ordered by decreasing unmapped nuc's"
		

		contig_id = np.array(array[:,0][sort_idx].astype(str))
		length = np.array(array[:,1][sort_idx].astype(int))
		gc = np.array(array[:,4][sort_idx].astype(float))
		n = np.array(array[:,3][sort_idx].astype(int))

		perc_uncov_sites = np.array(array[:,9][sort_idx]).astype(float)
		
		mask_for_cov = np.where(perc_uncov_sites <= 0.0)
		perc_uncov_sites[mask_for_cov] = 0.0

		print contig_id	
		print length
		print gc
		print n
		print perc_uncov_sites

		if plot_array == '':
			plot_array = np.array(array[:,9][sort_idx]).astype(float)
		else:
			plot_array = np.vstack((perc_uncov_sites, plot_array))
		
		print plot_array.shape
		#contig = np.arange(len(array))
		#idx = np.empty(len(array))
		#idx.fill(1)
		
	np.savetxt('out_test.txt', (plot_array.T), delimiter='\t', fmt='%.5f')
	outfile = "sorted_by_size"	
	x_label = "Contigs ordered by decreasing size"
		#plt.imshow((contig, perc_uncov_sites), cmap=get_cmap("Spectral"), interpolation=None)
	#extent = [1, columns, 1, rows]
	#plt.imshow(plot_array, cmap=get_cmap("jet"), aspect='auto', interpolation='none')
	plt.imshow(plot_array, cmap=get_cmap("gnuplot"), aspect='auto', interpolation='none')
	#plt.imshow(plot_array, cmap=get_cmap("afmhot_r"), aspect='auto', interpolation='none')
	#plt.imshow(plot_array, cmap=get_cmap("gnuplot2"), aspect='auto', interpolation='none')
	cbar = plt.colorbar()
	cbar.set_label('Unmapped nucleotides (%)',size=75, labelpad=100, rotation=270 )
	for t in cbar.ax.get_yticklabels():
		t.set_fontsize(75)
	#plt.set_xlabel("Contigs sorted by length", fontsize=40)
	#plt.set_ylabel("Dataset", fontsize=40)
	plt.yticks(np.arange(int(columns)), names, fontsize=75)
	plt.xticks( fontsize=75)
	plt.xlabel(x_label, fontsize=75)
	#plt.grid(which='minor', axis='y', color=white)
	plt.savefig(assembly_file + "." + outfile + ".png", format='png')
	plt.close()
	# Format
	

	# turn off the frame
	

	# Set the labels

	# note I could have used nba_sort.columns but made "labels" instead

	# rotate the

if __name__ == '__main__':
	assembly_file, cas_files = get_input()
	assembly_data = parse_assembly(assembly_file)
	data = {}
	for cas_file in cas_files:
		name = cas_file[0:9]
		print name
		cas_data = parse_cas(cas_file)
		data[name]= merge_data(assembly_data, cas_data)

	# id, sites, ATGC, N, GC, cov0, cov, cov2, cov3, 0cov/ATGC
	#print_to_file(data, 9)
	heatmap(data)
	