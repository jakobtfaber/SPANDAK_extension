import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
import csv
sys.path.insert(0, os.path.abspath('../extractor'))

def read_data(database="database.csv"):
	
	#Read in -tentative 'database'- .csv
	filepaths = pd.read_csv(database)
	filpaths = filepaths.iloc[:,0]
	rawpaths = filepaths.iloc[1, :][1:]
	fieldnames = filepaths.columns[1:]
	csvs = []

	for csv in os.listdir("../pipeline_playground/"):
		if csv.endswith("4chan.csv"):
			csvs.append(csv)

	return filepaths, filpaths, rawpaths, fieldnames, csvs


def parse_spandak(csvs):

	for csv in csvs:
		
		#Read in .csv produced by SPANDAK (on FRB121102 as an example)
		cands = pd.read_csv(csv)

		#Isolate B-valued candidates (probable if not using ML functionality)
		B_idx = cands[cands.iloc[:, :]['Category']=='B'].index.values

		#Find filenames for B candidtaes
		files = [[i for i in cands.loc[:, :]['filename']][b] for b in B_idx]

		#Find DMs for B candidates
		DMs = [[i for i in cands.loc[:, :]['DM']][b] for b in B_idx]

		#Find time constraints
		#Time widths
		time_widths = [[i for i in cands.loc[:, :]['WIDTH']][b] for b in B_idx]

		#Time centers
		parse_center = [[i for i in [j.split('_') for j in cands.iloc[:, :]['PNGFILE']]][b][2] for b in B_idx]
		time_stamps = [np.float(m[:-3]) for m in parse_center]

		#Time limits
		start_times = [time_stamps[b]-time_widths[b] for b in np.arange(len(B_idx))]
		end_times = [time_stamps[b]+time_widths[b] for b in np.arange(len(B_idx))]

		return B_idx, files, DMs, time_widths, start_times, end_times

def extract_auto(filpaths, rawpaths, fieldnames, B_idx, files, DMs, time_widths, start_times, end_times):

	#Parse and Form Raw Voltage Extraction Commands

	extract_run_commands = []

	for B in np.arange(len(B_idx)):
		for raw in np.arange(len(rawpaths)):
			extract_run = 'python ' + '../extractor/extract_blocks.py ' + rawpaths[raw] \
			 + ' ' + 'blc' + str(fieldnames[raw][3:]) + files[B][33:-25] + ' ' \
			 + str(start_times[B]) + ' ' + str(end_times[B]) \
			 + ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/rawfiles/' \
			 + str(start_times[B]) + '_' + str(end_times[B])
			extract_run_commands.append(extract_run)

	return extract_run_commands


def splice_auto(filpaths, rawpaths, fieldnames, B_idx, files, DMs, time_widths, start_times, end_times):

	#Parse and Form Raw File Splicing Commands

	splicer_run_commands = []

	for B in np.arange(len(B_idx)):	
		splicer_run = 'python ' + '../extractor/splicer_raw.py ' + 'datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/rawfiles/' \
		+ str(start_times[B]) + '_' + str(end_times[B]) + '/ ' + '9 ' + 'spliced' +  files[0][33:-25]
		splicer_run_commands.append(splicer_run)

	return splicer_run_commands

#def cdd_auto(filpaths, rawpaths, fieldnames, B_idx, files, DMs, time_widths, start_times, end_times):
#
#	#Specify CDD parameters
#
#	cepoch = 58178
#	output = files[0][57:-25]
#	polar = 4
#	phasebin = 32768
#	p = 0
#	chan = 3712
#	samples = 558 #number of MB
#
#	#Parse and Form Coherent Dispersion Commands
#
#	cdd_run_commands = []
#
#	for B in np.arange(len(B_idx)):
#		cdd_run = 'dspsr ' + '-T ' + str(time_widths[B]) + ' -c ' + str(time_widths[B]) \
#		+ ' -cepoch ' + str(cepoch) + ' -O ' + str(output) + '_' + str(start_times[B]) + '_' \
#		+ str(end_times[B]) + ' -D ' + str(DMs[B]) + ' -d ' + str(polar) + ' -b ' \
#		+ str(phasebin) + ' -K -F ' + str(chan) + ':D ' + '-p 0 ' + '-U ' + str(samples) \
#		+ ' ../pipeline_playground/rawfiles/' + str(start_times[B]) + '_' + str(end_times[B]) + ' spliced' +  files[0][33:-25]
#		cdd_run_commands.append(cdd_run)
#
#	return cdd_run_commands



def main():
	
	filepaths, filpaths, rawpaths, fieldnames, csvs = read_data()
	B_idx, files, DMs, time_widths, start_times, end_times = parse_spandak(csvs)
	extract_run_commands = extract_auto(filpaths, rawpaths, fieldnames, B_idx, files, DMs, time_widths, start_times, end_times)
	splicer_run_commands = splice_auto(filpaths, rawpaths, fieldnames, B_idx, files, DMs, time_widths, start_times, end_times)
	#cdd_run_commands = cdd_auto(filpaths, rawpaths, fieldnames, B_idx, files, DMs, time_widths, start_times, end_times)

	#print(extract_run_commands[0])
	#Extract Raw Voltages

	for erc in extract_run_commands:
		os.system(erc)
#
	##Splice Raw Files Into Contiguous Raw File
#
	for src in splicer_run_commands:
		os.system(src)
	#
	##Coherently Dedisperse Raw File With DSPSR
#
	#for cdd in cdd_run_commands:
	#	os.system(cdd)

main()









