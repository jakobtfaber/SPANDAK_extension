import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
import csv
sys.path.insert(0, os.path.abspath('../extractor'))

#with open('database.csv', 'r') as infile:
#	reader = csv.DictReader(infile)
#	fieldnames = reader.fieldnames[1:]

#Read in -tentative 'database'- .csv
filepaths = pd.read_csv("database.csv")
fieldnames = filepaths.columns[1:]
print(fieldnames)
filpaths = filepaths.iloc[:,0]
rawpaths = filepaths.iloc[1, :][1:]
csvs = []

for csv in os.listdir("../pipeline_playground/"):
	if csv.endswith("csv.csv"):
		csvs.append(csv)


def extractncdd(filpaths, rawpaths, fieldnames, csvs):

	for csv in csvs:
		#Read in .csv produced by SPANDAK (on FRB121102 as an example)
		cands = pd.read_csv(csv)

		#Isolate B-valued candidates (probable if not using ML functionality)
		B_vals = cands[cands.iloc[:, :]['Category']=='B'].index.values

		#Find filenames for B candidtaes
		files = [[i for i in cands.loc[:, :]['filename']][k] for k in B_vals]

		#Find DMs for B candidates
		DMs = [[i for i in cands.loc[:, :]['DM']][k] for k in B_vals]

		#Find time constraints
		#Time widths
		time_widths = [[i for i in cands.loc[:, :]['WIDTH']][k] for k in B_vals]

		#Time centers
		parse_center = [[i for i in [j.split('_') for j in cands.iloc[:, :]['PNGFILE']]][k][2] for k in B_vals]
		time_stamps = [np.float(m[:-3]) for m in parse_center]

		#Time limits
		start_times = [time_stamps[i]-time_widths[i] for i in np.arange(len(B_vals))]
		end_times = [time_stamps[i]+time_widths[i] for i in np.arange(len(B_vals))]

		#Parse and Form Raw Voltage Extraction Commands

		extract_run_commands = []

		for B in np.arange(len(B_vals)):
			for raw in np.arange(len(rawpaths)):
				extract_run = 'python ' + '../extractor/extract_blocks.py ' + rawpaths[raw] \
				 + ' ' + 'blc' + str(fieldnames[raw][3:]) + files[B][33:-25] + ' ' \
				 + str(start_times[B]) + ' ' + str(end_times[B]) \
				 + ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/rawfiles/' \
				 + str(start_times[B]) + '_' + str(end_times[B])
				extract_run_commands.append(extract_run)

		#Parse and Form Raw File Splicing Commands

		splicer_run_commands = []

		for B in np.arange(len(B_vals)):	
			splicer_run = 'python ' + '../extractor/raw_splicer.py ' + '../pipeline_playground/rawfiles/' \
			+ str(start_times[B]) + '_' + str(end_times[B]) + '/ ' + '9 ' + 'spliced' +  files[0][33:-25]
			splicer_run_commands.append(splicer_run)

		cepoch = 58178
		output = files[0][57:-25]
		polar = 4
		phasebin = 32768
		p = 0
		chan = 3712
		samples = 558 #number of MB

		#Parse and Form Coherent Dispersion Commands

		cdd_run_commands = []

		for B in np.arange(len(B_vals)):
			cdd_run = 'dspsr ' + '-T ' + str(time_widths[B]) + ' -c ' + str(time_widths[B]) \
			+ ' -cepoch ' + str(cepoch) + ' -O ' + str(output) + '_' + str(start_times[B]) + '_' \
			+ str(end_times[B]) + ' -D ' + str(DMs[B]) + ' -d ' + str(polar) + ' -b ' \
			+ str(phasebin) + ' -K -F ' + str(chan) + ':D ' + '-p 0 ' + '-U ' + str(samples) \
			+ ' ../pipeline_playground/rawfiles/' + str(start_times[B]) + '_' + str(end_times[B]) + ' spliced' +  files[0][33:-25]
			cdd_run_commands.append(cdd_run)

	return extract_run_commands, splicer_run_commands, cdd_run_commands



def main():

	extract_run_commands, splicer_run_commands, cdd_run_commands = extractncdd(filpaths, rawpaths, fieldnames, csvs)

	#Extract Raw Voltages

	for erc in extract_run_commands:
		os.system(erc)

	#Splice Raw Files Into Contiguous Raw File

	for src in splicer_run_commands:
		os.system(src)
	
	#Coherently Dedisperse Raw File With DSPSR

	for cdd in cdd_run_commands:
		os.system(cdd)

main()









