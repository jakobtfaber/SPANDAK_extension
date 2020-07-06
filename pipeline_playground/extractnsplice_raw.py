import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
import csv
sys.path.insert(0, os.path.abspath('../extractor'))

with open('database.csv', 'r') as infile:
	reader = csv.DictReader(infile)
	fieldnames = reader.fieldnames[1:]
#print(fieldnames[0][3:])

#Read in -tentative 'database'- .csv
filepaths = pd.read_csv("database.csv")
filpaths = filepaths.iloc[:,0]
rawpaths = filepaths.iloc[1, :][1:]
csvs = []

for csv in os.listdir("../pipeline_playground/"):
	if csv.endswith("csv.csv"):
		csvs.append(csv)

#print(csvs)
#print(rawpaths)

#print(filpaths)

#rawpaths = filepaths.iloc[:,1]

#print(rawpaths[0][1])

def run_extract(filpaths, rawpaths, fieldnames, csvs):

	for csv in csvs:
		#Read in .csv produced by SPANDAK (on FRB121102 as an example)
		cands = pd.read_csv(csv)
		#print(cands)

		#Isolate B-valued candidates (probable without ML)
		B_vals = cands[cands.iloc[:, :]['Category']=='B'].index.values
		#print(B_vals)

		#Find filenames
		files = [[i for i in cands.loc[:, :]['filename']][k] for k in B_vals]
		#Find time constraints
		#Widths
		time_widths = [[i for i in cands.loc[:, :]['WIDTH']][k] for k in B_vals]
		#print(time_widths)

		#Centers
		parse_center = [[i for i in [j.split('_') for j in cands.iloc[:, :]['PNGFILE']]][k][2] for k in B_vals]
		time_stamps = [np.float(m[:-3]) for m in parse_center]
		#print(time_stamps)

		#Limits
		start_times = [time_stamps[i]-time_widths[i] for i in np.arange(len(B_vals))]
		#print(start_times)
		end_times = [time_stamps[i]+time_widths[i] for i in np.arange(len(B_vals))]
		#print(end_times)

		extract_run_commands = []

		for B in np.arange(len(B_vals)):
			for raw in np.arange(len(rawpaths)):
				extract_run = 'python ' + '../extractor/extract_blocks.py ' + rawpaths[raw] + '' + 'blc' + str(fieldnames[raw][3:]) + files[B][33:-25] + ' ' + str(start_times[B]) + ' ' + str(end_times[B]) + ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/rawfiles'
				extract_run_commands.append(extract_run)


	return extract_run_commands



def main():

	extract_run_commands = run_extract(filpaths, rawpaths, fieldnames, csvs)
	print(extract_run_commands)
	#os.system(extract_run)

main()









