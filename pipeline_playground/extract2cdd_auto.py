import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
import csv
sys.path.insert(0, os.path.abspath('../extractor'))

def read_data(database="database.csv"):
	
	#Read in -tentative 'database'- .csv
	#csv_dir = "/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/SPANDAK_121102_csvs"
	csv_dir = "/Users/jakobfaber/Documents/spandak_extended/SPANDAK_extension/pipeline_playground/SPANDAK_121102_csvs/"
	filepaths = pd.read_csv(database)
	filpaths = filepaths.iloc[:,0]
	rawpaths = filepaths.iloc[1, :][1:]
	fieldnames = filepaths.columns[1:]
	csvs = ['spliced_guppi_57991_49905_DIAG_FRB121102_0011.gpuspec.0001.8.4chan.csv']

	#for csv in os.listdir(str(csv_dir)):
	#	if csv.endswith("49905_DIAG_FRB121102_0011.gpuspec.0001.8.4chan.csv"):
	#		csvs.append(csv)

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

		#Find source name
		sourcename = [[i for i in cands.loc[:, :]['SourceName']][b] for b in B_idx]

		#Find time constraints
		#Time widths
		time_widths = [[i for i in cands.loc[:, :]['WIDTH']][b] for b in B_idx]

		#Time centers
		parse_center = [[i for i in [j.split('_') for j in cands.iloc[:, :]['PNGFILE']]][b][2] for b in B_idx]
		time_stamps = [np.float(m[:-3]) for m in parse_center]

		#Calculate dispersion delay in s
		band = [[i for i in cands.loc[:, :]['BANDWIDTH']][b] for b in B_idx]
		tau_disp = []
		for B in np.arange(len(B_idx)):	
			tau_disp.append((4.149e+3) * DMs[B] * (band[B]**(-2)))

		#Time limits - widths adjusted to 3 * dispersion delay (i.e., tau_disp)
		#start_times = [time_stamps[b]-time_widths[b] for b in np.arange(len(B_idx))]
		start_times = [time_stamps[b]-(3 * tau_disp[b]) for b in np.arange(len(B_idx))]
		#end_times = [time_stamps[b]+time_widths[b] for b in np.arange(len(B_idx))]
		end_times = [time_stamps[b]+(3 * tau_disp[b]) for b in np.arange(len(B_idx))]

		return B_idx, files, DMs, sourcename, time_widths, start_times, end_times, tau_disp

def extract_auto(rawpaths, fieldnames, B_idx, files, start_times, end_times):

	#Parse and Form Raw Voltage Extraction Commands

	extract_run_commands = []

	for B in np.arange(len(B_idx)):
		for raw in np.arange(len(rawpaths)):
			extract_run = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/extract_blocks.py ' + rawpaths[raw] \
			 + ' ' + 'blc' + str(fieldnames[raw][3:]) + files[B][33:-25] + ' ' \
			 + str(start_times[B]) + ' ' + str(end_times[B]) \
			 + ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/SPANDAK_121102_raws/' \
			 + str(start_times[B]) + '_' + str(end_times[B])
			extract_run_commands.append(extract_run)

	return extract_run_commands


def splice_auto(B_idx, files, start_times, end_times):

	#Parse and Form Raw File Splicing Commands

	splicer_run_commands = []

	for B in np.arange(len(B_idx)):	
		splicer_run = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/splicer_raw.py ' \
		+ 'datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/SPANDAK_121102_raws/' \
		+ str(start_times[B]) + '_' + str(end_times[B]) + '_no7' '/ ' + '2 ' + 'spliced' +  files[0][33:-25]
		splicer_run_commands.append(splicer_run)

	return splicer_run_commands

def gen_par(sourcename, B_idx, DMs):

	par_file = []
	source_int = [str(sub.split('_')[1])[3:] for sub in sourcename]
	#str([int(s) for s in source_int[B].split() if s.isdigit])
	for B in np.arange(len(B_idx)):
		par_txt = 'PSR  ' + str(source_int[B]) + '\n' \
		+ 'RAJ  ' + '05:31:58.70' '\n' + 'DECJ  ' + '+33:08:52.5' + '\n' \
		+ 'C The following is (1+v/c)*f_topo where f_topo is 1/(2048*dt)' + '\n' \
		+ 'F0  ' + '47.680639719963075' + '\n' + 'DM  ' + str(DMs[B]) + '\n' + 'PEPOCH  ' \
		+ '57962.373622685185186' + '\n' + 'CLOCK  ' + 'UNCORR'
		par_file.append(par_txt)
		
	return source_int, par_file


def cdd_fits_auto(B_idx, files, par_fil_paths):

	#Specify CDD parameters

	cepoch = 58178
	output = files[0][57:-25]
	polar = 4
	phasebin = 2048
	p = 0
	chan = 5376
	samples = 1024 #number of MB

	#Parse and Form Coherent Dispersion Commands

	cdd_run_commands = []

	for B in np.arange(len(B_idx)):
		cdd_run = 'dspsr ' + '-U ' + str(samples) + ' -F ' + str(chan) + ':D ' \
		+ ' -K ' + ' -d ' + str(polar) + ' -b  ' + str(phasebin) + ' -E ' \
		+ par_fil_paths[B] + ' -s -a psrfits -e fits ' + 'datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/SPANDAK_121102_raws/' \
		+ str(start_times[B]) + '_' + str(end_times[B]) + '_no7/' + 'spliced' +  files[0][33:-25]
		cdd_run_commands.append(cdd_run)

	return cdd_run_commands


def main():
	
	filepaths, filpaths, rawpaths, fieldnames, csvs = read_data()
	B_idx, files, DMs, sourcename, time_widths, start_times, end_times, tau_disp = parse_spandak(csvs)
	extract_run_commands = extract_auto(rawpaths, fieldnames, B_idx, files, start_times, end_times)
	splicer_run_commands = splice_auto(B_idx, files, start_times, end_times)
	source_int, par_file = gen_par(sourcename, B_idx, DMs)
	
	#print('Par Files: ', par_file)
	#print('Dispersion Delay: ', tau_disp)
	#print('Start Time: ', start_times)
	#print('End Time: ', end_times)
	#print(extract_run_commands[0])
	#Extract Raw Voltages

	for erc in extract_run_commands:
		print('Extract Raw Commands: ', erc)
		#os.system(erc)
#
	##Splice Raw Files Into Contiguous Raw File
#
	for src in splicer_run_commands:
		print('Splice Raw Commands :', src)
		#os.system(src)

	par_fil_paths = []
	for B in np.arange(len(B_idx)):
		par_fil = 'datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/parfiles/' + 'FRB_' + str(source_int[B]) + '_' + str(B_idx[B]) + '.par'
		#par_fil = '/Users/jakobfaber/Documents/spandak_extended/SPANDAK_extension/pipeline_playground/parfiles/' + 'FRB_' + str(source_int[B]) + '_' + str(B_idx[B]) + '.par'
		par_fil_paths.append(par_fil)
		par = open(par_fil, "w")
		par.write(par_file[B])
		#par.close()


	cdd_run_commands = cdd_fits_auto(B_idx, files, par_fil_paths)

	#
	##Coherently Dedisperse Raw File With DSPSR
#
	for cdd in cdd_run_commands:
		print('Coherent Dedisp Commands: ', cdd)
		#os.system(cdd)

main()









