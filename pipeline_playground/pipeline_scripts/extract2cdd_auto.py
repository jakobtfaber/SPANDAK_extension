import numpy as np
import matplotlib.pyplot as plt
import fnmatch
import sys
import os
import pandas as pd
import csv
sys.path.insert(0, os.path.abspath('../extractor'))

def read_data(database, csv):
	
	#Read in -tentative 'database'- .csv

	filepaths = pd.read_csv(str(database))
	filpaths = filepaths.iloc[:,0]
	rawpaths = filepaths.iloc[1, :][1:]
	fieldnames = filepaths.columns[1:]
	#csvs = ['spliced_guppi_57991_49905_DIAG_FRB121102_0011.gpuspec.0001.8.4chan.csv']
	#csvs = ['/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/R3_csvs/AG_FRB180916_0005.rawspec.0001.csv']
	csvs = [str(csv)]
	
	#for csv in os.listdir(str(csv_dir)):
	#	if csv.endswith('.csv'):
	#		csvs.append(csv)

	return filepaths, filpaths, rawpaths, fieldnames, csvs


def parse_spandak(csvs):

	for csv in csvs:
		
		#Read in .csv produced by SPANDAK (on FRB121102 as an example)
		cands = pd.read_csv(csv)

		#Isolate B-valued candidates (probable if not using ML functionality)
		B_idx_nofil = cands[cands.iloc[:, :]['Category']=='B'].index.values

		#Find DMs for B candidates
		DMs_nofil = [[i for i in cands.loc[:, :]['DM']][b] for b in \
					B_idx_nofil]

		#Find DM-thresholded B-valued candidates
		B_idx = [B_idx_nofil[i] for i in np.arange(len(B_idx_nofil)) if 500 < DMs_nofil[i] < 700]

		#Find DMs for DM-thresholded B-valued candidates
		DMs = [[i for i in cands.loc[:, :]['DM']][b] for b in \
					B_idx]

		#Find filenames for B candidtaes
		files = [[i for i in cands.loc[:, :]['filename']][b] for b in \
					B_idx]

		#Find source name
		sourcename = [[i for i in cands.loc[:, :]['SourceName']][b] \
						for b in B_idx]

		#Find time constraints
		#Time widths
		time_widths = [[i for i in cands.loc[:, :]['WIDTH']][b] for b in \
						B_idx]

		#Find ime centers
		parse_center = [[i for i in [j.split('_') for j in \
		cands.iloc[:, :]['PNGFILE']]][b][2] for b in B_idx]
		time_stamps = [np.float(m[:-3]) for m in parse_center]

		#Calculate dispersion delay in s
		band = [[i for i in cands.loc[:, :]['BANDWIDTH']][b] for b in \
					B_idx]
		#band = 400
		tau_disp = []
		for B in np.arange(len(B_idx)):	
			tau_disp.append((4.149e+3) * DMs[B] * ((band[B]*4)**(-2)))

		#Calculate start and end times for raw voltage extraction

		#Time limits - widths adjusted to 3 * dispersion delay (i.e., tau_disp)
		#start_times = [time_stamps[b]-time_widths[b] for b in np.arange(len(B_idx))]
		#start_times = [time_stamps[b]-(2 * tau_disp[b]) for b in np.arange(len(B_idx))]
		start_times = [time_stamps[b]-(2*tau_disp[b]) for b in \
						np.arange(len(B_idx))]
		#end_times = [time_stamps[b]+time_widths[b] for b in np.arange(len(B_idx))]
		#end_times = [time_stamps[b]+(2 * tau_disp[b]) for b in np.arange(len(B_idx))]
		end_times = [time_stamps[b]+(50*tau_disp[b]) for b in \
					np.arange(len(B_idx))]

		return B_idx, files, DMs, sourcename, time_widths, time_stamps, \
				start_times, end_times, tau_disp, csv

def extract_auto(rawpaths, fieldnames, B_idx, files, filepaths, \
				start_times, end_times, csv):

	#Parse and Form Raw Voltage Extraction Commands

	#Store extraction commands in list
	extract_run_commands = []

	#extract_run_commands = {}
	#extract_run_commands['0'] = []
	#extract_run_commands['1'] = []
	#extract_run_commands['2'] = []
	#extract_run_commands['3'] = []
	#extract_run_commands['all'] = []

	'''Group together candidates with start times within a standard... 
	deviation of 0.2...works for now, can change..'''

	st_sorted = np.sort(start_times)
	gaps = [y - x for x, y in zip(st_sorted[:-1], st_sorted[1:])]
	sd = np.std(gaps)
	group_starts = [[[start_times.index(st_sorted[0]), st_sorted[0]]]]
	for t in st_sorted[1:]:
		if (t - group_starts[-1][-1][-1]) > 0.2:
			group_starts.append([])
		group_starts[-1].append([start_times.index(t), t])

	'''Store candidate idxs into subcategories within dictionary depending
	based on the subbands in which they were flagged by SPANDAK'''

	sub_cands = {}
	sub_cands['cand_0'] = []
	sub_cands['cand_1'] = []
	sub_cands['cand_2'] = []
	sub_cands['cand_3'] = []
	sub_cands['all'] = []
        
#
	for B in np.arange(len(B_idx)):
		if 'diced_0' in files[B]:
			sub_cands['cand_0'].append(B)
		elif 'diced_1' in files[B]:
			sub_cands['cand_1'].append(B)
		elif 'diced_2' in files[B]:
			sub_cands['cand_2'].append(B)
		elif 'diced_3' in files[B]:
			sub_cands['cand_3'].append(B)
		else:
			sub_cands['all'].append(B)

	plot_bands = {}
	plot_bands['3.8_5.1'] = []
	plot_bands['3.8_6.4'] = []
	plot_bands['3.8_7.7'] = []
	plot_bands['5.3_6.4'] = []
	plot_bands['5.3_7.7'] = []
	plot_bands['5.3_9'] = []
	plot_bands['6.6_7.7'] = []
	plot_bands['6.6_9'] = []
	plot_bands['7.9_9'] = []
	plot_bands['3.8_9'] = []

	for cgroup in group_starts:
		if len(cgroup) == 3:
			if cgroup[0][0] in sub_cands['all']:
				plot_bands['3.8_9'].append([cgroup[0][0], \
					start_times[cgroup[0][0]], end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3'] and cgroup[1][0] \
				in sub_cands['cand_2'] and \
					cgroup[2][0] in sub_cands['cand_1']:
				plot_bands['3.8_7.7'].append([[cgroup[0][0], cgroup[1][0], \
					cgroup[2][0]], start_times[cgroup[2][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3'] and cgroup[1][0] in \
					sub_cands['cand_2'] and \
					cgroup[2][0] in sub_cands['cand_0']:
				plot_bands['3.8_9'].append([[cgroup[0][0], cgroup[1][0], \
					cgroup[2][0]], start_times[cgroup[2][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3'] and cgroup[1][0] in \
					sub_cands['cand_1'] and \
					cgroup[2][0] in sub_cands['cand_0']:
				plot_bands['3.8_9'].append([[cgroup[0][0], cgroup[1][0], \
					cgroup[2][0]], start_times[cgroup[2][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_2'] and cgroup[1][0] in \
					sub_cands['cand_1'] and \
					cgroup[2][0] in sub_cands['cand_0']:
				plot_bands['5.3_9'].append([[cgroup[0][0], cgroup[1][0], \
					cgroup[2][0]], start_times[cgroup[2][0]], \
					end_times[cgroup[0][0]]])
		if len(cgroup) == 2:
			if cgroup[0][0] in sub_cands['all']:
				plot_bands['3.8_9'].append([cgroup[0][0], \
					start_times[cgroup[0][0]], end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3'] and cgroup[1][0] in \
				sub_cands['cand_2']:
				plot_bands['3.8_6.4'].append([[cgroup[0][0], cgroup[1][0]], \
					start_times[cgroup[1][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3'] and cgroup[1][0] in \
				sub_cands['cand_1']:
				plot_bands['3.8_7.7'].append([[cgroup[0][0], cgroup[1][0]], \
					start_times[cgroup[1][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3'] and cgroup[1][0] in \
				sub_cands['cand_0']:
				plot_bands['3.8_9'].append([[cgroup[0][0], cgroup[1][0]], \
					start_times[cgroup[1][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_2'] and cgroup[1][0] in \
				sub_cands['cand_1']:
				plot_bands['5.3_7.7'].append([[cgroup[0][0], cgroup[1][0]], \
					start_times[cgroup[1][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_2'] and cgroup[1][0] in \
				sub_cands['cand_0']:
				plot_bands['5.3_9'].append([[cgroup[0][0], cgroup[1][0]], \
					start_times[cgroup[1][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_1'] and cgroup[1][0] in \
				sub_cands['cand_0']:
				plot_bands['6.6_9'].append([[cgroup[0][0], cgroup[1][0]], \
					start_times[cgroup[1][0]], \
					end_times[cgroup[0][0]]])
		else:
			if cgroup[0][0] in sub_cands['all']:
				plot_bands['3.8_9'].append([cgroup[0][0], start_times[cgroup[0][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_3']:
				plot_bands['3.8_5.1'].append([cgroup[0][0], start_times[cgroup[0][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_2']:
				plot_bands['5.3_6.4'].append([cgroup[0][0], start_times[cgroup[0][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_1']:
				plot_bands['6.6_7.7'].append([cgroup[0][0], start_times[cgroup[0][0]], \
					end_times[cgroup[0][0]]])
			elif cgroup[0][0] in sub_cands['cand_0']:
				plot_bands['7.9_9'].append([cgroup[0][0], start_times[cgroup[0][0]], \
					end_times[cgroup[0][0]]])
#   
	##for B in np.arange(len(B_idx)):
	#for B in sub_cands['cand_0']:
	#	#if 'diced_0' in files[B]:
	#	for raw0 in np.arange(0, 7):
	#		extract_run_0 = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/extract_blocks.py ' \
	#		+ rawpaths[raw0] + ' ' + 'blc' + str(fieldnames[raw0][3:]) + files[B][75:-25] + ' ' \
	#		+ str(start_times[B]) + ' ' + str(end_times[B]) \
	#		+ ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/SPANDAK_121102_raws/' \
	#		+ str(start_times[B]) + '_' + str(end_times[B]) + '_7.9_9/'
	#		extract_run_commands['0'].append(extract_run_0)
	#for B in sub_cands['cand_1']:
	##elif 'diced_1' in files[B]:
	#	for raw1 in np.arange(7, 14):
	#		extract_run_1 = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/extract_blocks.py ' \
	#		+ rawpaths[raw1] + ' ' + 'blc' + str(fieldnames[raw1][3:]) + files[B][75:-25] + ' ' \
	#		+ str(start_times[B]) + ' ' + str(end_times[B]) \
	#		+ ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/SPANDAK_121102_raws/' \
	#		+ str(start_times[B]) + '_' + str(end_times[B]) + '_6.6_7.7/'
	#		extract_run_commands['1'].append(extract_run_1)
	#for B in sub_cands['cand_2']:
	##elif 'diced_2' in files[B]:
	#	for raw2 in np.arange(14, 21):
	#		extract_run_2 = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/extract_blocks.py ' \
	#		+ rawpaths[raw2] + ' ' + 'blc' + str(fieldnames[raw2][3:]) + files[B][75:-25] + ' ' \
	#		+ str(start_times[B]) + ' ' + str(end_times[B]) \
	#		+ ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/SPANDAK_121102_raws/' \
	#		+ str(start_times[B]) + '_' + str(end_times[B]) + '_5.3_6.4/'
	#		extract_run_commands['2'].append(extract_run_2)
	#for B in sub_cands['cand_3']:
	##elif 'diced_3' in files[B]:
	#	for raw3 in np.arange(21, 29):
	#		extract_run_3 = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/extract_blocks.py ' \
	#		+ rawpaths[raw3] + ' ' + 'blc' + str(fieldnames[raw3][3:]) + files[B][75:-25] + ' ' \
	#		+ str(start_times[B]) + ' ' + str(end_times[B]) \
	#		+ ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/SPANDAK_121102_raws/' \
	#		+ str(start_times[B]) + '_' + str(end_times[B]) + '_3.8_5.1/'
	#		extract_run_commands['3'].append(extract_run_3)
	for B in sub_cands['all']:
	#for B in np.arange(len(B_idx)):
		for raw in np.arange(len(rawpaths)):
			extract_run = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/extract_blocks.py ' \
			+ rawpaths[raw] + ' ' + 'blc' + str(fieldnames[raw][3:]) + files[B][33:-25] + ' '  \
			+ str(start_times[B]) + ' ' + str(end_times[B]) \
			+ ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/' +  str(start_times[B]) + '_' \
			+ str(end_times[B]) + '_3.8_9/'
			extract_run_commands.append(extract_run)
	
	#for B in sub_cands['all']:
	#	for raw in np.arange(len(rawpaths)):
	#		extract_run = 'python ' + 'datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/extract_blocks.py ' \
	#		+ rawpaths[raw] + ' ' + 'blc7' + str(fieldnames[raw][4:]) + filepaths.iloc[:,0][1][5:-22] + ' ' \
	#		+ str(start_times[B]) + ' ' + str(end_times[B]) \
	#		+ ' /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/0003_raws/' \
	#		+  str(start_times[B]) + '_' + str(end_times[B]) + '_3.8_9/'
	#		extract_run_commands.append(extract_run)

	return extract_run_commands, sub_cands, plot_bands


def splice_auto(sub_cands, files, start_times, end_times, ex_raws_path):

	#Parse and Form Raw File Splicing Commands

	splicer_run_commands = []

	#121102
	l = 33
	r = -25
	#r3
	#l
	#r
	

	for B in sub_cands['all']:	
		splicer_run = 'python ' + '/datax/scratch/jfaber/SPANDAK_extension/extractor/splicer_raw.py ' \
		+  str(ex_raws_path) + str(start_times[B]) + '_' + str(end_times[B]) \
		+ '_3.8_9/ ' + '2 ' + 'spliced' +  files[B][l:r] + str(start_times[B]) + '_' + str(end_times[B]) + '.raw'
		splicer_run_commands.append(splicer_run)

	return splicer_run_commands

def gen_par(sourcename, B_idx, DMs):

	par_file = []
	source_int = [str(sub.split('_')[1])[3:] for sub in sourcename]

	#121102 rot freq
	F0 = 47.680639719963075
	#r3 rot freq
	#F0 = 12

	#str([int(s) for s in source_int[B].split() if s.isdigit])
	for B in np.arange(len(B_idx)):
		par_txt = 'PSR  ' + str(source_int[B]) + '\n' \
		+ 'RAJ  ' + '05:31:58.70' '\n' + 'DECJ  ' + '+33:08:52.5' + '\n' \
		+ 'C The following is (1+v/c)*f_topo where f_topo is 1/(2048*dt)' + '\n' \
		+ 'F0  ' + str(F0) + '\n' + 'DM  ' + str(DMs[B]) + '\n' + 'PEPOCH  ' \
		+ '57962.373622685185186' + '\n' + 'CLOCK  ' + 'UNCORR'
		par_file.append(par_txt)
		
	return source_int, par_file


def cdd_auto(sub_cands, files, par_fil_paths, start_times, end_times):

	#Specify CDD parameters

	cepoch = 58178
	#121102
	output = files[0][57:-25]
	#r3
	#output = files[0][...]
	polar = 4
	phasebin = 2048
	p = 0
	chan = 6080
	samples = 1024 #number of MB

	#121102
	l = 33
	r = -25

	#Parse and Form Coherent Dispersion Commands

	cdd_run_commands = []

	for B in sub_cands['all']:
		cdd_run = 'dspsr ' + '-U ' + str(samples) + ' -F ' + str(chan) + ':D ' \
		+ ' -K ' + ' -d ' + str(polar) + ' -b  ' + str(phasebin) + ' -E ' \
		+ par_fil_paths[B] + ' -s -a psrfits -e fits ' \
		+ '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/bursts/' \
		+ str(start_times[B]) + '_' + str(end_times[B]) + '_3.8_9/raws/' + 'spliced' \
		+  files[B][l:r] + str(start_times[B]) + '_' + str(end_times[B]) + '.raw'
		cdd_run_commands.append(cdd_run)

	return cdd_run_commands

if __name__ == "__main__":

	database = sys.argv[1]
	csv = sys.argv[2]
	#database="database_r3.csv"
	
	filepaths, filpaths, rawpaths, fieldnames, csvs = read_data(database, csv)
	#print('Filepaths', filepaths)
	#print('Filpaths ', filepaths.iloc[:,0][1])
	#print('Fils ', filpaths)
	B_idx, files, DMs, sourcename, time_widths, time_stamps, start_times, end_times, \
		tau_disp, csv = parse_spandak(csvs)
	#print('B Index ', B_idx)
	#print('Start times ', start_times)
#	print('Tau ', tau_disp)
	#print('Start time ', start_times)
	#print('End time ', end_times)
	ex_raws_path = '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/'
	extract_run_commands, sub_cands, plot_bands = extract_auto(rawpaths, \
		fieldnames, B_idx, files, filepaths, start_times, end_times, csv)
	splicer_run_commands = splice_auto(sub_cands, files, start_times, end_times, ex_raws_path)
	source_int, par_file = gen_par(sourcename, B_idx, DMs)
	
	#for b in sub_cands['all']:
	#	print(time_stamps[b])
	#print('Sub cands', sub_cands)
	#print('Plot Bands', plot_bands)
	#print('Start Times', start_times)
	#print('Par Files: ', par_file)
	#print('Dispersion Delay: ', tau_disp)
	#print('Start Time: ', start_times)
	#print('End Time: ', end_times)
	#print(extract_run_commands[0])
	#print("Raw Voltage Paths: ", rawpaths)

	#Extract Raw Voltages
	#for B in sub_cands['all']:	
#		os.system('mkdir /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/' + str(start_times[B]) + '_' + str(end_times[B]) + '_7.9_9')
#		os.system('mkdir /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/' + str(start_times[B]) + '_' + str(end_times[B]) + '_6.6_7.7')
#		os.system('mkdir /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/' + str(start_times[B]) + '_' + str(end_times[B]) + '_5.3_6.4')
#		os.system('mkdir /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/' + str(start_times[B]) + '_' + str(end_times[B]) + '_3.8_5.1')
		#os.system('mkdir /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/R3/' + str(csv[13:17]) + '/' + str(start_times[B]) + '_' + str(end_times[B]) + '_3.8_9')
		#os.system('mkdir /datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/' + str(start_times[B]) + '_' + str(end_times[B]) + '_3.8_9')
	#print(extract_run_commands[0])
	#for erc in extract_run_commands:		
	#for k,v in extract_run_commands.items():
		#print('Extract 1: ', extract_run_commands['1'])
#		for erc in extract_run_commands['1']:
		#print("Extraction Commands: ", erc)
		#os.system(erc)
#
	##Splice Raw Files Into Contiguous Raw File
#
	#for src in splicer_run_commands:
	#	print('Splice Raw Commands :', src)
		#os.system(src)

	par_fil_paths = []
	for B in np.arange(len(B_idx)):
		par_fil = '/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/pars/' + 'FRB_' + str(start_times[B]).split('.')[0] + '.par'
		#par_fil = '/Users/jakobfaber/Documents/spandak_extended/SPANDAK_extension/pipeline_playground/parfiles/' + 'FRB_' + str(source_int[B]) + '_' + str(B_idx[B]) + '.par'
		par_fil_paths.append(par_fil)
		#par = open(par_fil, "w")
		#par.write(par_file[B])
		#par.close()
	#print(par_fil_paths)
#
	cdd_run_commands = cdd_auto(sub_cands, files, par_fil_paths, start_times, end_times)

	#Coherently Dedisperse
	for cdd in cdd_run_commands:
		#os.system('/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/' + cdd.split('/')[14] + '/fits/ ' + cdd)
		#print('Coherent Dedisp Commands: ', cdd)
		os.chdir(r'/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/FRB121102/bursts/' + str(cdd.split('/')[16]) + '/fits')
		print('DSPSR Output Funnelling Into: ' + os.getcwd())
		os.system(cdd)
		print('Coherent Dedispersion Complete')

	#database = 'pac -wp . -u fits'
	#os.system(database)
	#print('Database created')
	#fluxcal = 'fluxcal -i 15 -d database -c fluxcal.cfg'
	#os.system(fluxcal)
	#print(Flux & Pol calibration initiated)
	#calib = 'pac -x -d database.txt ' + pulse_fits
	#os.system(calib)
	#print(Calibration complete)










