import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd

#Read in -tentative 'database'- .csv
filepaths = pd.read_csv("database.csv")

filpaths = filepaths.iloc[:,0]

print(filpaths)

def run_spandak(filpaths, hidm=2000, lodm=100):

	run_commands = []
	
	for fil in filpaths:
<<<<<<< HEAD
		spandak_run = 'SPANDAK ' + '--fil ' + fil + ' --hidm ' + str(hidm) + ' --lodm ' + str(lodm) + ' --dorfi ' + '--logs=' + fil[:-4] + '.csv'
=======
		spandak_run = 'SPANDAK ' + '--fil ' + fil + ' --hidm ' + str(hidm) + ' --lodm ' + str(lodm) + ' --dorfi ' + '--logs=' + fil[26:-4] + '.csv'
>>>>>>> da66f6e97c873160cc72541b4b2760220a7f5e56
		
		run_commands.append(spandak_run)

	return run_commands


def main():
<<<<<<< HEAD
	
=======

>>>>>>> da66f6e97c873160cc72541b4b2760220a7f5e56
	run_commands = run_spandak(filpaths)

	for run in run_commands:
		#os.system(run)
		print(run)

main()
