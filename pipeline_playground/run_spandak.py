import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd

#Read in -tentative 'database'- .csv
filepaths = pd.read_csv("database.csv")

filpaths = filepaths.iloc[:,0]

def run_spandak(filpaths, hidm=2000, lodm=100):

	run_commands = []
	
	for fil in filpaths:
		spandak_run = 'SPANDAK ' + '--fil ' + fil + ' --hidm ' + str(hidm) + ' --lodm ' + str(lodm) + ' --dorfi ' + '--logs=' + fil[:-4] + '.csv'
		
		run_commands.append(spandak_run)

	return run_commands


def main():
	
	run_commands = run_spandak(filpaths)

	for run in run_commands:
		os.system(run)

main()
