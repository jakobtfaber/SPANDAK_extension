import matplotlib
matplotlib.use('Agg')
import numpy as np
import psrchive
import pylab
import sys
import os
from tempfile import TemporaryFile

#sys.path.insert(0, os.path.abspath('/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/61.4627973333_67.0552026667_fits/'))

def main():
	#for fits in '*.fits':
	with open		
		npar_file = 'pulse_120390656' + '.npy'
		npar_file = TemporaryFile()
		arch = psrchive.Archive_load('/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/61.4627973333_67.0552026667_fits/pulse_120390656.fits')
		arch.dedisperse()
		arch.remove_baseline()
		arch.convert_state('Stokes')
		data = arch.get_data()
		np.save(data[:, 0, :, :].mean(0))

main()

