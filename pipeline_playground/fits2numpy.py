import matplotlib
matplotlib.use('Agg')
import numpy as np
import psrchive
import pylab
import sys
import os

#sys.path.insert(0, os.path.abspath('/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/61.4627973333_67.0552026667_fits/'))

def fits2numpy():
	for fits in '/datax/scratch/jfaber/fits':
		npar_file = fits[:-4] + '.npy'
		with open(npar, 'wb') as npar_file:		
			arch = psrchive.Archive_load('/datax/scratch/jfaber/SPANDAK_extension/pipeline_playground/61.4627973333_67.0552026667_fits/pulse_120390656.fits')
			arch.dedisperse()
			arch.remove_baseline()
			arch.convert_state('Stokes')
			data = arch.get_data()
			np.save(npar_file, data[:, 0, :, :].mean(0))
			print('Array Written...')

fits2numpy()

