import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import pandas as pd
sys.path.insert(0, os.path.abspath('../extractor'))

#Read in -tentative 'database'- .csv
filepaths = pd.read_csv("database.csv")

#Read in .csv produced by SPANDAK (on FRB121102 as an example)
cands = pd.read_csv("SPANDAK_guppi_57991_51723_DIAG_FRB121102_0012.csv")
print(cands)

#Isolate B-valued candidates (probable without ML)
B_vals = cands[cands.iloc[:, :]['Category']=='B'].index.values
#print(B_vals)

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
print(start_times)
end_times = [time_stamps[i]+time_widths[i] for i in np.arange(len(B_vals))]
print(end_times)



