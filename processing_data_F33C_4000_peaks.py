"""
This datafile lies at the core of part of the analysis:
It reads in the raw measured data recorded with the spectrum analyzer and extracts the individual peaks, calculates the peak heights and the frequencies at which they are located.
"""

import glob
import stlabutils
import matplotlib.pyplot as plt
import pickle
from src.class_peakfinding import Peakfinder
from src.plot_dnSdfn import Plot_raw_inputoutput

numcurrent = 4000
# measurement data
myfile = glob.glob(
    'data_raw/F32_megameasurement3/F33C*_vs_dF_{}/*.dat'.format(numcurrent))
myfile = myfile[0]
print(myfile)
data = stlabutils.readdata.readdat(myfile)
Iset = data[0]['Is (A)'][0]
# Plotting the raw data
Plot_raw_inputoutput(data)

# extracting peak heights and calculating sensitivities
myFinder = Peakfinder({'xkey': 'dF (Hz)'})
thedf = myFinder.GetFreqsAndPeaks(data)
thedf = thedf.set_index('dF (Hz)')

# exporting for later use
pickle.dump(thedf, open(
    "data_processed/sensitivity_estimation/data_F33C_{}_peaks.pkl".format(numcurrent), "wb"))
