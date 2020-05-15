"""
This datafile lies at the core of part of the analysis:
It reads in the raw measured data recorded with the spectrum analyzer and extracts the individual peaks, calculates the peak heights and the frequencies at which they are located.
"""
import os
import glob
import stlabutils
import matplotlib.pyplot as plt
import pickle
from src.class_peakfinding import Peakfinder
from src.plot_dnSdfn import Plot_raw_inputoutput

# numcurrent = 4000

for numcurrent in [0, 500, 1000, 2000, 3000, 4000, 5000, 6000, 6500, 6599, 6699, 6799, 6900, 7000, 7100, 7200, 7300]:
    exportfile = "data_processed/sensitivity_estimation/data_F33C_{}_peaks.pkl".format(
        round(numcurrent, -1))
    print(exportfile)
    if not os.path.isfile(exportfile):
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
        pickle.dump(thedf, open(exportfile, "wb"))
