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
for numcurrent in [0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 6599, 6699, 6799, 6900, 7000, 7100, 7200, 7300, 7499, 7699]:
    exportfile = "data_processed/sensitivity_estimation/data_F33E_{}_peaks.pkl".format(
        round(numcurrent, -1))
    print(exportfile)
    if not os.path.isfile(exportfile):
        # measurement data
        myfile = glob.glob(
            'data_raw/F32_megameasurement3/F33E_*_SA_vs_Is_0to10uA_200mVpp_vs_pwr_{}/*.dat'.format(numcurrent))
        myfile = myfile[0]
        print(myfile)
        data = stlabutils.readdata.readdat(myfile)
        # Plotting the raw data
        Plot_raw_inputoutput(data, ykey='Carrier Power (dBm)',
                             xlabel='Carrier Power (dBm)', xscale=1)

        # extracting peak heights and calculating sensitivities
        myFinder = Peakfinder({'xkey': 'Carrier Power (dBm)'})
        thedf = myFinder.GetFreqsAndPeaks(data)
        thedf = thedf.set_index('Carrier Power (dBm)')

        # exporting for later use
        pickle.dump(thedf, open(exportfile, "wb"))
