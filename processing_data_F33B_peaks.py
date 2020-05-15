# -*- coding: utf-8 -*-
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

"""
This datafile lies at the core of part of the analysis:
It reads in the raw measured data recorded with the spectrum analyzer and extracts the individual peaks, calculates the peak heights and the frequencies at which they are located.
"""

import glob
import stlabutils
import matplotlib.pyplot as plt
import pickle
from src.class_peakfinding import Peakfinder

myfile = sorted(glob.glob(
    'data_raw/F32_megameasurement3/F33B*/*.dat'))
myfile = myfile[0]
print(myfile)
data = stlabutils.readdata.readdat(myfile)
print(data[0].head())

for i, block in enumerate(data[::-4]):
    if i == 0:
        first_Is = block['Is (A)'].iloc[0]/1e-6
    plt.plot((block['Frequency (Hz)']-block['Frequency (Hz)'][len(block['Frequency (Hz)'])//2]-50*i)/1e3+0.9,
             block['Spectrum (dBm)']-50*i+450)  # , label=r'{:.1f}µA'.format(block['Is (A)'].iloc[0]/1e-6))
# plt.legend()
plt.title('Raw measured data, shifted for visibility. Iset from {}uA to {}uA'.format(
    block['Is (A)'].iloc[0]/1e-6, first_Is))
plt.tight_layout()
# plt.savefig('plots/data_F33B_peaks.png', bbox_to_inches='tight')
plt.show()
plt.close()

myFinder = Peakfinder({'xkey': 'Is (A)'})
thedf = myFinder.GetFreqsAndPeaks(data)

# Exporting for later use
pickle.dump(thedf, open(
    "data_processed/sensitivity_estimation/data_F33B_peaks.pkl", "wb"))
