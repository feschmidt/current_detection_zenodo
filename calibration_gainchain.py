# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

r"""
Calculating the gain chain from spectrum analyzer background

We can calibrate gain of our amplification chain (including cable losses) by comparing the power arriving at the sample to the noise floor of our spectrum analyzer.
"""

import stlabutils
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
import pickle
import pandas as pd

showall = True

# Loading in the data
myfolders = 'data_raw/F40_noisefloor_3/'
myfiles = sorted(glob.glob(myfolders+'1Hz_0dB/*/*.dat'))
allampsoff = [x for x in myfiles if 'allampsoff' in x][0]
MITEQsonHEMToff = [x for x in myfiles if 'HEMToffMITEQson' in x][0]
allampson = [x for x in myfiles if 'allampson' in x][0]

data_alloff = stlabutils.readdata.readdat(allampsoff)[0]
data_HEMToff = stlabutils.readdata.readdat(MITEQsonHEMToff)[0]
data_allon = stlabutils.readdata.readdat(allampson)[0]

data_alloff.head()

data_HEMToff.head()

data_allon.head()

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,4))
plt.sca(ax1)
plt.plot(data_alloff['Frequency (Hz)'],
         data_alloff['Spectrum (dBm)'], label='allampsoff')
plt.plot(data_HEMToff['Frequency (Hz)'],
         data_HEMToff['Spectrum (dBm)'], label='MITEQsonHEMToff')
plt.plot(data_allon['Frequency (Hz)'],
         data_allon['Spectrum (dBm)'], label='allampson')
plt.legend()
plt.ylabel('Spectrum (dBm)')
plt.title('Comparing the noise floor: We are HEMT limited')
plt.legend()
plt.sca(ax2)
diffdata = 10**(data_allon['Spectrum (dBm)']/10)/10**(data_HEMToff['Spectrum (dBm)']/10)
plt.plot(data_HEMToff['Frequency (Hz)'],
         10*np.log10(diffdata),
         label='allampson')
plt.ylabel('Power difference (dB)')
plt.title('Diffmean = {}dB'.format(10*np.log10(np.mean(diffdata))))
plt.tight_layout()
plt.show()
plt.close()

pkldump = {'allon':data_allon,'HEMToff':data_HEMToff}
pickle.dump(pkldump,open('data_final/SM_HEMT.pkl','wb'))

# The above plot shows us, that our measurement is HEMT-limited: Turning on the HEMT adds background noise.
# This means that the added noise from the HEMT is significantly larger than the added MITEQ noise.

noisefloor_dBm = data_allon['Spectrum (dBm)']
noisefloor_mW = 10**(noisefloor_dBm/10)
noisefloor_mW_avg = np.mean(noisefloor_mW)
noisefloor_dBm_avg = 10*np.log10(noisefloor_mW_avg)

print("Measurement noisefloor (dBm):", noisefloor_dBm_avg)

# Comparing this with the noise power of the HEMT
kB = scipy.constants.k
Te = 2
ifbw = data_allon['RBW (Hz)'][0]
P_noise_HEMT = 10*np.log10(kB*Te/1e-3)+10*np.log10(ifbw/1)
print("Noise power of HEMT (dBm):", P_noise_HEMT)

atten_sample_HEMT = pickle.load(
    open("data_processed/sensitivity_estimation/power_in.pkl", "rb"))['atten_sample_HEMT (dB)']
print("Attenuation between sample and HEMT (dB): ", atten_sample_HEMT)

gain = noisefloor_dBm_avg - P_noise_HEMT - atten_sample_HEMT
print("Gain of amplifier chain (dB):", gain)

pkldict = {'Noisefloor (dBm)': noisefloor_dBm_avg,
           'P_noise_HEMT (dBm)': P_noise_HEMT, 'Gain (dB)': gain}
exportpath = "data_processed/sensitivity_estimation/gaindata.pkl"
print('Exporting the following data to', exportpath)
print(pkldict)
check = input('OK? [y/n]')
if check == 'y':
    pickle.dump(pkldict, open(
        exportpath, "wb"))
    print('Data exported')
else:
    print('Not exporting')
