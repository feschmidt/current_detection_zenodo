# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

r"""
HEMT thermal noise as calibration standard

We can calibrate the powers and attenuations in our circuit by using the thermal noise temperature of the HEMT and a VNA reference measurement as calibration standards.
The HEMT noise power is given by
\begin{align}
P_{\rm noise,HEMT}=10\log_{10}\left(\frac{\k_{\rm B}T_{\rm e}}{\SI{1}{mW}}\right) + 10\log_{10}\left( \frac{\rm IFBW}{\SI{1}{Hz}} \right)
\end{align}
From our reference measurement, we can average over the (to voltage) linearized S21 measurement background far from resonance, as this value will give us the background signal-to-noise ratio.
"""

import stlabutils
import glob
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants
import pickle


def func_PdBmtoV(PindBm, Z0=50):
    return np.sqrt(2*Z0*10**((PindBm-30)/10))


plotall = True
overview_plot = True

# Reading in the raw data
myfile = sorted(glob.glob(
    'data_raw/F32_megameasurement3/F33A*/*.dat'))
myfile = myfile[0]
print(myfile)
data = stlabutils.readdata.readdat(myfile)
mymtx = stlabutils.framearr_to_mtx(
    data[:-1], 'S21dB (dB)', xkey='Frequency (Hz)', ykey='Is (A)')
mymtx.applystep('rotate_ccw')

wbval = (0.1, 0.1)
cmap = 'RdBu'
lims = np.percentile(mymtx.pmtx.values, (wbval[0], 100 - wbval[1]))
vmin = lims[0]
vmax = lims[1]
extents = mymtx.getextents()

# Plotting the raw data
# Full overview
plt.imshow(mymtx.pmtx, aspect='auto', cmap=cmap, extent=(
    extents[0]/1e-6, extents[1]/1e-6, extents[2]/1e9, extents[3]/1e9), vmin=vmin, vmax=vmax)
plt.xlim(0, 8)
cbar = plt.colorbar()
cbar.set_label('S11dB (dB)')
plt.ylabel('Frequency (GHz)')
plt.xlabel('Bias current (uA)')
plt.title('Full overview')
plt.tight_layout()
if plotall:
    plt.show()
plt.close()

# Zoom into the traces for SNR extraction
[plt.plot(data[ii]['Frequency (Hz)'], data[ii]['S21dB (dB)'])
 for ii in range(len(data))]
# plt.xlim(7.44e9,7.46e9)
# plt.ylim(-50.5,-48.5)
plt.xlim(7.34e9, 7.365e9)
plt.ylim(-51, -49)
plt.title('All raw data traces')
if plotall:
    plt.show()
plt.close()

# trimming data
freqs = [x['Frequency (Hz)'] for x in data]
# we lost 20dB at the directional coupler after the MITEQs. This doesn't change anything though.
S21corr = [x['S21dB (dB)']+20 for x in data]

# idx = [np.logical_and(x<=7.46e9,x>=7.44e9) for x in freqs]
idx = [np.logical_and(x <= 7.365e9, x >= 7.33e9) for x in freqs]

freqs = [x[y] for x, y in zip(freqs, idx)]
S21corr = [x[y] for x, y in zip(S21corr, idx)]

# flattening data
flat_freqs = np.array([item for sublist in freqs for item in sublist])
flat_S21corr = np.array([item for sublist in S21corr for item in sublist])
plt.plot(flat_freqs, flat_S21corr, '.')
plt.title('{} datapoints'.format(len(flat_S21corr)))
if plotall:
    plt.show()
plt.close()

# Linearizing to Volts
flat_S21corr = 10**(flat_S21corr/20)

# Statistics on raw data
sigma_out = np.std(flat_S21corr)
S_out = np.mean(flat_S21corr)

# removing background as check, and statistics
m, b = np.polyfit(flat_freqs, flat_S21corr, 1)
flat_S21corr_noslope = np.array(flat_S21corr)-(m*np.array(flat_freqs)+b)
sigma_out_noslope = np.std(flat_S21corr_noslope)
S_out_noslope = np.mean(flat_S21corr_noslope)


if overview_plot:
    # Overview plot
    fig, (ax1, ax2) = plt.subplots(2, 2, figsize=(12, 8))

    plt.sca(ax1[0])
    plt.plot(flat_freqs, flat_S21corr, '.')
    # plt.plot(flat_freqs,m*np.array(flat_freqs)+b,lw=2)
    plt.axhline(S_out, c='k')
    plt.axhline(S_out+sigma_out, c='C1', ls='--', lw=2)
    plt.axhline(S_out-sigma_out, c='C1', ls='--', lw=2)
    # plt.title('Raw data with fit')
    plt.title(r'$P_{\rm out}='+r'({:.2e}\pm{:.2e})$ W, npoints={}'.format(S_out,
                                                                          sigma_out, len(flat_S21corr)))

    plt.sca(ax1[1])
    plt.plot(flat_freqs, flat_S21corr_noslope, '.')
    # plt.plot(flat_freqs,[0]*len(flat_freqs),lw=2)
    plt.axhline(S_out_noslope, c='k')
    plt.axhline(S_out_noslope+sigma_out_noslope, c='C1', ls='--', lw=2)
    plt.axhline(S_out_noslope-sigma_out_noslope, c='C1', ls='--', lw=2)
    # plt.title('Slope removed shows flat line')
    plt.title(r'$P_{\rm out}='+r'({:.2e}\pm{:.2e})$ W, npoints={}'.format(
        S_out_noslope, sigma_out_noslope, len(flat_S21corr_noslope)))

    plt.sca(ax2[0])
    n, bins, patches = plt.hist(flat_S21corr, 100)
    plt.axvline(S_out, c='k', ls='-', lw=2)
    plt.axvline(S_out+sigma_out, c='C1', ls='--', lw=2)
    plt.axvline(S_out-sigma_out, c='C1', ls='--', lw=2)
    plt.title('Histogram over raw data')

    plt.sca(ax2[1])
    n, bins, patches = plt.hist(flat_S21corr_noslope, 100)
    plt.axvline(S_out_noslope, c='k', ls='-', lw=2)
    plt.axvline(S_out_noslope+sigma_out_noslope, c='C1', ls='--', lw=2)
    plt.axvline(S_out_noslope-sigma_out_noslope, c='C1', ls='--', lw=2)
    plt.title('Histogram over flat line data')

    plt.tight_layout()
    plt.show()
    plt.close()

ifbw = 1000
navg = 1

# here we're using the standard deviation of the signal with removed background
SNRout = abs(S_out/(sigma_out_noslope*np.sqrt(navg)))
print("SNR in linear units: ", SNRout)
SNRout_dB = 20*np.log10(SNRout)
print("SNR in dB: ", SNRout_dB)
kB = scipy.constants.k
Te = 2
sigma_in_per_df = kB*Te

P_noise_HEMT = 10*np.log10(sigma_in_per_df/1e-3)+10*np.log10(ifbw/1)
print("Noise power HEMT (dBm): ", P_noise_HEMT)

P_in_HEMT = P_noise_HEMT+SNRout_dB
print("Signal power arriving at HEMT (dBm): ", P_in_HEMT)

atten_sample_HEMT = 2
print("Attenuation between sample and HEMT (dB): ", atten_sample_HEMT)

P_in_sample = P_in_HEMT+atten_sample_HEMT
print("Power arriving at sample (dBm): ", P_in_sample)
V_in_sample = func_PdBmtoV(P_in_sample)
print("Voltage arriving at sample (V): ", V_in_sample)

Pout_VNA = data[0]['Power (dBm)'][0]
print("Output power VNA (dBm): ", Pout_VNA)

total_input_attenuation = P_in_HEMT - Pout_VNA
print("Total input attenuation VNA (dB): ", total_input_attenuation)

# Note that the signal generator has 65dB less attenuation on its input line
diff_attenuation = 65
total_input_attenuation_SG = total_input_attenuation + diff_attenuation
print("Total input attenuation SG (dB): ", total_input_attenuation_SG)

# Exporting for later use
pkldict = {
    'Pout_VNA (dBm)': Pout_VNA,
    'SNRout_dB (dB)': SNRout_dB,
    'P_noise_HEMT (dBm)': P_noise_HEMT,
    'P_in_sample (dBm)': P_in_sample,
    'V_in_sample (V)': V_in_sample,
    'atten_sample_HEMT (dB)': atten_sample_HEMT,
    'P_in_HEMT (dBm)': P_in_HEMT,
    'input_att_VNA (dB)': total_input_attenuation,
    'input_att_SG (dB)': total_input_attenuation_SG,
    'diff_attenuation (dB)': diff_attenuation
}
exportpath = "data_processed/sensitivity_estimation/power_in.pkl"
print('Exporting the following data to', exportpath)
print(pkldict)
check = input('Do you want to export this data? [y/n]')
if check == 'y':
    pickle.dump(pkldict, open(
        exportpath, "wb"))
    print('Data exported')
elif check == 'n':
    print('Not exporting')
else:
    raise ValueError("Invalid reply! Has to be 'y' or 'n'!")
