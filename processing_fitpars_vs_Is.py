"""
This file generates the pickle fitpars_vs_Is.pkl which contains the raw and fitted data of the reference measurement F33A*.dat for our megameasurement
"""

import numpy as np
import glob
import stlabutils
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from src.algo_peakheights import S11_fitting_full, calculate_dSdf, calculate_measurement, calculate_expected, calculate_nth_sideband
from src.plot_dnSdfn import plot_peakamp_sens_abs, plot_four_overview, plot3_2Dmaps, plotall_four_overview


myfile = sorted(glob.glob(
    'data_raw/F32_megameasurement3/F33A*/*.dat'))
myfile = myfile[0]
print(myfile)
data = stlabutils.readdata.readdat(myfile)
data = data[:-1]  # last line was aborted
print(data[0].head())

testplot_raw = True
plotfourovw = True
plotallovw = False
plot_unclipped = True
plot_clipped = True
finalplotovw = True

if testplot_raw:
    mymtx = stlabutils.framearr_to_mtx(
        data, 'S21dB (dB)', xkey='Frequency (Hz)', ykey='Is (A)')
    # mymtx = stlabutils.framearr_to_mtx(data,'S21dB (dB)',xkey='Frequency (Hz)',ykey='Iset (A)')
    mymtx.applystep('rotate_ccw')

    wbval = (0.1, 0.1)
    cmap = 'RdBu'
    lims = np.percentile(mymtx.pmtx.values, (wbval[0], 100 - wbval[1]))
    vmin = lims[0]
    vmax = lims[1]
    extents = mymtx.getextents()

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    plt.sca(ax1)
    plt.imshow(mymtx.pmtx, aspect='auto', cmap=cmap, extent=(
        extents[0]/1e-6, extents[1]/1e-6, extents[2]/1e9, extents[3]/1e9), vmin=vmin, vmax=vmax)
    plt.xlim(0, 8)
    cbar = plt.colorbar()
    cbar.set_label('S11dB (dB)')
    plt.ylabel('Frequency (GHz)')
    plt.xlabel('Bias current (uA)')

    plt.sca(ax2)
    [plt.plot(block['Frequency (Hz)']/1e9, block['S21dB (dB)'])
     for block in data]
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('S21dB (dB)')
    plt.title('Resonance shifts down for higher bias current')
    plt.tight_layout()
    plt.show()
    plt.close()

# Extracting resonance frequency and S11 line by line
S11res_vs_I, mydfs = S11_fitting_full(data)

# Clipping unreasonably large resonance frequencies (failed fits)
f0sraw = S11res_vs_I['f0 (Hz)']
israw = S11res_vs_I['Is (A)']
idx = np.where(np.array(f0sraw) <= 7.44e9)[0]
S11res_vs_I_clip = S11res_vs_I.iloc[idx]
mydfs_clip = [mydfs[x] for x in idx]

# doing a bunch of testplots in between
if plotfourovw:
    plot_four_overview(mydfs_clip[0], key='S11', unit='()')

if plotallovw:
    plotall_four_overview(mydfs_clip, key='S11', unit='()')

if plot_unclipped:
    plot3_2Dmaps(mydfs, 'S11abs ()', 'S11fitabs ()',
                 'S11fitnobgabs ()', supertitle='Unclipped data')

if plot_clipped:
    plot3_2Dmaps(mydfs_clip, 'S11abs ()', 'S11fitabs ()', 'S11fitnobgabs ()',
                 supertitle='Clipped data. Scale of x-axis is now incorrect (clipped data is simply not displayed)')

# Generating the dict to be exported
fitpars = pd.DataFrame({'Is (A)': abs(S11res_vs_I_clip['Is (A)']), 'f0 (Hz)': abs(S11res_vs_I_clip['f0 (Hz)']),
                        'Qint ()': [S11res_vs_I_clip.iloc[x]['fitpars']['Qint'].value for x in range(len(S11res_vs_I_clip))],
                        'Qext ()': [S11res_vs_I_clip.iloc[x]['fitpars']['Qext'].value for x in range(len(S11res_vs_I_clip))],
                        'theta (rad)': [S11res_vs_I_clip.iloc[x]['fitpars']['theta'].value for x in range(len(S11res_vs_I_clip))],
                        'df0 (Hz)': [S11res_vs_I_clip.iloc[x]['fitpars']['f0'].stderr for x in range(len(S11res_vs_I_clip))],
                        'dQint ()': [S11res_vs_I_clip.iloc[x]['fitpars']['Qint'].stderr for x in range(len(S11res_vs_I_clip))],
                        'dQext ()': [S11res_vs_I_clip.iloc[x]['fitpars']['Qext'].stderr for x in range(len(S11res_vs_I_clip))],
                        'dtheta (rad)': [S11res_vs_I_clip.iloc[x]['fitpars']['theta'].stderr for x in range(len(S11res_vs_I_clip))]})
print(fitpars.head())
fitpars = fitpars.set_index('Is (A)')
if finalplotovw:
    _, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
    plt.sca(ax1)
    plt.plot(fitpars.axes[0].values/1e-6, fitpars['f0 (Hz)'] /
             1e9, '.', label='f0 (GHz)')
    plt.legend()
    plt.sca(ax2)
    plt.plot(fitpars.axes[0].values/1e-6, fitpars['Qint ()'] /
             1e3, '.', label='Qint (1k)')
    plt.plot(fitpars.axes[0].values/1e-6, fitpars['Qext ()'] /
             1e3, '.', label='Qext (1k)')
    plt.legend()
    plt.sca(ax3)
    plt.plot(fitpars.axes[0].values/1e-6,
             fitpars['theta (rad)'], '.', label=r'$\theta$ (rad)')
    plt.legend()
    for ax in [ax1, ax2, ax3]:
        ax.set_xlabel('Is (uA)')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    # plt.savefig('plots/sensitivity_estimation/S11_fitresults.png',bbox_to_inches='tight')
    plt.show()
    plt.close()


# Exporting for later use
pickle.dump(mydfs_clip, open(
    "data_processed/sensitivity_estimation/S11_vs_Is.pkl", "wb"))
pickle.dump(S11res_vs_I_clip, open(
    "data_processed/sensitivity_estimation/S11res_vs_Is.pkl", "wb"))
pickle.dump(fitpars, open(
    "data_processed/sensitivity_estimation/fitpars_vs_Is.pkl", "wb"))
