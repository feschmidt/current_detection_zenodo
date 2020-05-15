import glob
import stlabutils
from src.algo_peakheights import S11_fitting_full
import numpy as np
import pandas as pd
import pickle

# alternatively, this is the data with the furthest achieved tuning
myfile = glob.glob(
    'data_raw/F30_50Hz_noise_debugging/F30_2019_03_22_14.09.44_VNA_vs_current_10uAalloffcleanracks/*.dat')
myfile = myfile[0]
myfile
data = stlabutils.readdata.readdat(myfile)
data[0].head()


S11res_vs_I, mydfs = S11_fitting_full(data[:-7], xvariable='Iset (A)')
# Clipping unreasonably large resonance frequencies (failed fits)
f0sraw = S11res_vs_I['f0 (Hz)']
israw = S11res_vs_I['Iset (A)']
idx = np.where(np.array(f0sraw) <= 7.44e9)[0]
S11res_vs_I_clip = S11res_vs_I.iloc[idx]
mydfs_clip = [mydfs[x] for x in idx]
# Generating the dict to be exported
fitpars = pd.DataFrame({'Is (A)': abs(S11res_vs_I_clip['Iset (A)']), 'f0 (Hz)': abs(S11res_vs_I_clip['f0 (Hz)']),
                        'Qint ()': [S11res_vs_I_clip.iloc[x]['fitpars']['Qint'].value for x in range(len(S11res_vs_I_clip))],
                        'Qext ()': [S11res_vs_I_clip.iloc[x]['fitpars']['Qext'].value for x in range(len(S11res_vs_I_clip))],
                        'theta (rad)': [S11res_vs_I_clip.iloc[x]['fitpars']['theta'].value for x in range(len(S11res_vs_I_clip))]})
f0_of_I = fitpars.set_index('Is (A)')
f0_of_I.head()

pickle.dump(f0_of_I,
            open("data_processed/sensitivity_estimation/fitpars_derivatives_vs_Is_F30.pkl", "wb"))
pickle.dump((S11res_vs_I_clip, mydfs_clip),
            open("data_processed/sensitivity_estimation/fitpars_derivatives_vs_Is_F30_fulldata.pkl", "wb"))
