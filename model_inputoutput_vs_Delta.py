# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

"""
This file calculates and plots the measured peak heights together with the ones calculated from input-output theory, as a function of fixed bias current and varying detuning Delta.
"""

# %%
import glob
import stlabutils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.constants
from scipy.interpolate import interp1d
import pickle
from src.algo_inputoutput import first_order_coeffs, second_order_coeffs, third_order_coeffs
from src.algo_peakheights import calculate_measurement
from src.plot_dnSdfn import Plot_model_inputoutput
from src.model_currentbias import f0, lossrate

# %%
fixcurrent = 4e-6
numcurrent = round(fixcurrent/1e-9)
saveall = True

# %%
# Loading in data and parameters
power_gain = pickle.load(
    open("data_processed/sensitivity_estimation/gaindata.pkl", "rb"))
gain = power_gain['Gain (dB)']

power_in = pickle.load(
    open("data_processed/sensitivity_estimation/power_in.pkl", "rb"))
V0 = power_in['V_in_sample (V)']

fitpars = pickle.load(
    open("data_processed/sensitivity_estimation/fitpars_vs_Is.pkl", "rb"))

# %%
# Setting up the model
# Load fitparameters for model
popt = pickle.load(open(
    'data_processed/sensitivity_estimation/popt_simulation_JJarrays.pkl', 'rb'))
Iset_meas = fitpars.axes[0]
Iset = np.linspace(min(Iset_meas), max(Iset_meas), 401)
f0_of_I = f0(Iset, *popt, nJJ=1)
G1_of_I = np.gradient(f0_of_I, Iset)
G2_of_I = np.gradient(G1_of_I, Iset)
G3_of_I = np.gradient(G2_of_I, Iset)

# and now interpolation
G1 = interp1d(Iset, G1_of_I)(fixcurrent)
G2 = interp1d(Iset, G2_of_I)(fixcurrent)
G3 = interp1d(Iset, G3_of_I)(fixcurrent)

# %%
# measurement data
meas_data = pickle.load(open(
    "data_processed/sensitivity_estimation/data_F33C_{}_peaks.pkl".format(numcurrent), "rb"))

# %%
# we use the exponential fit from investigating_lossrates.ipynb
lossratedict = pickle.load(
    open('data_processed/sensitivity_estimation/lossrate_expfit.pkl', 'rb'))
kint = lossrate(
    fixcurrent, lossratedict['ki0 (Hz)'], lossratedict['kialpha (Hz)'], lossratedict['kiix (A)'])
kext = lossrate(
    fixcurrent, lossratedict['ke0 (Hz)'], lossratedict['kealpha (Hz)'], lossratedict['keix (A)'])
ktot = kint+kext

# %%
Sin = V0
iac = meas_data['Imodpp (A)'].iloc[0]/2
Omega = meas_data['Modfrec (Hz)'].iloc[0]
Delta = meas_data.axes[0].values
Delta_model = np.linspace(min(Delta), max(Delta), 401)

# %%
# Comparing with measurement
meas_dF1, meas_signal_p1, meas_sensitivity_p1, meas_signal_m1, meas_sensitivity_m1 = calculate_measurement(
    meas_data, order=1, xkey='dF (Hz)')
meas_dF2, meas_signal_p2, meas_sensitivity_p2, meas_signal_m2, meas_sensitivity_m2 = calculate_measurement(
    meas_data, order=2, xkey='dF (Hz)')
meas_dF3, meas_signal_p3, meas_sensitivity_p3, meas_signal_m3, meas_sensitivity_m3 = calculate_measurement(
    meas_data, order=3, xkey='dF (Hz)')

# %%
"""
************************************************************************************
************************************************************************************
************************************************************************************
First order calculation
************************************************************************************
************************************************************************************
************************************************************************************
"""

# %%
# calculating the first order peak coefficients heights
cm1, cp1 = first_order_coeffs(
    iac, Omega, G1, kext, ktot, Sin, Delta_model)
Sm1 = np.sqrt(kext)*cm1
Sp1 = np.sqrt(kext)*cp1
Sm1_W = abs(Sm1)**2/(2*50)
Sm1_dBm = 10*np.log10(Sm1_W)+30+gain

idx_measmax = meas_signal_m1.values.argmax()

Plot_model_inputoutput((Delta-Delta[idx_measmax])/1e3, meas_signal_m1,
                       Delta_model/1e3, Sm1_dBm,
                       porder=1, corder=1, savefig=saveall, xvariable='dF', xlabel='Detuning (kHz)', ktot=ktot/1e3)

# %%
# save for figure
figure3_df = {'xmeas': Delta/1e3, 'ymeas': meas_signal_m1, 'ymeaslabel': 'data',
              'xtheo': Delta_model/1e3, 'ytheo': Sm1_dBm, 'ytheolabel': 'model',
              'xlabel': 'Detuning (kHz)', 'ylabel': 'Amplitude (dBm)'}
pickle.dump(figure3_df, open('data_final/fig3_panel_df.pkl', 'wb'))

# %%
"""
************************************************************************************
************************************************************************************
************************************************************************************
Second order calculation
************************************************************************************
************************************************************************************
************************************************************************************
"""

# calculating the second order peak coefficients heights
cm1, cp1, cm2, cp2 = second_order_coeffs(
    iac, Omega, G1, G2, kext, ktot, Sin, Delta_model)

Sm1 = np.sqrt(kext)*cm1
Sm1_W = abs(Sm1)**2/(2*50)
Sm1_dBm = 10*np.log10(Sm1_W)+30+gain

Sp1 = np.sqrt(kext)*cp1
Sp1_W = abs(Sp1)**2/(2*50)
Sp1_dBm = 10*np.log10(Sp1_W)+30+gain

Sm2 = np.sqrt(kext)*cm2
Sm2_W = abs(Sm2)**2/(2*50)
Sm2_dBm = 10*np.log10(Sm2_W)+30+gain
Sp2 = np.sqrt(kext)*cp2
Sp2_W = abs(Sp2)**2/(2*50)
Sp2_dBm = 10*np.log10(Sp2_W)+30+gain

# %%
Plot_model_inputoutput(Delta/1e3, meas_signal_m1,
                       Delta_model/1e3, Sm1_dBm, porder=1, corder=2, savefig=saveall, xvariable='dF', xlabel='Detuning (kHz)', ktot=ktot/1e3)

# %%
Plot_model_inputoutput(Delta/1e3, meas_signal_m2,
                       Delta_model/1e3, Sm2_dBm, porder=2, corder=2, savefig=saveall,
                       xvariable='dF', xlabel='Detuning (kHz)', ktot=ktot/1e3)

# %%
# save for figure
figure3M_df = {'xmeas': (Delta-Delta[idx_measmax])/1e3, 'ymeas': meas_signal_m2, 'ymeaslabel': 'data',
               'xtheo': Delta_model/1e3, 'ytheo': Sm2_dBm, 'ytheolabel': 'model',
               'xlabel': 'Detuning (kHz)', 'ylabel': 'Amplitude (dBm)'}
pickle.dump(figure3M_df, open('data_final/SM_fig3_panel_df.pkl', 'wb'))

# %%
"""
************************************************************************************
************************************************************************************
************************************************************************************
Third order calculation
************************************************************************************
************************************************************************************
************************************************************************************
"""

# calculating the second order peak coefficients heights
cm1, cp1, cm2, cp2, cm3, cp3 = third_order_coeffs(
    iac, Omega, G1, G2, G3, kext, ktot, Sin, Delta_model)

Sm1 = np.sqrt(kext)*cm1
Sm1_W = abs(Sm1)**2/(2*50)
Sm1_dBm = 10*np.log10(Sm1_W)+30+gain
Sp1 = np.sqrt(kext)*cp1
Sp1_W = abs(Sp1)**2/(2*50)
Sp1_dBm = 10*np.log10(Sp1_W)+30+gain

Sm2 = np.sqrt(kext)*cm2
Sm2_W = abs(Sm2)**2/(2*50)
Sm2_dBm = 10*np.log10(Sm2_W)+30+gain
Sp2 = np.sqrt(kext)*cp2
Sp2_W = abs(Sp2)**2/(2*50)
Sp2_dBm = 10*np.log10(Sp2_W)+30+gain

Sm3 = np.sqrt(kext)*cm3
Sm3_W = abs(Sm3)**2/(2*50)
Sm3_dBm = 10*np.log10(Sm3_W)+30+gain
Sp3 = np.sqrt(kext)*cp3
Sp3_W = abs(Sp3)**2/(2*50)
Sp3_dBm = 10*np.log10(Sp3_W)+30+gain

# %%
Plot_model_inputoutput(Delta/1e3, meas_signal_m1,
                       Delta_model/1e3, Sm1_dBm, porder=1, corder=3, savefig=saveall,
                       xvariable='dF', xlabel='Detuning (kHz)', ktot=ktot/1e3)

# %%
Plot_model_inputoutput(Delta/1e3, meas_signal_m2,
                       Delta_model/1e3, Sm2_dBm, porder=2, corder=3, savefig=saveall,
                       xvariable='dF', xlabel='Detuning (kHz)', ktot=ktot/1e3)

# %%
Plot_model_inputoutput(Delta/1e3, meas_signal_m3,
                       Delta_model/1e3, Sm3_dBm, porder=3, corder=3, savefig=saveall,
                       xvariable='dF', xlabel='Detuning (kHz)', ktot=ktot/1e3)
