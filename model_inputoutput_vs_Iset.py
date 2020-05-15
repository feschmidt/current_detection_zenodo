# -*- coding: utf-8 -*-
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

"""
This file calculates and plots the measured peak heights together with the ones calculated from input-output theory, as a function of varying bias current.
"""

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

"""
************************************************************************************
************************************************************************************
************************************************************************************
Calculating the peak height using the input-output formalism
************************************************************************************
************************************************************************************
************************************************************************************
"""

saveall = True

# Loading in data and parameters
power_gain = pickle.load(
    open("data_processed/sensitivity_estimation/gaindata.pkl", "rb"))
gain = power_gain['Gain (dB)']

power_in = pickle.load(
    open("data_processed/sensitivity_estimation/power_in.pkl", "rb"))

# Load fitparameters for model
popt = pickle.load(open(
    'data_processed/sensitivity_estimation/popt_simulation_JJarrays.pkl', 'rb'))

# Comparing with measurement
meas_data = pickle.load(
    open("data_processed/sensitivity_estimation/data_F33B_peaks.pkl", "rb"))
meas_curr1, meas_signal_p1, meas_sensitivity_p1, meas_signal_m1, meas_sensitivity_m1 = calculate_measurement(
    meas_data, order=1)
meas_curr2, meas_signal_p2, meas_sensitivity_p2, meas_signal_m2, meas_sensitivity_m2 = calculate_measurement(
    meas_data, order=2)
meas_curr3, meas_signal_p3, meas_sensitivity_p3, meas_signal_m3, meas_sensitivity_m3 = calculate_measurement(
    meas_data, order=3)

# Defining the individual terms to match the variable names in the formulas, and create the model
Iset_meas = meas_data['Is (A)']
Iset = np.linspace(min(Iset_meas), max(Iset_meas), 401)
f0_of_I = f0(Iset, *popt, nJJ=1)
G1_of_I = np.gradient(f0_of_I, Iset)
G2_of_I = np.gradient(G1_of_I, Iset)
G3_of_I = np.gradient(G2_of_I, Iset)
Sin = power_in['V_in_sample (V)']
iac = meas_data['Imodpp (A)'].iloc[0]/2
Omega = meas_data['Modfrec (Hz)'].iloc[0]

# instead of interpolating the loss rates, we use the exponential fit from investigating_lossrates.ipynb
lossratedict = pickle.load(
    open('data_processed/sensitivity_estimation/lossrate_expfit.pkl', 'rb'))
# there's no need for interpolation, but this way we don't have to change the rest of the code
ki_fit = interp1d(Iset, lossrate(
    Iset, lossratedict['ki0 (Hz)'], lossratedict['kialpha (Hz)'], lossratedict['kiix (A)']))
ke_fit = interp1d(Iset, lossrate(
    Iset, lossratedict['ke0 (Hz)'], lossratedict['kealpha (Hz)'], lossratedict['keix (A)']))
ktot_fit = ki_fit(Iset)+ke_fit(Iset)

"""
************************************************************************************
************************************************************************************
************************************************************************************
First order calculation
************************************************************************************
************************************************************************************
************************************************************************************
"""

# calculating the first order peak coefficients heights
cm1, cp1 = first_order_coeffs(
    iac, Omega, G1_of_I, ke_fit(0), ktot_fit[0], Sin, Delta=0)
cm1_fit, cp1_fit = first_order_coeffs(
    iac, Omega, G1_of_I, ke_fit(Iset), ktot_fit, Sin, Delta=0)
Sm1 = np.sqrt(ke_fit(0))*cm1
Sm1_W = abs(Sm1)**2/(2*50)
Sm1_dBm = 10*np.log10(Sm1_W)+30+gain
Sm1_fit = np.sqrt(ke_fit(Iset))*cm1_fit
Sm1_W_fit = abs(Sm1_fit)**2/(2*50)
Sm1_dBm_fit = 10*np.log10(Sm1_W_fit)+30+gain

# Final plot
Plot_model_inputoutput(meas_curr1/1e-6, meas_signal_m1,
                       Iset/1e-6, Sm1_dBm, Is_model2=Iset/1e-6, Signal_model2=Sm1_dBm_fit, porder=1, corder=1, savefig=saveall)
# save for figure
figure3_Is = {'xmeas': meas_curr1/1e-6, 'ymeas': meas_signal_m1, 'ymeaslabel': 'data',
              'xtheo': Iset/1e-6, 'ytheo': Sm1_dBm, 'ytheolabel': 'model ki fix',
              'xtheo2': Iset/1e-6, 'ytheo2': Sm1_dBm_fit, 'ytheolabel2': 'model ki varies',
              'xlabel': r'Bias current (µA)', 'ylabel': 'Amplitude (dBm)'}
pickle.dump(figure3_Is, open('data_final/fig3_panel_Is.pkl', 'wb'))

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
_, _, cm2, cp2 = second_order_coeffs(
    iac, Omega, G1_of_I, G2_of_I, ke_fit(0), ktot_fit[0], Sin, Delta=0)
_, _, cm2_fit, cp2_fit = second_order_coeffs(
    iac, Omega, G1_of_I, G2_of_I, ke_fit(Iset), ktot_fit, Sin, Delta=0)

Sm2 = np.sqrt(ke_fit(0))*cm2
Sm2_W = abs(Sm2)**2/(2*50)
Sm2_dBm = 10*np.log10(Sm2_W)+30+gain

Sm2_fit = np.sqrt(ke_fit(Iset))*cm2_fit
Sm2_W_fit = abs(Sm2_fit)**2/(2*50)
Sm2_dBm_fit = 10*np.log10(Sm2_W_fit)+30+gain

# +
# Final plot
Plot_model_inputoutput(meas_curr2/1e-6, meas_signal_m2,
                       Iset/1e-6, Sm2_dBm, Is_model2=Iset/1e-6, Signal_model2=Sm2_dBm_fit, porder=2, corder=2, savefig=saveall)
# Final_plot_inputoutput(Iset/1e-6, Sm1, Sp1, Sm1_dBm, meas_curr1/1e-6,
#                        meas_signal_p1, corder=2, porder=1)
# Final_plot_inputoutput(Iset/1e-6, Sm2, Sp2, Sp2_dBm, meas_curr2/1e-6,
#                        meas_signal_p2, corder=2, porder=2)

# save for figure
figure3M_Is = {'xmeas': meas_curr2/1e-6, 'ymeas': meas_signal_m2, 'ymeaslabel': 'data',
               'xtheo': Iset/1e-6, 'ytheo': Sm2_dBm, 'ytheolabel': 'model ki fix',
               'xtheo2': Iset/1e-6, 'ytheo2': Sm2_dBm_fit, 'ytheolabel2': 'model ki varies',
               'xlabel': r'Bias current (µA)', 'ylabel': 'Amplitude (dBm)'}
pickle.dump(figure3M_Is, open('data_final/SM_fig3_panel_Is.pkl', 'wb'))
# -

"""
************************************************************************************
************************************************************************************
************************************************************************************
Third order calculation
************************************************************************************
************************************************************************************
************************************************************************************
"""

# calculating the third order peak coefficients heights
_, _, _, _, cm3, cp3 = third_order_coeffs(
    iac, Omega, G1_of_I, G2_of_I, G3_of_I, ke_fit(0), ktot_fit[0], Sin, Delta=0)
_, _, _, _, cm3_fit, cp3_fit = third_order_coeffs(
    iac, Omega, G1_of_I, G2_of_I, G3_of_I, ke_fit(Iset), ktot_fit, Sin, Delta=0)

Sm3 = np.sqrt(ke_fit(0))*cm3
Sm3_W = abs(Sm3)**2/(2*50)
Sm3_dBm = 10*np.log10(Sm3_W)+30+gain
Sm3_fit = np.sqrt(ke_fit(Iset))*cm3_fit
Sm3_W_fit = abs(Sm3_fit)**2/(2*50)
Sm3_dBm_fit = 10*np.log10(Sm3_W_fit)+30+gain

# Final plot
Plot_model_inputoutput(meas_curr3/1e-6, meas_signal_m3,
                       Iset/1e-6, Sm3_dBm, Is_model2=Iset/1e-6, Signal_model2=Sm3_dBm_fit, porder=3, corder=3, savefig=saveall)

# Final_plot_inputoutput(Iset/1e-6, Sm1, Sp1, Sm1_dBm, meas_curr1/1e-6,
#                        meas_signal_p1, corder=3, porder=1)
# Final_plot_inputoutput(Iset/1e-6, Sm2, Sp2, Sp2_dBm, meas_curr2/1e-6,
#                        meas_signal_p2, corder=3, porder=2)
# Final_plot_inputoutput(Iset/1e-6, Sm3, Sp3, Sp3_dBm, meas_curr3/1e-6,
#                        meas_signal_p3, corder=3, porder=3)
