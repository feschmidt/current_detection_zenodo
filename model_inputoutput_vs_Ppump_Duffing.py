"""
This file calculates and plots the measured peak heights together with the ones calculated from input-output theory, as a function of fixed bias current and verying pump power.
It finally generates the overview plots
"""

import glob
import stlabutils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
import pickle
from src.algo_inputoutput import first_order_coeffs, second_order_coeffs, third_order_coeffs, nphot_to_a0
from src.algo_peakheights import calculate_measurement, func_PdBmtoV, func_cmpx_to_dBm
from src.plot_dnSdfn import Plot_model_inputoutput
from src.model_currentbias import f0, lossrate
from src.algo_inputoutput import c_a0

"""
************************************************************************************
************************************************************************************
************************************************************************************
Calculating the peak height using the input-output formalism
at fixed current for varying pump power
************************************************************************************
************************************************************************************
************************************************************************************
"""

saveall = True
showall = True

fixcurrent = 4e-6
numcurrent = round(fixcurrent/1e-9)

# Loading in data and parameters
power_gain = pickle.load(
    open("data_processed/sensitivity_estimation/gaindata.pkl", "rb"))
gain = power_gain['Gain (dB)']

power_in = pickle.load(
    open("data_processed/sensitivity_estimation/power_in.pkl", "rb"))

fitpars = pickle.load(
    open("data_processed/sensitivity_estimation/fitpars_vs_Is.pkl", "rb"))
f0_of_Ifix = fitpars.loc[fixcurrent]

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
f0fit0 = f0(fixcurrent, *popt, nJJ=1)
G1 = interp1d(Iset, G1_of_I)(fixcurrent)
G2 = interp1d(Iset, G2_of_I)(fixcurrent)
G3 = interp1d(Iset, G3_of_I)(fixcurrent)

# we use the exponential fit from investigating_lossrates.ipynb
lossratedict = pickle.load(
    open('data_processed/sensitivity_estimation/lossrate_expfit.pkl', 'rb'))
kint = lossrate(
    fixcurrent, lossratedict['ki0 (Hz)'], lossratedict['kialpha (Hz)'], lossratedict['kiix (A)'])
kext = lossrate(
    fixcurrent, lossratedict['ke0 (Hz)'], lossratedict['kealpha (Hz)'], lossratedict['keix (A)'])
ktot = kint+kext

# measurement data
meas_data = pickle.load(open(
    "data_processed/sensitivity_estimation/data_F33E_{}_peaks.pkl".format(numcurrent), "rb"))

# recalculating the input voltage
Pout_SG = meas_data.axes[0].values
Delta_PSG = Pout_SG+power_in['diff_attenuation (dB)']
P_in_sample = power_in['P_in_sample (dBm)']
V0 = func_PdBmtoV(P_in_sample+Delta_PSG)

# extracting the corrected intracavity photon number
simfile_min = glob.glob('Duffing/F_ResponsePowerDep_SingleTone_min.dat')[0]
simdata_min = stlabutils.readdata.readdat(simfile_min)

f_pump = interp1d(P_in_sample+Delta_PSG, [f0fit0]*len(V0))
# From the simulation, extracting actual photon number at pump frequency
nph_int = []
pwr = []
fsim_min = []

for simblock in simdata_min:
    #     simblock = simblock.set_index('Probe Frequency (Hz)')
    thepower = simblock['Probe power (dBm)'].iloc[0]
    nph = simblock['Intracavity photons ()'].values
    freq = simblock['Probe Frequency (Hz)'].values
    yy = simblock['S11dB (dB)'].values
    idx = yy.argmin()
    fsim_min.append(freq[idx])
    nph_int.append(interp1d(freq, nph)(f_pump(thepower)))
    pwr.append(thepower)
nph_int = np.array(nph_int)
pwr = np.array(pwr)
# converting intracavity photon number to voltage amplitude
alpha_0 = nphot_to_a0(nph_int, f_pump(pwr))

# Defining the individual terms to match the variable names in the formulas
iac = meas_data['Imodpp (A)'].iloc[0]/2
Omega = meas_data['Modfrec (Hz)'].iloc[0]
V0 = interp1d(P_in_sample+Delta_PSG, func_PdBmtoV(P_in_sample+Delta_PSG))
Sin = V0(pwr)

# Comparing with measurement
meas_pow1, meas_signal_p1, meas_sensitivity_p1, meas_signal_m1, meas_sensitivity_m1 = calculate_measurement(
    meas_data, order=1, xkey='Carrier Power (dBm)')
meas_pow2, meas_signal_p2, meas_sensitivity_p2, meas_signal_m2, meas_sensitivity_m2 = calculate_measurement(
    meas_data, order=2, xkey='Carrier Power (dBm)')
meas_pow3, meas_signal_p3, meas_sensitivity_p3, meas_signal_m3, meas_sensitivity_m3 = calculate_measurement(
    meas_data, order=3, xkey='Carrier Power (dBm)')

# corrections compared to non-Duffing analysis
# Using the same cavity parameters as for the simulation
beta = simdata_min[0]['Beta (Hz)'][0]
fres_Duffing = f0fit0+2*beta*nph_int
Delta_Duffing = f_pump(pwr)-fres_Duffing


func_PdBmtoV(P_in_sample+Delta_PSG).shape

"""
************************************************************************************
************************************************************************************
************************************************************************************
Comparing intracavity fields with and without Duffing
Comparing frequency detuning
************************************************************************************
************************************************************************************
************************************************************************************
"""

Delta_Duffing.shape

# intracavity fields
P2 = interp1d(P_in_sample+Delta_PSG, P_in_sample+Delta_PSG)
alpha0_0 = c_a0(kext, ktot, Sin=func_PdBmtoV(P2(pwr)),
                Delta=Delta_Duffing, a0_calc=False)  # will be complex
plt.plot(pwr, abs(alpha0_0), label='input-output')
plt.plot(pwr, alpha_0, label='input-output with Duffing')
plt.xlabel('Input power (dBm)')
plt.ylabel('Amplitude (V)')
plt.legend()
plt.title('Intracavity field with and without Duffing correction')
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
# plt.savefig('Duffing/comparison_intracavity_field.png', bbox_to_inches='tight')
plt.show()
plt.close()

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
cm1, _ = first_order_coeffs(
    iac, Omega, G1, kext, ktot, Sin, Delta=Delta_Duffing, a0_calc=alpha_0)
Sm1_dBm = func_cmpx_to_dBm(cm1, kext, gain)

# same but without detuning
cm1_nodet, _ = first_order_coeffs(
    iac, Omega, G1, kext, ktot, Sin, Delta=0, a0_calc=alpha_0)
Sm1_dBm_nodet = func_cmpx_to_dBm(cm1_nodet, kext, gain)

Plot_model_inputoutput(P_in_sample+Delta_PSG, meas_signal_m1, pwr, Sm1_dBm, Is_model2=pwr, Signal_model2=Sm1_dBm_nodet,
                       porder=1, corder=1, savefig=saveall, showfig=showall,
                       xvariable='Ppump_Duffing', xlabel='Power at sample (dBm)',
                       label1='model', label2='model on resonance')

# save for figure
figure3_Ppump = {'xmeas': P_in_sample+Delta_PSG, 'ymeas': meas_signal_m1, 'ymeaslabel': 'data',
                 'xtheo': pwr, 'ytheo': Sm1_dBm, 'ytheolabel': 'model',
                 'xtheo2': pwr, 'ytheo2': Sm1_dBm_nodet, 'ytheolabel2': 'model on resonance',
                 'xlabel': 'Power at sample (dBm)', 'ylabel': 'Amplitude (dBm)'}
pickle.dump(figure3_Ppump, open('data_final/fig3_panel_Ppump.pkl', 'wb'))

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
    iac, Omega, G1, G2, kext, ktot, Sin, Delta=Delta_Duffing, a0_calc=alpha_0)
# same but without detuning
cm1_nodet, cp1_nodet, cm2_nodet, cp2_nodet_ = second_order_coeffs(
    iac, Omega, G1, G2, kext, ktot, Sin, Delta=0, a0_calc=alpha_0)

Sm1_dBm = func_cmpx_to_dBm(cm1, kext, gain)
Sm2_dBm = func_cmpx_to_dBm(cm2, kext, gain)
Sm2_dBm_nodet = func_cmpx_to_dBm(cm2_nodet, kext, gain)

Plot_model_inputoutput(P_in_sample+Delta_PSG, meas_signal_m1, pwr, Sm1_dBm, Is_model2=pwr, Signal_model2=Sm1_dBm_nodet,
                       porder=1, corder=2, savefig=saveall, showfig=showall,
                       xvariable='Ppump_Duffing', xlabel='Power at sample (dBm)')
Plot_model_inputoutput(P_in_sample+Delta_PSG, meas_signal_m2, pwr, Sm2_dBm, Is_model2=pwr, Signal_model2=Sm2_dBm_nodet,
                       porder=2, corder=2, savefig=saveall, showfig=showall,
                       xvariable='Ppump_Duffing', xlabel='Power at sample (dBm)')

# save for figure
figure3M_Ppump = {'xmeas': P_in_sample+Delta_PSG, 'ymeas': meas_signal_m2, 'ymeaslabel': 'data',
                  'xtheo': pwr, 'ytheo': Sm2_dBm, 'ytheolabel': 'model',
                  'xtheo2': pwr, 'ytheo2': Sm2_dBm_nodet, 'ytheolabel2': 'model on resonance',
                  'xlabel': 'Power at sample (dBm)', 'ylabel': 'Amplitude (dBm)'}
pickle.dump(figure3M_Ppump, open('data_final/SM_fig3_panel_Ppump.pkl', 'wb'))

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
    iac, Omega, G1, G2, G3, kext, ktot, Sin, Delta=Delta_Duffing, a0_calc=alpha_0)

Sm1_dBm = func_cmpx_to_dBm(cm1, kext, gain)
Sm2_dBm = func_cmpx_to_dBm(cm2, kext, gain)
Sm3_dBm = func_cmpx_to_dBm(cm3, kext, gain)

Plot_model_inputoutput(P_in_sample+Delta_PSG, meas_signal_m1, pwr, Sm1_dBm,
                       porder=1, corder=3, savefig=saveall, showfig=showall,
                       xvariable='Ppump_Duffing', xlabel='Power at sample (dBm)')
Plot_model_inputoutput(P_in_sample+Delta_PSG, meas_signal_m2, pwr, Sm2_dBm,
                       porder=2, corder=3, savefig=saveall, showfig=showall,
                       xvariable='Ppump_Duffing', xlabel='Power at sample (dBm)')
Plot_model_inputoutput(P_in_sample+Delta_PSG, meas_signal_m3, pwr, Sm3_dBm,
                       porder=3, corder=3, savefig=saveall, showfig=showall,
                       xvariable='Ppump_Duffing', xlabel='Power at sample (dBm)')
