# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

"""
This file calculates and plots the measured peak heights together with the ones calculated from input-output theory, as a function of fixed bias current and verying pump power.
It needs ```processing_dnfdIn_VNA_vs_power_vs_current.py``` to be run before and is the last file to be executed.
It finally generates the overview plots
This file is identical to `calculation_inputoutput_vs_Ppump.py`, except that here we correct the intracavity photon number by the value simulated using the Duffing model.
"""

import glob
import stlabutils
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
import pickle
from src.algo_inputoutput import first_order_coeffs, second_order_coeffs, third_order_coeffs, nphot_to_a0
from src.algo_peakheights import calculate_measurement, func_PdBmtoV, func_cmpx_to_dBm
from src.plot_dnSdfn import Plot_model_inputoutput
from src.model_currentbias import f0
from src.model_currentbias import G1 as G1func
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

Iac = 10e-9
Omega = 1e3

# Loading in data and parameters
power_gain = pickle.load(
    open("data_processed/sensitivity_estimation/gaindata.pkl", "rb"))
gain = power_gain['Gain (dB)']

simfiles = sorted(glob.glob('Duffing/highres_1JJ/*_1JJ_*.dat'))

fpumpfile = sorted(glob.glob('data_plots/pos2/F33_pwr_f0+0f1_mtx_frequency.pkl'))
fpumpdat = pickle.load(open(fpumpfile[0],'rb'))
fpumpf0s = fpumpdat.pmtx

fpumpf0s.values[:,0]

plt.imshow(fpumpf0s.values)

plt.plot(fpumpf0s.axes[1].values,fpumpf0s.values[10,:],'.-')
plt.axvline(7.3)
plt.axvline(7.5)
plt.xlim(6,8)

myf0func = interp1d(fpumpf0s.axes[1].values*1e-6,fpumpf0s.values[10,:])

extrashift = np.arange(-600e3,501e3,100e3)
extrashift

for theshift in extrashift:
    shiftid = '%+04dkHz' % (int(theshift/1e3))
    print(shiftid)
    
    for simfile_min in simfiles:
        print(simfile_min)
        simdata_min = stlabutils.readdata.readdat(simfile_min)

        fixcurrent = simdata_min[0]['Iset (A)'][0]
        numcurrent = round(int(fixcurrent/1e-9), -1)
        fpump = simdata_min[0]['f0 (Hz)'][0]+theshift # fitted resonance
    #     fpump = myf0func(fixcurrent) # measured resonance
        G1 = simdata_min[0]['G1 (Hz/A)'][0]
        kext = simdata_min[0]['kext (Hz)'][0]
        kint = simdata_min[0]['kint (Hz)'][0]
        ktot = kext+kint

        # From the simulation, extracting actual photon number at pump frequency
        nph_int = []
        nph_int_max = []
        pwr = []
        fsim_min = []

        for simblock in simdata_min:
            thepower = simblock['Probe power (dBm)'].iloc[0]
            nph = simblock['Intracavity photons ()'].values
            freq = simblock['Probe Frequency (Hz)'].values
            yy = simblock['S11dB (dB)'].values
            idx = yy.argmin()
            fsim_min.append(freq[idx])
            try:
                nph_int.append(interp1d(freq, nph)(fpump))
                nph_int_max.append(interp1d(freq, nph)(freq[idx]))
            except ValueError:
                print('Warning! fpump is outside of the simulated frequency range!')
                n1 = InterpolatedUnivariateSpline(freq, nph)(fpump)
                if n1 > 0:
                    nph_int.append(n1)
                else:
                    nph_int.append(0)
                n2 = InterpolatedUnivariateSpline(
                    freq, nph)(freq[idx])
                if n2 > 0:
                    nph_int_max.append(n2)
                else:
                    nph_int_max.append(0)
            pwr.append(thepower)
        nph_int = np.array(nph_int)
        nph_int_max = np.array(nph_int_max)
        pwr = np.array(pwr)
        # converting intracavity photon number to voltage amplitude
        alpha_0 = nphot_to_a0(nph_int, fpump)
        fsim_min = np.array(fsim_min)
        alpha_0_max = nphot_to_a0(nph_int_max, fsim_min)

        # Defining the individual terms to match the variable names in the formulas
        iac = 10e-9
        Omega = 1e3

        # corrections compared to non-Duffing analysis
        # Using the same cavity parameters as for the simulation
        beta = simdata_min[0]['Beta (Hz)'][0]
        fres_Duffing = fpump+2*beta*nph_int
        Delta_Duffing = fpump-fres_Duffing+theshift

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
            iac, Omega, G1, kext, ktot, Sin=0, Delta=Delta_Duffing, a0_calc=alpha_0)
        Sm1_dBm = func_cmpx_to_dBm(cm1, kext, gain)
        Sp1_dBm = func_cmpx_to_dBm(cp1, kext, gain)

        # same but without detuning
        cm1_nodet, cp1_nodet = first_order_coeffs(
            iac, Omega, G1, kext, ktot, Sin=0, Delta=0, a0_calc=alpha_0_max)
        Sm1_dBm_nodet = func_cmpx_to_dBm(cm1_nodet, kext, gain)
        Sp1_dBm_nodet = func_cmpx_to_dBm(cp1_nodet, kext, gain)

        # save for figure
        simpkl = {'xtheo': pwr, 'ytheo': Sm1_dBm, 'ytheolabel': 'model',
                  'xtheo2': pwr, 'ytheo2': Sm1_dBm_nodet, 'ytheolabel2': 'model on resonance',
                  'xtheop': pwr, 'ytheop': Sp1_dBm, 'ytheolabelp': 'modelp',
                  'xtheo2p': pwr, 'ytheo2p': Sp1_dBm_nodet, 'ytheolabel2p': 'modelp on resonance',
                  'xlabel': 'Power at sample (dBm)', 'ylabel': 'Amplitude (dBm)',
                  'Iset (A)': fixcurrent, 'Beta (Hz)': beta,
                  # carry over old values
                  'kext (Hz)': kext, 'kint (Hz)': kint, 'fpump (Hz)': fpump, 'G1 (Hz/A)': G1}
        exportpath = 'data_final/fig3_panel_Ppump_1JJv2detuningcheck_'+shiftid+'{:04d}.pkl'.format(
            numcurrent)
        pickle.dump(simpkl, open(exportpath, 'wb'))
        print('Exported keys:', simpkl.keys())
        print('Exported file:', exportpath)
        print()

for simfile_min in simfiles:
    print(simfile_min)
    simdata_min = stlabutils.readdata.readdat(simfile_min)

    fixcurrent = simdata_min[0]['Iset (A)'][0]
    numcurrent = round(int(fixcurrent/1e-9), -1)
    fpump = simdata_min[0]['f0 (Hz)'][0]+theshift # fitted resonance
#     fpump = myf0func(fixcurrent) # measured resonance
    G1 = simdata_min[0]['G1 (Hz/A)'][0]
    kext = simdata_min[0]['kext (Hz)'][0]
    kint = simdata_min[0]['kint (Hz)'][0]
    ktot = kext+kint

    # From the simulation, extracting actual photon number at pump frequency
    nph_int = []
    nph_int_max = []
    pwr = []
    fsim_min = []

    for simblock in simdata_min:
        thepower = simblock['Probe power (dBm)'].iloc[0]
        nph = simblock['Intracavity photons ()'].values
        freq = simblock['Probe Frequency (Hz)'].values
        yy = simblock['S11dB (dB)'].values
        idx = yy.argmin()
        fsim_min.append(freq[idx])
        try:
            nph_int.append(interp1d(freq, nph)(fpump))
            nph_int_max.append(interp1d(freq, nph)(freq[idx]))
        except ValueError:
            print('Warning! fpump is outside of the simulated frequency range!')
            n1 = InterpolatedUnivariateSpline(freq, nph)(fpump)
            if n1 > 0:
                nph_int.append(n1)
            else:
                nph_int.append(0)
            n2 = InterpolatedUnivariateSpline(
                freq, nph)(freq[idx])
            if n2 > 0:
                nph_int_max.append(n2)
            else:
                nph_int_max.append(0)
        pwr.append(thepower)
    nph_int = np.array(nph_int)
    nph_int_max = np.array(nph_int_max)
    pwr = np.array(pwr)
    # converting intracavity photon number to voltage amplitude
    alpha_0 = nphot_to_a0(nph_int, fpump)
    fsim_min = np.array(fsim_min)
    alpha_0_max = nphot_to_a0(nph_int_max, fsim_min)

    # Defining the individual terms to match the variable names in the formulas
    iac = 10e-9
    Omega = 1e3

    # corrections compared to non-Duffing analysis
    # Using the same cavity parameters as for the simulation
    beta = simdata_min[0]['Beta (Hz)'][0]
    fres_Duffing = fpump+2*beta*nph_int
    Delta_Duffing = fpump-fres_Duffing+theshift

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
        iac, Omega, G1, kext, ktot, Sin=0, Delta=Delta_Duffing, a0_calc=alpha_0)
    Sm1_dBm = func_cmpx_to_dBm(cm1, kext, gain)
    Sp1_dBm = func_cmpx_to_dBm(cp1, kext, gain)

    # same but without detuning
    cm1_nodet, cp1_nodet = first_order_coeffs(
        iac, Omega, G1, kext, ktot, Sin=0, Delta=0, a0_calc=alpha_0_max)
    Sm1_dBm_nodet = func_cmpx_to_dBm(cm1_nodet, kext, gain)
    Sp1_dBm_nodet = func_cmpx_to_dBm(cp1_nodet, kext, gain)

    # save for figure
    simpkl = {'xtheo': pwr, 'ytheo': Sm1_dBm, 'ytheolabel': 'model',
              'xtheo2': pwr, 'ytheo2': Sm1_dBm_nodet, 'ytheolabel2': 'model on resonance',
              'xtheop': pwr, 'ytheop': Sp1_dBm, 'ytheolabelp': 'modelp',
              'xtheo2p': pwr, 'ytheo2p': Sp1_dBm_nodet, 'ytheolabel2p': 'modelp on resonance',
              'xlabel': 'Power at sample (dBm)', 'ylabel': 'Amplitude (dBm)',
              'Iset (A)': fixcurrent, 'Beta (Hz)': beta,
              # carry over old values
              'kext (Hz)': kext, 'kint (Hz)': kint, 'fpump (Hz)': fpump, 'G1 (Hz/A)': G1}
    exportpath = 'data_final/fig3_panel_Ppump_1JJv2detuningcheck_'+shiftid+'{:04d}.pkl'.format(
        numcurrent)
    pickle.dump(simpkl, open(exportpath, 'wb'))
    print('Exported keys:', simpkl.keys())
    print('Exported file:', exportpath)
    print()


