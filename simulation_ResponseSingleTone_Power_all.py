"""
We here simulate the Duffing response of a driven resonator, where we choose the values for beta based on previous comparison of data and measurement.
"""

import pickle
import stlabutils
from matplotlib import pyplot as plt
import numpy as np
import scipy.constants
import pandas as pd
import time
from src.model_currentbias import Duffing_parameter
from src.model_currentbias import f0 as f0model
from src.model_currentbias import lossrate
from scipy.interpolate import interp1d
import os

# Constants

hbar = scipy.constants.hbar
pi = np.pi

# Cavity parameters
currs = np.arange(0, 7.7e-6, 0.1e-6)
currlist = [int(round(x * 1e9)) for x in currs]
fr, Lr, Ic = pickle.load(open(
    'data_processed/sensitivity_estimation/popt_simulation_JJarrays.pkl', 'rb'))
f0theo = f0model(currs, fr=fr, Lr=Lr, Ic=Ic)
G1theo = np.gradient(f0theo, currs)
# we use the exponential fit from investigating_lossrates.ipynb
lossratedict = pickle.load(
    open('data_processed/sensitivity_estimation/lossrate_expfit.pkl', 'rb'))
kinttheo = lossrate(
    currs, lossratedict['ki0 (Hz)'], lossratedict['kialpha (Hz)'], lossratedict['kiix (A)'])
kexttheo = lossrate(
    currs, lossratedict['ke0 (Hz)'], lossratedict['kealpha (Hz)'], lossratedict['keix (A)'])
ktottheo = kinttheo+kexttheo

P_out_VNA = np.linspace(-20, 30, 101)  # range from experiment
power_in = pickle.load(
    open("data_processed/sensitivity_estimation/power_in.pkl", "rb"))
attenuation = -power_in['P_in_sample (dBm)']
P_in_sample = P_out_VNA - attenuation
Probe_powers = P_in_sample
Probe_powers_lin = 10**(Probe_powers / 10) * 1e-3  # dBm to W

# sol = 'min'

betaexport = pickle.load(open('data_processed/Duffing/beta_export.pkl', 'rb'))
# print(min(currs), max(currs), min(betaexport['curr']), max(betaexport['curr']))
betaint = interp1d(betaexport['curr'], betaexport['betaexp'])(currs)

for kk, numcurrent in enumerate(currlist):

    start_time = time.time()

    mypath = 'Duffing/highres_1JJ/'
    prefix = "F"
    idstring = "ResponsePowerDep_SingleTone_1JJ_{:04d}".format(numcurrent)
    filepath = mypath+prefix+'_'+idstring+".dat"
    print(filepath)
    if True:  # not os.path.isfile(filepath):

        Probe_freqs = np.linspace(f0theo[kk] - 20e6, f0theo[kk] + 5e6, 1001)
        freqlength = len(Probe_freqs)

        for i, (P_p_log, P_p_lin) in enumerate(zip(Probe_powers, Probe_powers_lin)):

            # On-chip photon flux
            n_p = P_p_lin / hbar / 2 / pi / Probe_freqs

            # Intracavity photon number
            a = 4 * pi**2 * betaint[kk]**2
            b = -2 * 2 * pi * (Probe_freqs - f0theo[kk]) * 2 * pi * betaint[kk]
            c = (2 * pi * (Probe_freqs - f0theo[kk])
                 )**2 + (2 * pi * ktottheo[kk])**2 / 4
            d = -2 * pi * kexttheo[kk] * n_p
            coeffs = [[a, x2, x3, x4]
                      for (x2, x3, x4) in zip(b, c, d)]

            allroots = [np.roots(coeff) for coeff in coeffs]
            allroots = [[i.real if i.imag == 0 else np.nan for i in roots]
                        for roots in allroots]
            alpha_0 = np.array(
                [np.nanmin(alpha) for alpha in allroots])

            # Calculate response
            Delta = Probe_freqs - f0theo[kk]
            phi = np.arctan(-2 *
                            (Delta - betaint[kk] * alpha_0) / ktottheo[kk])
            S11 = 1 - np.sqrt(2 * pi * kexttheo[kk] * alpha_0) / np.sqrt(
                n_p) * np.exp(-1j * phi)

            mydf = pd.DataFrame({
                'Probe Frequency (Hz)': Probe_freqs,
                'Probe power (W)': [P_p_lin]*freqlength,
                'Probe power (dBm)': [P_p_log]*freqlength,
                'S11dB (dB)': 20*np.log10(abs(S11)),
                'S11ph (rad)': np.angle(S11),
                'Intracavity photons ()': alpha_0,
                'Beta (Hz)': [betaint[kk]]*freqlength,
                'Delta (Hz)': Delta,
                'Phi (rad)': phi,
                'Iset (A)': currs[kk],
                'kint (Hz)': kinttheo[kk],
                'kext (Hz)': kexttheo[kk],
                'f0 (Hz)': f0theo[kk],
                'G1 (Hz/A)': G1theo[kk]
            })

            if i == 0:
                myfile = stlabutils.newfile(
                    prefix,
                    idstring,
                    mypath=mypath,
                    usedate=False,
                    usefolder=False,
                    colnames=list(mydf),
                    git_id=False)

            stlabutils.saveframe(myfile, mydf)
            stlabutils.utils.metagen.fromarrays(
                myfile,
                Probe_freqs,
                Probe_powers[0:i + 1],
                xtitle='Probe Frequency (Hz)',
                ytitle='Probe power (dBm)',
                colnames=list(mydf))

        myfile.close()
        elapsed_time = time.time() - start_time
        print('Elapsed time: {:.2f}s'.format(elapsed_time))
