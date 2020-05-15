import pickle
import stlabutils
from matplotlib import pyplot as plt
import numpy as np
import scipy.constants
import pandas as pd
from src.model_currentbias import f0, lossrate

# Sweep parameters

# P_start = -100  # dBm
# P_stop = -150  # dBm

# P_p_lin_start = 10**((P_start)/10)/1000  # W
# P_p_lin_stop = 10**((P_stop)/10)/1000  # W

# Probe_powers = np.linspace(P_start, P_stop, 201)
# Probe_powers_log = 10*np.log10(Probe_powers)+30
# Probe_powers_log = 10*np.log10(Probe_powers)+30

# Constants

hbar = scipy.constants.hbar
pi = np.pi

fixcurrent = 4e-6
# Cavity parameters

# we use the exponential fit from investigating_lossrates.ipynb
lossratedict = pickle.load(
    open('data_processed/sensitivity_estimation/lossrate_expfit.pkl', 'rb'))
kint = lossrate(
    fixcurrent, lossratedict['ki0 (Hz)'], lossratedict['kialpha (Hz)'], lossratedict['kiix (A)'])
kext = lossrate(
    fixcurrent, lossratedict['ke0 (Hz)'], lossratedict['kealpha (Hz)'], lossratedict['keix (A)'])
ktot = kint+kext
fr, Lr, Ic = pickle.load(open(
    'data_processed/sensitivity_estimation/popt_simulation_JJarrays.pkl', 'rb'))
f_0 = f0(fixcurrent, fr=fr, Lr=Lr, Ic=Ic)

P_out_VNA = np.linspace(-20, 30, 101)  # range from experiment

attenuation = 129  # in dB
P_in_sample = P_out_VNA - attenuation
Probe_powers = P_in_sample

k = kint + kext
beta = -2500

Probe_freqs = np.linspace(f_0 - 12e6, f_0 + 3e6, 1001)

colnames = [
    'Probe Frequency (Hz)', 'S11re ()', 'S11im ()', 'S11dB (dB)',
    'S11ph (rad)', 'S11abs ()', 'Intracavity photons ()', 'Probe power (W)', 'Probe power (dBm)',
    'Beta (Hz)'
]
for sol in ['min', 'max', 'med']:
    # sol = 'min'
    myfile = stlabutils.newfile("F",
                                "ResponsePowerDep_SingleTone_" + sol,
                                mypath='Duffing/',
                                usedate=False,
                                usefolder=False,
                                colnames=colnames,
                                git_id=False)

    for i, P_p_log in enumerate(Probe_powers):

        for ii, f_p in enumerate(Probe_freqs):
            P_p_lin = 10**(P_p_log / 10) * 1e-3  # dBm to W
            # On-chip Probe Power

            # On-chip photon flux
            n_p = P_p_lin / hbar / 2 / pi / f_p

            # Intracavity photon number
            a = 4 * pi**2 * beta**2
            b = -2 * 2 * pi * (f_p - f_0) * 2 * pi * beta
            c = (2 * pi * (f_p - f_0))**2 + (2 * pi * k)**2 / 4
            d = -2 * pi * kext * n_p
            coeff = [a, b, c, d]
            a_1, a_2, a_3 = np.roots(coeff)

            # Replace complex roots by -1
            roots = [a_1, a_2, a_3]
            roots = [i.real if i.imag == 0 else -1 for i in roots]

            # Choose branch:
            # sol=='min': find smallest real and positive solution for the low-amplitude branch
            # sol=='max': find largest real and positive solution for the high-amplitude branch
            # sol=='med': find intermediate branch value
            if sol == 'min':
                alpha_0 = min([alpha for alpha in roots if alpha > 0])
            elif sol == 'max':
                alpha_0 = max([alpha for alpha in roots if alpha > 0])
            elif sol == 'med':
                alpha_0 = np.median([alpha for alpha in roots if alpha > 0])
            # print(roots, alpha_0)

            # Calculate response
            Delta = f_p - f_0
            phi = np.arctan(-2 * (Delta - beta * alpha_0) / k)
            S11 = 1 - np.sqrt(2 * pi * kext * alpha_0) / np.sqrt(n_p) * np.exp(
                -1j * phi)

            line = [
                f_p,
                np.real(S11),
                np.imag(S11),
                20 * np.log10(np.abs(S11)),
                np.angle(S11),
                np.abs(S11),
                alpha_0,
                P_p_lin,
                P_p_log,
                beta
            ]
            stlabutils.writeline(myfile, line)

            # myfile.write('\n')
        myfile.write('\n')
        stlabutils.utils.metagen.fromarrays(myfile,
                                            Probe_freqs,
                                            Probe_powers[0:i + 1],
                                            xtitle='Probe Frequency (Hz)',
                                            ytitle='Probe power (dBm)',
                                            colnames=colnames)

    myfile.close()
