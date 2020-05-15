# Algorithms for calculating the expected sensitivity for the measured cavity response
# d<n><X>d<Y><n> denotes the <n>th partial derivative of quantity <X> with respect to quantity <Y>, i.e. $\partial^{n}_{Y} X$
# Extracting resonance fequency and S11 line by line from raw dataset (megameasurement-A)

try:
    from .algo_sensitivity import Sens
except ImportError:  # importing from a different directory
    from algo_sensitivity import Sens
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
import matplotlib.pyplot as plt
import stlabutils
import math


def func_cmpx_to_dBm(cmpx, kext, gain):
    Smpx = np.sqrt(kext)*cmpx  # input-output formalism
    return func_Smpx_to_dBm(Smpx, gain)


def func_Smpx_to_dBm(Smpx, gain):
    Smpx_W = abs(Smpx)**2/(2*50)  # P=V**2/(2*Z0)
    Smpx_dBm = 10*np.log10(Smpx_W)+30+gain  # P from W to dBm, plus gain
    return Smpx_dBm


def func_PdBmtoV(PindBm, Z0=50):
    # converts powers in dBm to Volts
    return np.sqrt(2*Z0*10**((PindBm-30)/10))


def algo_Ppeak(Vpeak, Z0=50):
    # converts volts to powers (both linear)
    return abs(Vpeak)**2/2/Z0


def algo_Vpeak0(V0, S11):
    # main drive: just the reflected voltage
    return V0 * np.sqrt(2*np.pi) * S11


def algo_Vpeak1(V0, I1, dSdf, dfdI, corr=2*np.sqrt(2*np.pi)*2):
    # first order side peak
    # FIXME: The correction factor is for testing purposes only!!
    return V0 * np.sqrt(2*np.pi) * I1 * dSdf * dfdI / corr


def algo_Vpeak2(V0, I1, dSdf, d2Sdf2, dfdI, d2fdI2, corr=2*np.sqrt(2*np.pi)*(2**2)):
    # second order side peak
    # FIXME: The correction factor is for testing purposes only!!
    return V0 * np.sqrt(2*np.pi) * I1**2 / 2 * (d2Sdf2 * dfdI**2 + dSdf * d2fdI2) / corr


def algo_Vpeak3(V0, I1, dSdf, d2Sdf2, d3Sdf3, dfdI, d2fdI2, d3fdI3, corr=2*np.sqrt(2*np.pi)*(2**3)):
    # third order side peak
    # FIXME: The correction factor is for testing purposes only!!
    return V0 * np.sqrt(2*np.pi) * 1 / 6 * I1**3 * (d3Sdf3 * dfdI**3 + 3 * d2Sdf2 * dfdI * d2fdI2 + dSdf * d3fdI3) / corr


def algo_dnSdfn(kext, kint, Delta, n):
    ktot = kext + kint
    return (1j)**n * 2**(n+1) * math.factorial(n) * kext / (ktot - 2j*Delta)**(n+1)


# Extracting resonance fequency and S11 line by line from raw dataset (megameasurement-A)
def S11_fitting_full(data, fitwidth=None, trimwidth=5, ftype='A', doplots=False, margin=51, xvariable='Is (A)', savefits=False):
    # #************************************************
    # fitwidth = None  # Width of trace to fit in units of detected peak width
    # trimwidth = 5.  # Width to trim around detected peak for background fit in units of detected peak width
    # ftype = 'A'  # A and -A are reflection models while B and -B are side coupled models
    # doplots = False  # Show debugging plots
    # margin = 51  # smoothing margin for peak detection
    # #************************************************

    # f0sfit = []
    # i0s = []
    # Qint = []
    # Qext = []
    # theta = []

    mydfs = []

    S11res_list = []

    for myline in data:
        xval = myline[xvariable][0]
        # i0s.append(i0)
        freqs = myline['Frequency (Hz)']
        S11complex = np.asarray(
            [a + 1j * b for a, b in zip(myline['S21re ()'], myline['S21im ()'])])  # S11 measurement
        # Do fit with some given parameters.  More options available.
        params, _, _, _ = stlabutils.S11fit(
            freqs, S11complex, ftype=ftype, doplots=doplots, trimwidth=trimwidth, fitwidth=fitwidth, margin=margin)
        fitresnobg = stlabutils.utils.S11fit.S11theo(
            freqs, params, ftype=ftype)  # S11 remove background
        fitres = stlabutils.S11func(
            freqs, params, ftype=ftype)  # S11 include background

        f0 = params['f0'].value

        # Part1: S11(f) vs I0 (xval)
        mydict = {'Frequency (Hz)': freqs, xvariable: xval, 'f0 (Hz)': f0,
                  'S11 ()': S11complex, 'S11phase (rad)': np.angle(S11complex), 'S11re ()': S11complex.real, 'S11im ()': S11complex.imag, 'S11abs ()': np.abs(S11complex),
                  'S11fit ()': fitres, 'S11fitphase (rad)': np.angle(fitres), 'S11fitre ()': fitres.real, 'S11fitim ()': fitres.imag, 'S11fitabs ()': np.abs(fitres),
                  'S11fitnobg ()': fitresnobg, 'S11fitnobgphase (rad)': np.angle(fitresnobg), 'S11fitnobgre ()': fitresnobg.real, 'S11fitnobgim ()': fitresnobg.imag, 'S11fitnobgabs ()': np.abs(fitresnobg)}
        mydfs.append(pd.DataFrame(mydict))

        if savefits:
            plt.plot(mydict['Frequency (Hz)'], mydict['S11abs ()'])
            plt.plot(mydict['Frequency (Hz)'], mydict['S11fitabs ()'])
            plt.title(xval)
            plt.show()
            plt.close()

        # Part2: Evaluating S11 on resonance
        try:
            S11_int = interp1d(freqs, S11complex)
            S11fit_int = interp1d(freqs, fitres)
            S11fitnobg_int = interp1d(freqs, fitresnobg)

            S11res_list.append({xvariable: xval, "f0 (Hz)": f0, "S11res ()": complex(S11_int(f0)),
                                "S11res_fit ()": complex(S11fit_int(f0)),
                                "S11res_fit_nobg ()": complex(S11fitnobg_int(f0)), 'fitpars': params})

        except ValueError:  # dirty trick to not having to deal with interpolation out of range
            S11_int = InterpolatedUnivariateSpline(freqs, S11complex)
            S11fit_int = InterpolatedUnivariateSpline(freqs, fitres)
            S11fitnobg_int = InterpolatedUnivariateSpline(freqs, fitresnobg)

            S11res_list.append({xvariable: xval, "f0 (Hz)": f0, "S11res ()": complex(S11_int(f0)),
                                "S11res_fit ()": complex(S11fit_int(f0)),
                                "S11res_fit_nobg ()": complex(S11fitnobg_int(f0)), 'fitpars': params})

    ds = S11res_list
    d = {}
    for k in ds[0].keys():
        d[k] = tuple(d[k] for d in ds)
    S11res_vs_I = pd.DataFrame(d)

    return S11res_vs_I, mydfs


def calculate_dSdf(f0_of_I, S11_of_I, key0='S11', unit0='()', key1='dSdf', unit1='(1/Hz)'):
    # calculates nth derivative of S11 w.r.t. bias current
    dSdf_vs_Is = []
    dSdfres_list = []

    for f0, myline in zip(f0_of_I, S11_of_I):
        i0 = myline['Is (A)'][0]
        freqs = myline['Frequency (Hz)']

        # Calculating the frequency-dependent derivative
        dSdf = np.gradient(myline[key0+' '+unit0], freqs)
        dSdffit = np.gradient(myline[key0+'fit '+unit0], freqs)
        dSdffitnobg = np.gradient(myline[key0+'fitnobg '+unit0], freqs)
        mydict = {'Frequency (Hz)': freqs, 'Is (A)': i0,
                  key1+' '+unit1: dSdf, key1+'re '+unit1: dSdf.real, key1+'im '+unit1: dSdf.imag, key1+'abs '+unit1: np.abs(dSdf),
                  key1+'fit '+unit1: dSdffit, key1+'fitre '+unit1: dSdffit.real, key1+'fitim '+unit1: dSdffit.imag, key1+'fitabs '+unit1: np.abs(dSdffit),
                  key1+'fitnobg '+unit1: dSdffitnobg, key1+'fitnobgre '+unit1: dSdffitnobg.real, key1+'fitnobgim '+unit1: dSdffitnobg.imag, key1+'fitnobgabs '+unit1: np.abs(dSdffitnobg)}
        dSdf_vs_Is.append(pd.DataFrame(mydict))

        # Evaluating the derivative on resonance
        dSdf_int = interp1d(freqs, dSdf)
        dSdffit_int = interp1d(freqs, dSdffit)
        dSdffitnobg_int = interp1d(freqs, dSdffitnobg)

        dSdfres_list.append({"Is (A)": i0, key1+'_data '+unit1: complex(dSdf_int(f0)),
                             key1+'_fit '+unit1: complex(dSdffit_int(f0)),
                             key1+'_fit_nobg '+unit1: complex(dSdffitnobg_int(f0))})

        ds = dSdfres_list
        d = {}
        for k in ds[0].keys():
            d[k] = tuple(d[k] for d in ds)
        dSdfres_vs_Is = pd.DataFrame(d)

    return dSdfres_vs_Is, dSdf_vs_Is

# Old and deprecated
# def calculate_first_sideband(V0, I1, dSdf, dfdI, Is_dSdf=[], Is_dfdI=[], testplot=True):

#     dSdf_int_real = interp1d(Is_dSdf, dSdf.values.real)
#     dSdf_int_imag = interp1d(Is_dSdf, dSdf.values.imag)
#     dfdI_int = interp1d(Is_dfdI, dfdI)

#     theo_currs = np.linspace(max([min(Is_dSdf), min(Is_dfdI)]), min(
#         [max(Is_dSdf), max(Is_dfdI)]), 401)

#     Vpeak = algo_Vpeak1(V0=V0, I1=I1, dSdf=dSdf_int_real(
#         theo_currs)+1j*dSdf_int_imag(theo_currs), dfdI=dfdI_int(theo_currs))
#     Ppeak = algo_Ppeak(Vpeak)

#     if testplot:

#         _, ax1 = plt.subplots()
#         plt.plot(Is_dSdf, abs(dSdf), '.-')
#         plt.ylabel('dS/df', color='C0')
#         plt.xlabel('Is (A)')
#         ax1.tick_params(axis='y', labelcolor='C0')
#         ax2 = ax1.twinx()
#         ax2.plot(Is_dfdI, dfdI, '.-C1')
#         plt.ylabel('df/dI', color='C1')
#         ax2.tick_params(axis='y', labelcolor='C1')
#         plt.title('Measured values')
#         plt.show()
#         plt.close()

#         _, ax1 = plt.subplots()
#         plt.plot(theo_currs, abs(dSdf_int_real(theo_currs) +
#                                  1j*dSdf_int_imag(theo_currs)), '.')
#         plt.ylabel('dS/df')
#         ax2 = ax1.twinx()
#         ax2.plot(theo_currs, dfdI_int(theo_currs), '.C1')
#         plt.ylabel('df/dI')
#         plt.xlabel('Is (A)')
#         plt.title('Interpolated values')
#         plt.show()
#         plt.close()

#         _, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
#         plt.sca(ax1)
#         plt.plot(theo_currs/1e-6, abs(Vpeak), label='Vpeak (V)')
#         plt.legend()
#         plt.sca(ax2)
#         plt.plot(theo_currs/1e-6, 10*np.log10(Ppeak*1000), label='Ppeak (dBm)')
#         plt.legend()
#         plt.suptitle('Expected peak height from theory')
#         plt.show()
#         plt.close()

#     return theo_currs, Ppeak


def calculate_nth_sideband(V0, I1, dSdf_dict, dfdI_dict, Is_dSdf_dict, Is_dfdI_dict, testplot=True, order='1'):

    dfdI = dfdI_dict['1']
    dSdf = dSdf_dict['1']
    Is_dSdf = Is_dSdf_dict['1']
    Is_dfdI = Is_dfdI_dict['1']
    dSdf_int_real = interp1d(Is_dSdf, dSdf.values.real)
    dSdf_int_imag = interp1d(Is_dSdf, dSdf.values.imag)
    dfdI_int = interp1d(Is_dfdI, dfdI)

    theo_currs1 = np.linspace(max([min(Is_dSdf), min(Is_dfdI)]), min(
        [max(Is_dSdf), max(Is_dfdI)]), 401)
    theo_currs = theo_currs1

    if len(dSdf_dict) == 1:
        Vpeak = algo_Vpeak1(V0=V0, I1=I1, dSdf=dSdf_int_real(
            theo_currs1)+1j*dSdf_int_imag(theo_currs1), dfdI=dfdI_int(theo_currs1))

    elif len(dSdf_dict) > 1:
        d2fdI2 = dfdI_dict['2']
        d2Sdf2 = dSdf_dict['2']
        Is_d2Sdf2 = Is_dSdf_dict['2']
        Is_d2fdI2 = Is_dfdI_dict['2']
        d2Sdf2_int_real = interp1d(Is_d2Sdf2, d2Sdf2.values.real)
        d2Sdf2_int_imag = interp1d(Is_d2Sdf2, d2Sdf2.values.imag)
        d2fdI2_int = interp1d(Is_d2fdI2, d2fdI2)

        theo_currs2 = np.linspace(max([min(theo_currs), min(Is_d2Sdf2), min(
            Is_d2fdI2)]), min([max(Is_d2Sdf2), max(Is_d2fdI2), max(theo_currs)]), 401)
        theo_currs = theo_currs2

        if len(dSdf_dict) == 2:
            Vpeak = algo_Vpeak2(V0=V0, I1=I1,
                                dSdf=dSdf_int_real(theo_currs2) +
                                1j*dSdf_int_imag(theo_currs2),
                                dfdI=dfdI_int(theo_currs2),
                                d2Sdf2=d2Sdf2_int_real(
                                    theo_currs2)+1j*d2Sdf2_int_imag(theo_currs2),
                                d2fdI2=d2fdI2_int(theo_currs2))

        elif len(dSdf_dict) == 3:
            d3fdI3 = dfdI_dict['3']
            d3Sdf3 = dSdf_dict['3']
            Is_d3Sdf3 = Is_dSdf_dict['3']
            Is_d3fdI3 = Is_dfdI_dict['3']
            d3Sdf3_int_real = interp1d(Is_d3Sdf3, d3Sdf3.values.real)
            d3Sdf3_int_imag = interp1d(Is_d3Sdf3, d3Sdf3.values.imag)
            d3fdI3_int = interp1d(Is_d3fdI3, d3fdI3)

            theo_currs3 = np.linspace(max([min(theo_currs2), min(Is_d3Sdf3), min(
                Is_d3fdI3)]), min([max(theo_currs2), max(Is_d3Sdf3), max(Is_d3fdI3)]), 401)
            theo_currs = theo_currs3

            Vpeak = algo_Vpeak3(V0=V0, I1=I1,
                                d3Sdf3=d3Sdf3_int_real(
                                    theo_currs3)+1j*d3Sdf3_int_imag(theo_currs3),
                                dfdI=dfdI_int(theo_currs3),
                                d2Sdf2=d2Sdf2_int_real(
                                    theo_currs3)+1j*d2Sdf2_int_imag(theo_currs3),
                                d2fdI2=d2fdI2_int(theo_currs3),
                                dSdf=dSdf_int_real(
                                    theo_currs3)+1j*dSdf_int_imag(theo_currs3),
                                d3fdI3=d3fdI3_int(theo_currs3))
        else:
            raise ValueError('Only works up to third order')

    Ppeak = algo_Ppeak(Vpeak)

    if testplot:

        _, ax1 = plt.subplots()
        plt.plot(Is_dSdf_dict[order], abs(dSdf_dict[order]), '.-')
        plt.ylabel('dS/df', color='C0')
        plt.xlabel('Is (A)')
        ax1.tick_params(axis='y', labelcolor='C0')
        ax2 = ax1.twinx()
        ax2.plot(Is_dfdI_dict[order], dfdI_dict[order], '.-C1')
        plt.ylabel('df/dI', color='C1')
        ax2.tick_params(axis='y', labelcolor='C1')
        plt.title('Measured values')
        plt.show()
        plt.close()

        # _, ax1 = plt.subplots()
        # plt.plot(theo_currs, abs(dSdf_int_real(theo_currs) +
        #                          1j*dSdf_int_imag(theo_currs)), '.')
        # plt.ylabel('dS/df')
        # ax2 = ax1.twinx()
        # ax2.plot(theo_currs, dfdI_int(theo_currs), '.C1')
        # plt.ylabel('df/dI')
        # plt.xlabel('Is (A)')
        # plt.title('Interpolated values')
        # plt.show()
        # plt.close()

        _, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
        plt.sca(ax1)
        plt.plot(theo_currs/1e-6, abs(Vpeak), label='Vpeak (V)')
        plt.legend()
        plt.sca(ax2)
        plt.plot(theo_currs/1e-6, 10*np.log10(Ppeak*1000), label='Ppeak (dBm)')
        plt.legend()
        plt.suptitle('Expected peak height from theory')
        plt.show()
        plt.close()

    return theo_currs, Ppeak


def calculate_expected(meas_data, theo_currs, Ppeak_theo, amplification):

    # fixed values
    exp_bg = meas_data['Background (dBm)'].iloc[0]
    exp_RBW = meas_data['RBW (Hz)'].iloc[0]
    exp_Imodpp = meas_data['Imodpp (A)'].iloc[0]

    theo_signal = 10*np.log10(Ppeak_theo*1000)
    amplification = 100  # careful, this one has a big influence! exact value?
    theo_spectrum = theo_signal + amplification
    theo_sensitivity = Sens(exp_Imodpp, theo_spectrum, exp_bg, exp_RBW)

    # sometimes current suddenly turns complex?
    return abs(theo_currs), theo_spectrum, theo_sensitivity


def calculate_measurement(meas_data, order=1, xkey='Is (A)'):

    # fixed values
    exp_bg = meas_data['Background (dBm)'].iloc[0]
    exp_RBW = meas_data['RBW (Hz)'].iloc[0]
    exp_Imodpp = meas_data['Imodpp (A)'].iloc[0]

    # FIXME: sometimes current suddenly gets a zero imaginary component?
    try:
        meas_curr = meas_data[xkey].real
    except KeyError:
        # if xvalues were assigned to DataFrame as index
        meas_curr = meas_data.axes[0].values
    meas_signal_px = meas_data['f0+{}f1 (dBm)'.format(int(order))]
    meas_sensitivity_px = Sens(exp_Imodpp, meas_signal_px, exp_bg, exp_RBW)

    if order > 0:
        meas_signal_mx = meas_data['f0-{}f1 (dBm)'.format(int(order))]
        meas_sensitivity_mx = Sens(exp_Imodpp, meas_signal_mx, exp_bg, exp_RBW)

        return meas_curr, meas_signal_px, meas_sensitivity_px, meas_signal_mx, meas_sensitivity_mx

    else:
        return meas_curr, meas_signal_px, meas_sensitivity_px
