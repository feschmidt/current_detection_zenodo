"""
Model of current bias cavity
"""
import numpy as np
import scipy.constants as sciconst
from scipy.interpolate import interp1d

echarge = sciconst.e
Phi0 = sciconst.Planck/(2*echarge)
Rk = sciconst.Planck/(echarge**2)
pi = np.pi


def lossrate(x, k0, alpha, ix):
    # exponential dependence of loss rates on bias current
    return k0 + alpha*np.exp(x/ix)


def Lj(I, Ic):
    # Josephson inductance assuming sinusoidal CPR
    return Phi0 / (2*pi*np.sqrt(Ic**2 - I**2))


def f0(I, fr, Lr, Ic, nJJ=1):
    # each JJ shifts the resonance frequency further downwards
    return fr / (1 + nJJ*Lj(I, Ic) / Lr)


def f0leads(I, fr, Lr, Ic, nJJ=1, Lg=27e-12):
    # each JJ shifts the resonance frequency further downwards
    # assuming additional lead inductance
    return fr / (1 + (nJJ*Lj(I, Ic)+Lg) / Lr)


def f0Lk(I, fk0, alpha, Istar):
    # absolute resonance frequency shift due to partial kinetic inductance
    # alpha = Lk/(Lk+Lg)
    return fk0 / np.sqrt(1 + alpha * I**2 / Istar**2)


def f0Lkshift(I, I1, I2=1):
    # relative frequency shift due to only kinetic inductance
    return -(1/2*(I/I1)**2+1/4*(I/I2)**4)


def Xpp(vph, Z0, l, Lj, w0):
    # source: general model, originating from Black Box Quantization (BBQ)
    numerator = 2*echarge**2 * Lj * w0**2
    denominator = (1 + Lj*w0**2 * l/(vph*Z0*np.sin(l*w0/vph)**2))**2
    return numerator/denominator


def Duffing_parameter(I, fr, Lr, Ic, Z0=50, nJJ=1):
    # source: Pogarzalek
    wr = 2*np.pi*fr
    return (np.pi**2 * wr * Z0 / Rk) * (nJJ*Lj(I, Ic) / Lr)**3


def Duffing_parameter_2(I, fr, Lr, Ic, Z0=50, nJJ=1):
    # source: Pogarzalek, but with participation ratio
    wr = 2*np.pi*fr
    return (np.pi**2 * wr * Z0 / Rk) * (nJJ*Lj(I, Ic) / (Lr + nJJ*Lj(I, Ic)))**3


def G1(I, fr, Lr, Ic, nJJ=1):
    if type(I) == float or type(
            I) == int:  # single bias current, need to evaluate here
        ii = np.linspace(-Ic, Ic, 1001)
        return interp1d(ii, np.gradient(f0(ii, fr, Lr, Ic, nJJ), ii))(I)
    else:  # array of bias currents
        return np.gradient(f0(I, fr, Lr, Ic, nJJ), I)


def G2(I, fr, Lr, Ic, nJJ=1):
    # only works for array of bias currents
    return np.gradient(G1(I, fr, Lr, Ic, nJJ), I)


def G3(I, fr, Lr, Ic, nJJ=1):
    # only works for array of bias currents
    return np.gradient(G2(I, fr, Lr, Ic, nJJ), I)


"""
# old and deprecated functions

def Duffing_parameter_array(I, fr, Lr, Ic, nJJ=1):
    # source: Master thesis BADWMI
    w0 = 2*np.pi*f0(I, fr, Lr, Ic, nJJ)
    p = nJJ*Lj(I, Ic) / (nJJ*Lj(I, Ic)+Lr)
    return -w0**2*p**3/(2*nJJ**2)


def alpha0(ke, ki, Sin, Delta):
    return -2 * np.sqrt(ke) * Sin / (ki + ke + 2j * Delta)


def S1out(I, fr, Lr, Ic, ke, ki, Sin, Delta, Omega, nJJ, Iac, calc_alpha0=False):
    if type(calc_alpha0) == bool and calc_alpha0 == False:  # low-power regime
        S1outplus = np.sqrt(ke) * alpha0(ke, ki, Sin, Delta) * G1(
            I, fr, Lr, Ic, nJJ) * Iac / (-1j*(ke + ki) + 2 * (Delta + Omega))
        S1outminus = np.sqrt(ke) * alpha0(ke, ki, Sin, Delta) * G1(
            I, fr, Lr, Ic, nJJ) * Iac / (-1j*(ke + ki) + 2 * (Delta - Omega))
    else:  # Duffing regime
        S1outplus = np.sqrt(ke) * calc_alpha0 * G1(
            I, fr, Lr, Ic, nJJ) * Iac / (-1j*(ke + ki) + 2 * (Delta + Omega))
        S1outminus = np.sqrt(ke) * calc_alpha0 * G1(
            I, fr, Lr, Ic, nJJ) * Iac / (-1j*(ke + ki) + 2 * (Delta - Omega))
    return S1outplus, S1outminus
"""
