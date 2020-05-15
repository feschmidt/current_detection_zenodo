# This is the algorithm to calculate the sensitivity

import numpy as np


def Sens(ipp, signal, background, RBW):
    """
    Returns the sensitivity in A/sqrt(Hz), given an input current in A, signal and background in dBm, and RBW in Hz
    WARNING: make sure that the ipp is NOT the peak-to-peak value, but rather the amplitude iac=ipp/2
    """
    return ipp/(np.sqrt(ENBW(RBW) * 10**((signal-background)/10)))


def ENBW(RBW, ftype='Rauscher'):
    """
    returns the corrected ENBW from the RBW
    factor calculated: see `Spectrum analysis theory.ipynb`
    factor Rauscher see: Rauscher et al., Fundamentals of spectrum analysis
    TODO: why does the measured result differ from literature?
    """
    factors = {'Rauscher': 1.065, 'calculated': 1.0553117067566966}
    return factors[ftype]*RBW


if __name__ == "__main__":
    bg = -90  # in dBm
    RBW = 5  # in Hz
    imodpp = 500e-9  # in A
    sig = -40  # in dBm

    sensitivity = Sens(imodpp, sig, bg, RBW)
    print(sensitivity)
