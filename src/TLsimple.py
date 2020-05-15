# Simple TL model for current bias cavity using lumped components

import numpy as np


def S11theo(w, w0, kint, kext):
    dw = w-w0
    return (kext-kint-2j*dw) / (kext+kint+2j*dw)


def w0(cs, c, l):
    return np.sqrt((2*c+cs)/(l*c*cs))


def Qint(w0, r, cs, c):
    return (c+cs)/(w0*r*c*cs)


def Qext(w0, Z0, cs, c):
    return w0*Z0*cs*(c+cs)/c


def kappa(w0, Q):
    return w0/Q


def Zin(cs, c, kint, w, w0):
    dw = w-w0
    return c/(cs*(c+cs))*1/(kint+2j*dw)
