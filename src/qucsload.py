import numpy as np


def loadQUCSparsweep(myfile):
    mydata = np.loadtxt(myfile, skiprows=1, delimiter=';')

    steps = sum(abs(np.diff(mydata[:, 1])) > 0)+1
    fulllen = mydata.shape[0]
    # numpars = mydata.shape[1]
    steplen = fulllen//steps

    freqs = mydata[:, 0].reshape((steps, steplen))
    currs = mydata[:, 1].reshape((steps, steplen))
    S11re = mydata[:, 2].reshape((steps, steplen))
    S11im = mydata[:, 3].reshape((steps, steplen))

    S11complex = S11re+1j*S11im
    S11dB = 20*np.log10(abs(S11complex))
    S11ph = np.angle(S11complex)

    myoutput = [freqs, currs, S11dB, S11ph]
    return myoutput


def loadQUCSsingle(myfile):
    mydata = np.loadtxt(myfile, skiprows=1, delimiter=';')

    steps = sum(abs(np.diff(mydata[:, 1])) > 0)+1
    fulllen = mydata.shape[0]
    # numpars = mydata.shape[1]
    steplen = fulllen//steps

    freqs = mydata[:, 0].reshape((steps, steplen))
    S11re = mydata[:, 1].reshape((steps, steplen))
    S11im = mydata[:, 2].reshape((steps, steplen))

    S11complex = S11re+1j*S11im
    S11dB = 20*np.log10(abs(S11complex))
    S11ph = np.angle(S11complex)

    myoutput = [freqs, S11dB, S11ph]
    return myoutput
