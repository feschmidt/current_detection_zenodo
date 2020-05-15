# This is the algorithm for finding the peaks

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def GetPeaks(freqs, spec, f1=1e3, df=100, figshow=False, figsave=False, nrange=np.arange(-3, 4), bgmethod='first'):
    """
    Load in data, return the frequencies and peakvalues and left over frequencies and background
    """

    f0 = freqs.iloc[len(freqs)//2]  # freqs[np.argmax(freqs)]
    dfs = []
    nrange = np.sort(nrange)  # to be sure, sort from low to high
    keys_freq = []

    # maybe here implement that the ranges are found automatically, in case the span is varied...
    for n in nrange:
        fnew = f0+n*f1
        dfs.append([fnew-df, fnew+df])
        keys_freq.append('f0{0:+d}f1'.format(n))

    if figshow:
        plt.plot(freqs, spec)
        [plt.axvspan(xmin=df[0], xmax=df[1], alpha=0.5, color='orange')
         for df in dfs]
        plt.grid(which='both')
        plt.title('Data and intervals')
        plt.show()
        plt.close()

    # find the indexes of frequencies outside the peak intervals
    indexes = []
    for thedf in dfs:
        idx = np.logical_not((freqs >= thedf[0]) & (freqs <= thedf[1]))
        indexes.append(idx)

    fullidx = np.logical_and.reduce(indexes)
    # notidx = np.logical_not(fullidx)

    # find the peaks as maximum values in each interval
    peakfreqs = []
    peakvals = []
    for idx in indexes:
        idx = np.logical_not(idx)
        newfreq = freqs[idx].values
        newspec = spec[idx].values
        mxdx = np.argmax(newspec)
        peakfreqs.append(newfreq[mxdx])
        peakvals.append(max(newspec))

    # Get the background
    background = GetBackground(spec[fullidx], bgmethod)
    if ((figshow) or (figsave != False)):
        plt.plot(freqs, spec)
        plt.plot(freqs[fullidx], spec[fullidx], '.')
        [plt.axvspan(xmin=df[0], xmax=df[1], alpha=0.5, color='orange')
         for df in dfs]
        plt.plot(peakfreqs, peakvals, 'o-')
        plt.grid(which='both')
        plt.title('Everything')
        if figshow:
            plt.show()
        if figsave != False:
            plt.savefig(figsave+'.png')
        plt.close()

    if figshow:
        plt.plot(freqs[fullidx], spec[fullidx], '.')
        plt.axhline(background[0], c='k')
        [plt.axvspan(xmin=df[0], xmax=df[1], alpha=0.5, color='orange')
         for df in dfs]
        plt.grid(which='both')
        plt.title('Background only')
        plt.show()
        plt.close()

    mydf = pd.DataFrame({'peakfreqs': peakfreqs, 'peakvals': peakvals,
                         'background': background[0], 'bgkey': background[1], 'keys_freq': keys_freq})

    return mydf


def GetBackground(spec, bgmethod='first', calcmethod='mean'):
    """
    Calculates the background using either the mean ('mean') or max value ('max') of the given spectrum ('first','mean') or by loading in a reference spectrum ('refdata').
    Returns (background,Key) where the Key is either True or False depending on whether the next time the bgmethod gets called should repeat the calculation or not.
    """

    if bgmethod == 'first':
        bgkey = False
    elif bgmethod == 'refdata':
        # TODO: to be implemented
        raise KeyError('refdata not yet implemented')
    elif bgmethod == 'mean':
        bgkey = True
    else:
        raise KeyError('Unknown bgmethod: '+bgmethod)

    if calcmethod == 'mean':
        thebg = mW2dBm(np.mean(dBm2mW(spec)))
    elif calcmethod == 'max':
        thebg = mW2dBm(np.max(dBm2mW(spec)))
    else:
        raise KeyError('Unknown calcmethod: '+calcmethod)

    return (thebg, bgkey)


def dBm2mW(dBmpow):
    return 10**(dBmpow/10)


def mW2dBm(mWpow):
    return 10*np.log10(mWpow)


if __name__ == "__main__":

    import stlabutils
    import glob
    import pandas as pd

    megafiles = glob.glob(
        '../../../../measurement_data/BlueFors/Felix/190222_DCbias_SQUID_2ndgen_currentbias/Switch_position_2/F31_megameasurement2/F31B_2019_04_01_19.22.24_SA_vs_Is_0to9uA_500mVpp/*.dat')
    myfile = megafiles[0]

    mydata = stlabutils.readdata.readdat(myfile)
    block = mydata[-1]

    spec = block['PSD (dBm)']
    freqs = block['Frequency (Hz)']

    result = GetPeaks(freqs, spec, f1=1e3, df=100, figshow=True)
