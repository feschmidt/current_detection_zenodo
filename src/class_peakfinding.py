# This is the class to automate the peakfinding for multiple blocks of data

import numpy as np
import pandas as pd
from .algo_peakfinding import GetPeaks, GetBackground


class Peakfinder():

    def __init__(self, classparams):

        self.xkey = 'Is (A)'
        self.ykey = 'Carrier Power (dBm)'
        self.speckey = 'Spectrum (dBm)'  # PSD or Spectrum
        self.RBW = 5
        self.df = 100
        self.filename = False
        self.nrange = np.arange(-3, 4)
        self.bgmethod = 'first'
        self.RBWconst = True
        for key, val in classparams.items():
            setattr(self, key, val)
        self.nrange = np.sort(self.nrange)  # to be sure, sort from low to high

        pass

    def GetFreqsAndPeaks(self, mydata, figshow=False, figsave=False, xval=np.nan, setdF=None):
        """
        input: spectrum analyzer trace
        output: dataframe with {xkey,freqs,peakvals,background} 
        """

        xvals, yvals = [], []
        Imodpp, Modfreq, RBW = [], [], []
        background = []
        temp_freqs = [[] for _ in self.nrange]
        temp_signs = [[] for _ in self.nrange]
        bgkey = True

        for i, block in enumerate(mydata):
            if not self.filename:
                xvals.append(block[self.xkey].iloc[0])
            else:
                xvals.append(xval)
            if setdF != None:
                block[self.ykey] = setdF[i]
            yvals.append(block[self.ykey].iloc[0])
            Imodpp.append(block['Imodpp (A)'].iloc[0])
            f1 = block['Modfrec (Hz)'].iloc[0]
            Modfreq.append(f1)
            try:
                RBW.append(block['RBW (Hz)'].iloc[0])
            except KeyError:
                # Note: The RBW was mantained constant during all measurements as can be seen in the measurement scripts
                RBW.append(5)
            try:
                spec = block[self.speckey]
            except KeyError:
                spec = block['PSD (dBm)']
            freqs = block['Frequency (Hz)']

            if figsave:
                thefigsave = figsave+'{:03d}'.format(i)
            else:
                thefigsave = False
            result = GetPeaks(freqs, spec, f1=f1,
                              df=self.df, figshow=figshow, figsave=thefigsave, nrange=self.nrange, bgmethod=self.bgmethod)

            # background is the same values for all peaks. It just gets returned once for each peak because of the dataframe structure
            bgkey = result['bgkey'].iloc[0]  # True or False
            thebg = result['background'].iloc[0]  # actual value
            if bgkey:
                background.append(thebg)
            else:
                if self.RBWconst == True:
                    # does the RBW change during the measurement? This will lead to changing backgrounds
                    if i == 0:
                        background.append(thebg)
                    else:
                        background.append(background[i-1])
                else:
                    # does the RBW change during the measurement? This will lead to changing backgrounds
                    background.append(thebg)

            fs = result['peakfreqs']  # frequency values
            ps = result['peakvals']  # peak values

            for i, (thef, thep) in enumerate(zip(fs, ps)):
                temp_freqs[i].append(thef)
                temp_signs[i].append(thep)

        if self.speckey == 'PSD (dBm)':
            self.speckey = 'Spectrum (dBm)'  # old files have the wrong label
        mydict = {self.xkey: xvals,
                  self.ykey: yvals,
                  'Imodpp (A)': Imodpp,
                  'Modfrec (Hz)': Modfreq,
                  'RBW (Hz)': RBW,
                  'Background (dBm)': background}
        for i, thekey in enumerate(result['keys_freq']):
            mydict[thekey+' (Hz)'] = temp_freqs[i]
            mydict[thekey+' (dBm)'] = temp_signs[i]

        mydf = pd.DataFrame(mydict)

        return mydf


if __name__ == "__main__":

    import stlabutils
    import glob
    import pandas as pd

    megafiles = glob.glob(
        '../../../../measurement_data/BlueFors/Felix/190222_DCbias_SQUID_2ndgen_currentbias/Switch_position_2/F31_megameasurement2/F31B_2019_04_01_19.22.24_SA_vs_Is_0to9uA_500mVpp/*.dat')
    myfile = megafiles[0]

    mydata = stlabutils.readdata.readdat(myfile)
    myFinder = Peakfinder({})

    thedf = myFinder.GetFreqsAndPeaks(mydata)
    print(thedf.head())
