# This is the algorithm to read in raw data and extract peaks

import glob
import pandas as pd
import stlabutils
from .class_peakfinding import Peakfinder
from .class_sensitivity import SensCalculator


class FileRunner():

    def __init__(self, classparams):
        self.pos = 'pos2'
        self.dev = 'singleJJ'
        for key, val in classparams.items():
            setattr(self, key, val)
        pass

    # change xkey to e.g. 'Carrier Power (dBm)' to get a `raw2peaks_power_dF` function
    def raw2peaks_current_dF(self, megafiles, figsave=False, savedI=True, xkey='Is (A)', saveddF=(True, [])):
        for myfile in megafiles:
            fileid = myfile.split('/')[-2]
            print(fileid)
            mydata = stlabutils.readdata.readdat(myfile)
            if savedI:
                myFinder = Peakfinder({'xkey': xkey, 'ykey': 'dF (Hz)'})
                if figsave:
                    figsave = 'data_processed/'+self.pos+'/.testplots/'+self.dev+fileid
                if saveddF[0] == True:
                    thedf = myFinder.GetFreqsAndPeaks(mydata, figsave=figsave)
                # if the datafile doesn't contain the `dF (Hz)` column, we can manually enter it here from looking at the measurement file
                else:
                    thedf = myFinder.GetFreqsAndPeaks(
                        mydata, figsave=figsave, setdF=saveddF[1])
            else:
                # need to get the current from filename
                isetnA = float(myfile.split('_')[-1].strip('.dat'))
                Iset = isetnA*1e-9
                myFinder = Peakfinder(
                    {'xkey': 'Is (A)', 'ykey': 'dF (Hz)', 'filename': True})
                if figsave:
                    figsave = 'data_processed/'+self.pos+'/.testplots/'+self.dev+fileid
                thedf = myFinder.GetFreqsAndPeaks(
                    mydata, xval=Iset, figsave=figsave)
            thedf.to_pickle('data_processed/'+self.pos+'/'+fileid+'.pkl')

    def raw2peaks_current_pwr(self, megafiles, figsave=False, xkey='Is (A)'):
        for myfile in megafiles:
            fileid = myfile.split('/')[-2]
            print(fileid)
            mydata = stlabutils.readdata.readdat(myfile)
            myFinder = Peakfinder({'xkey': xkey})
            if figsave:
                figsave = 'data_processed/'+self.pos+'/.testplots/'+self.dev+fileid
            thedf = myFinder.GetFreqsAndPeaks(
                mydata, figsave=figsave)
            thedf.to_pickle('data_processed/'+self.pos+'/'+fileid+'.pkl')

    def peaks2sens_current(self, pklfiles, figsave=True, xkey='dF'):
        for myfile in pklfiles:
            fileid = myfile.split('/')[-1].strip('.pkl')
            fileid = 'sensitivity_'+fileid
            print(fileid)
            mydata = pd.read_pickle(myfile)
            if xkey == 'dF':
                mySens = SensCalculator({'xkey': 'dF (Hz)'})
            else:
                mySens = SensCalculator({})
            if figsave:
                figsave = 'data_sensitivity/'+self.pos+'/.testplots /'+fileid
            sensdf = mySens.calculate(
                mydata, figsave=figsave)
            sensdf.to_pickle('data_sensitivity/'+self.pos+'/'+fileid+'.pkl')
