# This is the class to automate the peakfinding for multiple blocks of data

import pandas as pd
from .algo_sensitivity import Sens, ENBW
import matplotlib.pyplot as plt
import numpy as np
import re
import time
# to suppress SettingWithCopyWarning for placing background in params
pd.options.mode.chained_assignment = None


class SensCalculator():

    def __init__(self, classparams):
        self.xkey = 'Carrier Power (dBm)'
        for key, val in classparams.items():
            setattr(self, key, val)

        pass

    def calculate(self, mydata, figshow=False, figsave=False, change_rc=True):
        """
        input: processed dataframe with {params,freqs,peakvals,background}
        output: processed dataframe with {params,sensitivity}
        """

        params = mydata.loc[:, :'RBW (Hz)']
        rex = re.compile(r'f0..f1..dBm.')
        l = ",".join(list(mydata.axes[1]))
        columnnames = rex.findall(l)
        # WARNING: Iac should be amplitude, not peak-to-peak value!
        ipp = mydata['Imodpp (A)']/2
        background = mydata['Background (dBm)'].iloc[0]
        params['Background (dBm)'] = background
        RBW = mydata['RBW (Hz)']
        Iset = mydata['Is (A)'].iloc[0]
        mysens = {}
        if change_rc:
            plt.rcParams['figure.dpi'] = 160

        for i, column in enumerate(columnnames):
            sensitivity = Sens(ipp=ipp, signal=mydata[column],
                               background=background, RBW=RBW)
            mysens['Sensitivity ' +
                   column.replace('(dBm)', '(A/sqrt(Hz))')] = sensitivity  # store sensitivity
            mysens[column] = mydata[column]  # store amplitude
            if ((figshow) or (figsave != False)):
                plt.plot(mydata[self.xkey], sensitivity/1e-12, '.')
                plt.plot(mydata[self.xkey][sensitivity.idxmin()],
                         min(sensitivity)/1e-12, '*')
                plt.xlabel(self.xkey)
                plt.ylabel(r'Sensitivity ($pA/\sqrt{Hz}$)')
                plt.title(
                    column+', Is={:.3f}uA, min|S|={:.3f}'.format(Iset/1e-6, min(sensitivity/1e-12))+r'$pA /\sqrt{Hz}$')
                plt.ylim(0, 1000)
                plt.tight_layout()
                if figshow:
                    plt.show()
                # note: this will only save the figures with this column name!!
                if ((figsave != False) and ((column == 'f0-1f1 (dBm)') or (column == 'f0+1f1 (dBm)'))):
                    plt.savefig(
                        figsave+column+'{:03d}.png'.format(i), bbox_to_inches='tight')
                plt.close()

        # add the frequency values
        rex = re.compile(r'f0..f1..Hz.')
        l = ",".join(list(mydata.axes[1]))
        columnnames = rex.findall(l)
        for i, column in enumerate(columnnames):
            mysens[column] = mydata[column]  # store frequency

        processed_frame = pd.DataFrame(mysens)
        mydf = params.join(processed_frame)

        return mydf


if __name__ == "__main__":

    import stlabutils
    import glob
    import pandas as pd

    pklfiles = sorted(glob.glob('data_processed/pos2/F31E*.pkl'))
    pklfile = pklfiles[10]
    mydata = pd.read_pickle(pklfile)

    mySens = SensCalculator({})

    result = mySens.calculate(mydata, figshow=True)
    print(result.head())
