# This is the algorithm to read in sensitivity files and plot them

import numpy as np
from scipy.interpolate import interp1d
import pandas as pd
import matplotlib.pyplot as plt
import stlabutils
import pickle


class FilePlotter():

    def __init__(self, classparams):
        self.pos = 'pos2'
        self.id = 'F33'
        self.xkey = 'Is (A)'
        self.xplotkey = 'Bias current (uA)'
        self.xscale = 1e-6
        self.ykey = 'dF (Hz)'
        self.thepeaks = ['f0-3f1', 'f0-2f1', 'f0-1f1',
                         'f0+0f1', 'f0+1f1', 'f0+2f1', 'f0+3f1']
        self.dpi = 160
        self.transp = True
        self.bbox = 'tight'
        for key, val in classparams.items():
            setattr(self, key, val)

        if self.ykey == 'dF (Hz)':
            self.yplotkey = 'Detuning (kHz)'
            self.yscale = 1e3
            self.yunit = 'kHz'
            self.key = 'dF'
        elif self.ykey == 'Carrier Power (dBm)':
            self.yplotkey = self.ykey
            self.yscale = 1
            self.yunit = 'dBm'
            self.key = 'pwr'
        else:
            raise KeyError('unknown ykey')

        self.figdir = 'plots/'+self.pos+'/'+self.id+'_'+self.key+'_'
        self.pkldir = 'data_plots/'+self.pos+'/'+self.id+'_'+self.key+'_'
        self.setrcparams(self.dpi, self.transp, self.bbox)
        pass

    def setrcparams(self, dpi=160, transp=True, bbox='tight'):
        plt.rcParams['figure.dpi'] = dpi
        plt.rcParams['savefig.transparent'] = transp
        plt.rcParams['savefig.bbox'] = bbox

    def runall_sens2plots(self, pklfiles, mycmap='RdBu', processlist=None, clim=None, figsave=False, figshow=True, pklsave=True):
        for whichpeak in self.thepeaks:
            self.sens2plots(pklfiles, whichpeak, mycmap,
                            processlist, clim, figsave, figshow, pklsave)

    def sens2plots(self, pklfiles, whichpeak, mycmap='RdBu', processlist=None, clim=None, figsave=False, figshow=True, pklsave=True):
        if self.ykey == 'dF (Hz)':
            mydfs = []
            for myfile in pklfiles:
                mydata = pd.read_pickle(myfile)
                yvalue = mydata[self.ykey]
                xvalue = mydata[self.xkey].iloc[0]
                thesensitivity = mydata['Sensitivity ' +
                                        whichpeak+' (A/sqrt(Hz))']
                mydict = {}
                mydict['Sensitivity '+whichpeak +
                       ' (pA/sqrt(Hz))'] = thesensitivity/1e-12
                mydf = pd.DataFrame(mydict)
                mydf[self.yplotkey] = yvalue/self.yscale
                mydf[self.xplotkey] = xvalue/self.xscale
                mydfs.append(mydf)
            # Here we don't need interpolation because the dFs are constant

        elif self.ykey == 'Carrier Power (dBm)':
            xx = []
            carrpow = []
            sens = []

            for myfile in pklfiles:
                mydata = pd.read_pickle(myfile)
                yvalue = mydata[self.ykey]
                xvalue = mydata[self.xkey].iloc[0]
                thesensitivity = mydata['Sensitivity ' +
                                        whichpeak+' (A/sqrt(Hz))']
                xx.append(xvalue)
                carrpow.append(yvalue)
                sens.append(thesensitivity)

            # We need interpolation because the carrier powers are slightly different from each other
            minpwr = min([min(x) for x in carrpow])
            maxpwr = max([max(x) for x in carrpow])
            newpwrs = np.linspace(minpwr, maxpwr, len(carrpow[0])*2)

            # funcs = []
            mydfs = []
            mydict = {}
            for (thepower, ii, thesensitivity) in zip(carrpow, xx, sens):
                func = interp1d(thepower, thesensitivity,
                                bounds_error=False)
                mydict['Sensitivity '+whichpeak +
                       ' (pA/sqrt(Hz))'] = func(newpwrs)/1e-12
                mydf = pd.DataFrame(mydict)
                mydf[self.yplotkey] = newpwrs/self.yscale
                mydf[self.xplotkey] = ii/self.xscale
                mydfs.append(mydf)

        try:
            mymtx = stlabutils.utils.stlabdict.framearr_to_mtx(
                mydfs, key='Sensitivity '+whichpeak+' (pA/sqrt(Hz))', xkey=self.xplotkey, ykey=self.yplotkey)
        except:
            # if script got interrupted, the last dataframe might be incomplete
            mymtx = stlabutils.utils.stlabdict.framearr_to_mtx(
                mydfs[:-1], key='Sensitivity '+whichpeak+' (pA/sqrt(Hz))', xkey=self.xplotkey, ykey=self.yplotkey)

        values = mymtx.pmtx.values
        minloc = np.unravel_index(np.nanargmin(values), values.shape)
        axes = mymtx.pmtx.axes
        minIs = axes[1][minloc[1]]
        minYY = axes[0][minloc[0]]
        minS = values[minloc]

        if processlist:
            mymtx.applyprocesslist(processlist)

        x = axes[1].values
        y = axes[0].values
        X, Y = np.meshgrid(x, y)
        Z = mymtx.pmtx.values
        # plt.pcolormesh(X, Y, Z, cmap=mycmap, clim=clim)
        plt.imshow(mymtx.pmtx, extent=mymtx.getextents(),
                   aspect='auto', cmap=mycmap, clim=clim)
        plt.plot(minIs, minYY, '*', c='red', ms=15)
        plt.xlabel(self.xplotkey)
        plt.ylabel(self.yplotkey)
        cbar = plt.colorbar()
        if any('log' in x for x in processlist):
            cbar.set_label(
                r'log(Sensitivity) ($\log(\mathrm{pA}/\sqrt{\mathrm{Hz}})$)')
            plt.title(r'log(Sensitivity) for '+whichpeak+r', min|S| {:.3f}'.format(
                minS)+r' $\mathrm{pA}/\sqrt{\mathrm{Hz}}$'+r' @ {:.1f} $\mu$A, {:.1f}'.format(minIs, minYY)+' '+self.yunit)
        else:
            cbar.set_label(
                r'Sensitivity ($\mathrm{pA}/\sqrt{\mathrm{Hz}}$)')
            plt.title(r'Sensitivity for '+whichpeak+r', min|S| {:.3f}'.format(
                minS)+r' $\mathrm{pA}/\sqrt{\mathrm{Hz}}$'+r' @ {:.1f} $\mu$A, {:.1f}'.format(minIs, minYY)+' '+self.yunit)
        if figsave:
            plt.savefig(self.figdir+whichpeak+'.png')
        if figshow:
            plt.show()
        plt.close()

        if pklsave:
            mymtx.pmtx.to_pickle(self.pkldir+whichpeak+'.pkl')
            pickle_out = open(self.pkldir+whichpeak+'_mtx.pkl', 'wb')
            pickle.dump(mymtx, pickle_out)
            pickle_out.close()

    def runall_sens2plots_1D(self, pklfile, logscale=False, figsave=False, figshow=True, pklsave=True):
        for whichpeak in self.thepeaks:
            self.sens2plots_1D(pklfile, whichpeak, logscale,
                               figsave, figshow, pklsave)

    def sens2plots_1D(self, pklfile, whichpeak, logscale=False, figsave=False, figshow=True, pklsave=True):
        # for plotting a single linecut
        print('here')
        # FIXME: weird error here, which I cannot reproduce outside of this structure
        mydata = pd.read_pickle(pklfile)
        print('here2')
        xvalues = mydata[self.ykey]
        thesensitivity = mydata['Sensitivity ' +
                                whichpeak+' (A/sqrt(Hz))']
        # Here we don't need interpolation

        minX = np.nanargmin(thesensitivity)
        minY = np.nanmin(thesensitivity)

        if logscale:
            plt.plot(xvalues, np.log10(thesensitivity))
            plt.title(r'log(Sensitivity) for '+whichpeak+r', min|S| {:.3f}'.format(
                minY)+r' $\mathrm{pA}/\sqrt{\mathrm{Hz}}$'+r' @ {:.1f}'.format(minX)+' '+self.yunit)
        else:
            plt.plot(xvalues, thesensitivity)
            plt.title(r'Sensitivity for '+whichpeak+r', min|S| {:.3f}'.format(
                minY)+r' $\mathrm{pA}/\sqrt{\mathrm{Hz}}$'+r' @ {:.1f}'.format(minX)+' '+self.yunit)
        plt.plot(minX, minY, '*', c='yellow', ms=15)
        plt.xlabel(self.xplotkey)
        plt.ylabel(self.yplotkey)
        if figsave:
            plt.savefig(self.figdir+whichpeak+'.png')
        if figshow:
            plt.show()
        plt.close()

        if pklsave:
            pickle_out = open(self.pkldir+whichpeak+'_linecut.pkl', 'wb')
            pickle.dump([xvalues, thesensitivity], pickle_out)
            pickle_out.close()

    def runall_pwr2plots(self, pklfiles, mycmap='RdBu', processlist=None, clim=None, figsave=False, figshow=True, pklsave=True, unit='dBm', noint=False):
        for whichpeak in self.thepeaks:
            if not noint:
                self.pwr2plots(pklfiles, whichpeak, mycmap,
                               processlist, clim, figsave, figshow, pklsave, unit)
            else:
                self.pwr2plots_noint(pklfiles, whichpeak, mycmap,
                                     processlist, clim, figsave, figshow, pklsave, unit)

    def pwr2plots(self, pklfiles, whichpeak, mycmap='RdBu', processlist=None, clim=None, figsave=False, figshow=True, pklsave=True, unit='dBm'):
        if unit == 'dBm':
            amplitude = 'amplitude'
        elif unit == 'Hz':
            amplitude = 'frequency'
        else:
            raise KeyError('incorrect unit: either "dBm" or "Hz"!')
        if self.ykey == 'dF (Hz)':
            mydfs = []
            for myfile in pklfiles:
                mydata = pd.read_pickle(myfile)
                yvalue = mydata[self.ykey]
                xvalue = mydata[self.xkey].iloc[0]
                theamplitude = mydata[whichpeak+' ('+unit+')']
                mydict = {}
                mydict[whichpeak + ' ('+unit+')'] = theamplitude
                mydf = pd.DataFrame(mydict)
                mydf[self.yplotkey] = yvalue/self.yscale
                mydf[self.xplotkey] = xvalue/self.xscale
                mydfs.append(mydf)
            # Here we don't need interpolation because the dFs are constant

        elif self.ykey == 'Carrier Power (dBm)':
            xx = []
            carrpow = []
            amps = []

            for myfile in pklfiles:
                mydata = pd.read_pickle(myfile)
                yvalue = mydata[self.ykey]
                xvalue = mydata[self.xkey].iloc[0]
                theamplitude = mydata[whichpeak+' ('+unit+')']
                xx.append(xvalue)
                carrpow.append(yvalue)
                amps.append(theamplitude)

            # We need interpolation because the carrier powers are slightly different from each other
            minpwr = min([min(x) for x in carrpow])
            maxpwr = max([max(x) for x in carrpow])
            newpwrs = np.linspace(minpwr, maxpwr, len(carrpow[0])*2)

            # funcs = []
            mydfs = []
            mydict = {}
            for (thepower, ii, theamplitude) in zip(carrpow, xx, amps):
                func = interp1d(thepower, theamplitude,
                                bounds_error=False)
                mydict[whichpeak + ' ('+unit+')'] = func(newpwrs)
                mydf = pd.DataFrame(mydict)
                mydf[self.yplotkey] = newpwrs/self.yscale
                mydf[self.xplotkey] = ii/self.xscale
                mydfs.append(mydf)

        try:
            print(mydfs[0].keys())
            mymtx = stlabutils.utils.stlabdict.framearr_to_mtx(
                mydfs, key=whichpeak+' ('+unit+')', xkey=self.xplotkey, ykey=self.yplotkey)
        except:
            # if script got interrupted, the last dataframe might be incomplete
            mymtx = stlabutils.utils.stlabdict.framearr_to_mtx(
                mydfs[:-1], key=whichpeak+' ('+unit+')', xkey=self.xplotkey, ykey=self.yplotkey)

        values = mymtx.pmtx.values
        maxloc = np.unravel_index(np.nanargmax(values), values.shape)
        axes = mymtx.pmtx.axes
        minIs = axes[1][maxloc[1]]
        minYY = axes[0][maxloc[0]]
        minS = values[maxloc]

        if processlist:
            mymtx.applyprocesslist(processlist)

        x = axes[1].values
        y = axes[0].values
        X, Y = np.meshgrid(x, y)
        Z = mymtx.pmtx.values
        # plt.pcolormesh(X, Y, Z, cmap=mycmap, clim=clim)
        plt.imshow(mymtx.pmtx, extent=mymtx.getextents(),
                   aspect='auto', cmap=mycmap, clim=clim)
        plt.plot(minIs, minYY, '*', c='red', ms=15)
        plt.xlabel(self.xplotkey)
        plt.ylabel(self.yplotkey)
        cbar = plt.colorbar()
        # if any('log' in x for x in processlist):
        #     cbar.set_label(
        #         r'log(Sensitivity) ($\log(\mathrm{pA}/\sqrt{\mathrm{Hz}})$)')
        #     plt.title(r'log(Sensitivity) for '+whichpeak+r', min|S| {:.3f}'.format(
        #         minS)+r' $\mathrm{pA}/\sqrt{\mathrm{Hz}}$'+r' @ {:.1f} $\mu$A, {:.1f}'.format(minIs, minYY)+' '+self.yunit)
        # else:
        cbar.set_label(
            r''+amplitude+' ('+unit+')')
        plt.title(r''+amplitude+' for '+whichpeak+r', max {:.3f}'.format(
            minS)+r' '+unit+''+r' @ {:.1f} $\mu$A, {:.1f}'.format(minIs, minYY)+' '+self.yunit)
        if figsave:
            plt.savefig(self.figdir+whichpeak+'_'+amplitude+'.png')
        if figshow:
            plt.show()
        plt.close()

        if pklsave:
            mymtx.pmtx.to_pickle(self.pkldir+whichpeak+'_'+amplitude+'.pkl')
            pickle_out = open(self.pkldir+whichpeak +
                              '_mtx_'+amplitude+'.pkl', 'wb')
            pickle.dump(mymtx, pickle_out)
            pickle_out.close()

    def pwr2plots_noint(self, pklfiles, whichpeak, mycmap='RdBu', processlist=None, clim=None, figsave=False, figshow=True, pklsave=True, unit='dBm'):
        if unit == 'dBm':
            amplitude = 'amplitude'
        elif unit == 'Hz':
            amplitude = 'frequency'
        else:
            raise KeyError('incorrect unit: either "dBm" or "Hz"!')
        # if self.ykey == 'dF (Hz)':
        mydfs = []
        for myfile in pklfiles:
            mydata = pd.read_pickle(myfile)
            yvalue = mydata[self.ykey]
            xvalue = mydata[self.xkey].iloc[0]
            theamplitude = mydata[whichpeak+' ('+unit+')']
            mydict = {}
            mydict[whichpeak + ' ('+unit+')'] = theamplitude
            mydf = pd.DataFrame(mydict)
            mydf[self.yplotkey] = yvalue/self.yscale
            mydf[self.xplotkey] = xvalue/self.xscale
            mydfs.append(mydf)
            # No interpolation as it should have been from the start. Just adding this instead of rewriting the original function to not mess with the rest of the code

        try:
            print(mydfs[0].keys())
            mymtx = stlabutils.utils.stlabdict.framearr_to_mtx(
                mydfs, key=whichpeak+' ('+unit+')', xkey=self.xplotkey, ykey=self.yplotkey)
        except:
            # if script got interrupted, the last dataframe might be incomplete
            mymtx = stlabutils.utils.stlabdict.framearr_to_mtx(
                mydfs[:-1], key=whichpeak+' ('+unit+')', xkey=self.xplotkey, ykey=self.yplotkey)

        values = mymtx.pmtx.values
        maxloc = np.unravel_index(np.nanargmax(values), values.shape)
        axes = mymtx.pmtx.axes
        minIs = axes[1][maxloc[1]]
        minYY = axes[0][maxloc[0]]
        minS = values[maxloc]

        if processlist:
            mymtx.applyprocesslist(processlist)

        x = axes[1].values
        y = axes[0].values
        X, Y = np.meshgrid(x, y)
        Z = mymtx.pmtx.values
        # plt.pcolormesh(X, Y, Z, cmap=mycmap, clim=clim)
        plt.imshow(mymtx.pmtx, extent=mymtx.getextents(),
                   aspect='auto', cmap=mycmap, clim=clim)
        plt.plot(minIs, minYY, '*', c='red', ms=15)
        plt.xlabel(self.xplotkey)
        plt.ylabel(self.yplotkey)
        cbar = plt.colorbar()
        # if any('log' in x for x in processlist):
        #     cbar.set_label(
        #         r'log(Sensitivity) ($\log(\mathrm{pA}/\sqrt{\mathrm{Hz}})$)')
        #     plt.title(r'log(Sensitivity) for '+whichpeak+r', min|S| {:.3f}'.format(
        #         minS)+r' $\mathrm{pA}/\sqrt{\mathrm{Hz}}$'+r' @ {:.1f} $\mu$A, {:.1f}'.format(minIs, minYY)+' '+self.yunit)
        # else:
        cbar.set_label(
            r''+amplitude+' ('+unit+')')
        plt.title(r''+amplitude+' for '+whichpeak+r', max {:.3f}'.format(
            minS)+r' '+unit+''+r' @ {:.1f} $\mu$A, {:.1f}'.format(minIs, minYY)+' '+self.yunit)
        if figsave:
            plt.savefig(self.figdir+whichpeak+'_'+amplitude+'.png')
        if figshow:
            plt.show()
        plt.close()

        if pklsave:
            mymtx.pmtx.to_pickle(self.pkldir+whichpeak+'_'+amplitude+'.pkl')
            pickle_out = open(self.pkldir+whichpeak +
                              '_mtx_'+amplitude+'.pkl', 'wb')
            pickle.dump(mymtx, pickle_out)
            pickle_out.close()

    def runall_noisefloor2plots(self, pklfiles, mycmap='RdBu', processlist=None, clim=None, figsave=False, figshow=True, pklsave=True):
        if self.ykey == 'dF (Hz)':
            mydfs = []
            for myfile in pklfiles:
                mydata = pd.read_pickle(myfile)
                yvalue = mydata[self.ykey]
                xvalue = mydata[self.xkey].iloc[0]
                noisefloor = mydata['Background (dBm)']
                mydict = {}
                mydict['Background (dBm)'] = noisefloor
                mydf = pd.DataFrame(mydict)
                mydf[self.yplotkey] = yvalue/self.yscale
                mydf[self.xplotkey] = xvalue/self.xscale
                mydfs.append(mydf)
            # Here we don't need interpolation because the dFs are constant

        elif self.ykey == 'Carrier Power (dBm)':
            xx = []
            carrpow = []
            noises = []

            for myfile in pklfiles:
                mydata = pd.read_pickle(myfile)
                yvalue = mydata[self.ykey]
                xvalue = mydata[self.xkey].iloc[0]
                noisefloor = mydata['Background (dBm)']
                xx.append(xvalue)
                carrpow.append(yvalue)
                noises.append(noisefloor)

            # We need interpolation because the carrier powers are slightly different from each other
            minpwr = min([min(x) for x in carrpow])
            maxpwr = max([max(x) for x in carrpow])
            newpwrs = np.linspace(minpwr, maxpwr, len(carrpow[0])*2)

            # funcs = []
            mydfs = []
            mydict = {}
            for (thepower, ii, theamplitude) in zip(carrpow, xx, noises):
                func = interp1d(thepower, theamplitude,
                                bounds_error=False)
                mydict['Background (dBm)'] = func(newpwrs)
                mydf = pd.DataFrame(mydict)
                mydf[self.yplotkey] = newpwrs/self.yscale
                mydf[self.xplotkey] = ii/self.xscale
                mydfs.append(mydf)

        try:
            print(mydfs[0].keys())
            mymtx = stlabutils.utils.stlabdict.framearr_to_mtx(
                mydfs, key='Background (dBm)', xkey=self.xplotkey, ykey=self.yplotkey)
        except:
            # if script got interrupted, the last dataframe might be incomplete
            mymtx = stlabutils.utils.stlabdict.framearr_to_mtx(
                mydfs[:-1], key='Background (dBm)', xkey=self.xplotkey, ykey=self.yplotkey)

        values = mymtx.pmtx.values
        maxloc = np.unravel_index(np.nanargmax(values), values.shape)
        axes = mymtx.pmtx.axes
        minIs = axes[1][maxloc[1]]
        minYY = axes[0][maxloc[0]]
        minS = values[maxloc]

        if processlist:
            mymtx.applyprocesslist(processlist)

        x = axes[1].values
        y = axes[0].values
        X, Y = np.meshgrid(x, y)
        Z = mymtx.pmtx.values
        # plt.pcolormesh(X, Y, Z, cmap=mycmap, clim=clim)
        plt.imshow(mymtx.pmtx, extent=mymtx.getextents(),
                   aspect='auto', cmap=mycmap, clim=clim)
        plt.plot(minIs, minYY, '*', c='red', ms=15)
        plt.xlabel(self.xplotkey)
        plt.ylabel(self.yplotkey)
        cbar = plt.colorbar()
        # if any('log' in x for x in processlist):
        #     cbar.set_label(
        #         r'log(Sensitivity) ($\log(\mathrm{pA}/\sqrt{\mathrm{Hz}})$)')
        #     plt.title(r'log(Sensitivity) for '+whichpeak+r', min|S| {:.3f}'.format(
        #         minS)+r' $\mathrm{pA}/\sqrt{\mathrm{Hz}}$'+r' @ {:.1f} $\mu$A, {:.1f}'.format(minIs, minYY)+' '+self.yunit)
        # else:
        cbar.set_label(
            r'Noise floor (dBm)')
        plt.title(r'Noise floor, max {:.3f}'.format(
            minS)+r' dBm'+r' @ {:.1f} $\mu$A, {:.1f}'.format(minIs, minYY)+' '+self.yunit)
        if figsave:
            plt.savefig(self.figdir+'_noisefloor.png')
        if figshow:
            plt.show()
        plt.close()

        if pklsave:
            mymtx.pmtx.to_pickle(self.pkldir+'_noisefloor.pkl')
            pickle_out = open(self.pkldir+'_mtx_noisefloor.pkl', 'wb')
            pickle.dump(mymtx, pickle_out)
            pickle_out.close()
