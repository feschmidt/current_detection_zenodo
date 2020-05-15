# %% [markdown]
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# %%

# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# %%
from src.model_currentbias import f0 as f0model
from src.model_currentbias import Lj as Ljmodel
import lmfit
import matplotlib.pyplot as plt
import glob
import stlabutils
from src.algo_peakheights import S11_fitting_full
import numpy as np
import pandas as pd
import pickle

# alternatively, this is the data with the furthest achieved tuning
myfile = glob.glob(
    'data_raw/F30_50Hz_noise_debugging/F30_2019_03_22_14.09.44_VNA_vs_current_10uAalloffcleanracks/*.dat')
myfile = myfile[0]
myfile
data = stlabutils.readdata.readdat(myfile)


# %%


# 18.8s runtime

S11res_vs_I, mydfs = S11_fitting_full(
    data[:-7], xvariable='Iset (A)', savefits=True)
# Clipping unreasonably large resonance frequencies (failed fits)
f0sraw = S11res_vs_I['f0 (Hz)']
israw = S11res_vs_I['Iset (A)']
idx = np.where(np.array(f0sraw) <= 7.44e9)[0]
S11res_vs_I_clip = S11res_vs_I.iloc[idx]
mydfs_clip = [mydfs[x] for x in idx]
# Generating the dict to be exported
fitpars = pd.DataFrame({'Is (A)': abs(S11res_vs_I_clip['Iset (A)']), 'f0 (Hz)': abs(S11res_vs_I_clip['f0 (Hz)']),
                        'Qint ()': [S11res_vs_I_clip.iloc[x]['fitpars']['Qint'].value for x in range(len(S11res_vs_I_clip))],
                        'Qext ()': [S11res_vs_I_clip.iloc[x]['fitpars']['Qext'].value for x in range(len(S11res_vs_I_clip))],
                        'theta (rad)': [S11res_vs_I_clip.iloc[x]['fitpars']['theta'].value for x in range(len(S11res_vs_I_clip))]})
f0_of_I = fitpars.set_index('Is (A)')
f0_of_I.head()

# pickle.dump(f0_of_I,
#             open("data_processed/sensitivity_estimation/fitpars_derivatives_vs_Is_F30.pkl", "wb"))


# %%


f0_of_I.tail()
# we can see that the last row is a bad fit. Checking the testplot above reveals that this row should be omitted


# %%


# f0_of_I = pickle.load(open(pklfilename, "rb"))
Iset = f0_of_I.axes[0].values

meascurr = np.array(Iset)
measf0 = f0_of_I['f0 (Hz)'].values
plt.plot(meascurr, measf0, 'o')
plt.grid()


# %%


fres = f0_of_I['f0 (Hz)'].values
kext = fres/f0_of_I['Qext ()'].values
kint = fres/f0_of_I['Qint ()'].values
ktot = kext+kint


# %%


fres[0], kext[0], kint[0], ktot[0]


# %%


testplot = True
testplot2 = True

# here we omit the bad fit
meascurr = meascurr[:-1]
measf0 = measf0[:-1]
max(meascurr)


# %%


fr0 = 7.5602e9
Lr0 = 4e-9
Ic0 = 9e-6
nJJ = 1


# %%


mymodel = lmfit.Model(f0model)
params = mymodel.make_params(fr=fr0, Lr=Lr0, Ic=Ic0, nJJ=nJJ)
params['nJJ'].vary = False
result_Icvary = mymodel.fit(measf0, params, I=meascurr)
params['Ic'].vary = False
result_Icfix = mymodel.fit(measf0, params, I=meascurr)


# %%


result_Icvary.params


# %%


popt = [result_Icvary.params['fr'].value,
        result_Icvary.params['Lr'].value, result_Icvary.params['Ic'].value]


# %%


plt.scatter(meascurr, measf0, label='data', edgecolors='k', facecolors='None')
# plt.plot(meascurr,result_Icfix.init_fit,label='initial guess')
plt.plot(meascurr, result_Icfix.best_fit, label='Ic={:.3f}uA'.format(
    result_Icfix.params['Ic'].value/1e-6))
plt.plot(meascurr, result_Icvary.best_fit, label='Ic={:.3f}uA'.format(
    result_Icvary.params['Ic'].value/1e-6))
plt.xlim(0, Ic0)
plt.legend()
mytxt = r'fr={:.3f} GHz, Lr={:.2f} nH, Ic={:.2f} uA'.format(
    popt[0]/1e9, popt[1]/1e-9, popt[2]/1e-6)
plt.text(5e-7, f0model(meascurr, *popt)[-1], mytxt)
# plt.savefig('figures/best_fit.png')
plt.show()
plt.close()


# %%
plt.scatter(meascurr, measf0, label='data', edgecolors='k', facecolors='None')
ii = np.linspace(min(meascurr), max(meascurr), 837)
ytheo = f0model(ii, fr=result_Icvary.params['fr'].value,
                Lr=result_Icvary.params['Lr'].value,
                Ic=result_Icvary.params['Ic'].value,
                nJJ=result_Icvary.params['nJJ'].value)
plt.plot(ii, ytheo)
plt.xlim(min(ii), max(ii))


# %%
fig2data = {'x1': meascurr/1e-6, 'y1': measf0/1e9, 'label': 'data',
            'xtheo': ii/1e-6, 'ytheo': ytheo/1e9, 'label': 'fit',
            'xlabel': 'Bias current ($\mu$A)', 'ylabel': 'Frequency (GHz)',
           'params':result_Icvary.params}
pickle.dump(fig2data, open('data_final/fig2_panel_f0fit.pkl', 'wb'))


# %%
print('Initial guess:', popt)
print('Ic vary:', result_Icvary.best_values)
print('Ic fix:', result_Icfix.best_values)


# %%
result_Icvary.params


# %%


pickle.dump(popt, open(
    'data_processed/sensitivity_estimation/popt_simulation_JJarrays.pkl', 'wb'))
data2export = (Iset, f0model(Iset, *popt))
pickle.dump(data2export, open(
    'data_processed/sensitivity_estimation/popt_f0_vs_Iset.pkl', 'wb'))


# %%


G1 = np.gradient(ytheo/1e6, ii/1e-6)
plt.plot(ii/1e-6, G1)
plt.xlabel('Bias current (uA)')
plt.ylabel('df0/dI (MHz/uA)')


# %%


max(abs(G1))


# %%
