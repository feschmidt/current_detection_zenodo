
# coding: utf-8

# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# We here show that both the internal and external loss rates of our device can be described by the following function:
#
# $\begin{align}
# \kappa(I) = \kappa(0) + \alpha \exp(I/I^*)
# \end{align}$
#
# Based on `model_inputoutput_vs_Iset.py`

# In[1]:


from src.model_currentbias import f0 as f0model
import numpy as np
import matplotlib.pyplot as plt
import pickle
import lmfit
from src.model_currentbias import lossrate


# In[8]:


# Cavity parameters assuming varying values (further down indicated with "_fit")
fitpars = pickle.load(
    open("data_processed/sensitivity_estimation/fitpars_vs_Is.pkl", "rb"))
Iset = fitpars.axes[0].values
Qint = fitpars['Qint ()'].values
Qext = fitpars['Qext ()'].values
dQint = fitpars['dQint ()'].values
dQext = fitpars['dQext ()'].values

plt.errorbar(Iset, Qint, yerr=dQint, fmt='o',
             markerfacecolor='none', label='Qint')
plt.errorbar(Iset, Qext, yerr=dQext, fmt='o',
             markerfacecolor='none', label='Qext')
plt.legend()
plt.grid()
plt.title('Qint drops for higher bias currents')
# plt.savefig('input-output formalism/model_Qint.png')
plt.show()
plt.close()


# In[19]:


f0 = fitpars['f0 (Hz)'].values
df0 = fitpars['df0 (Hz)'].values
kappaint = f0/Qint
kappaext = f0/Qext
dkappaint = (df0/f0 + dQint/Qint)*kappaint
dkappaext = (df0/f0 + dQext/Qext)*kappaext

plt.errorbar(Iset, kappaint/1e3, yerr=dkappaint/1e3,
             fmt='o', markerfacecolor='none', label='kint (MHz)')
plt.errorbar(Iset, kappaext/1e3, yerr=dkappaext/1e3,
             fmt='o', markerfacecolor='none', label='kext (MHz)')
plt.legend()
plt.grid()
plt.title('kint grows for higher bias currents')
# plt.savefig('input-output formalism/model_kint.png')
plt.show()
plt.close()


# In[25]:


mymodel = lmfit.Model(lossrate)
kiparams = mymodel.make_params(k0=kappaint[0], alpha=1, ix=1)
kiresult = mymodel.fit(kappaint/1e3, kiparams, x=Iset/1e-6)
keparams = mymodel.make_params(k0=kappaext[0], alpha=1, ix=1)
keresult = mymodel.fit(kappaext/1e3, keparams, x=Iset/1e-6)


# In[26]:


kiresult.params


# In[27]:


keresult.params


# In[28]:


plt.errorbar(Iset/1e-6, kappaint/1e3, yerr=dkappaint/1e3,
             fmt='o', markerfacecolor='none', label='kint (MHz)')
plt.plot(Iset/1e-6, kiresult.best_fit, 'k')  # ,label='exp fit')
plt.errorbar(Iset/1e-6, kappaext/1e3, yerr=dkappaext/1e3,
             fmt='o', markerfacecolor='none', label='kext (MHz)')
plt.plot(Iset/1e-6, keresult.best_fit, 'k')  # ,label='exp fit')
plt.ylabel('Loss rate (kHz)')
plt.xlabel('Bias current (uA)')
plt.legend()
plt.show()
plt.close()


# In[29]:


mydict = {'xmeas': Iset/1e-6, 'ki': kappaint/1e3, 'kitheo': kiresult.best_fit,
          'ke': kappaext/1e3, 'ketheo': keresult.best_fit,
          'dke': dkappaext/1e3, 'dki': dkappaint/1e3,
          'ylabel': 'Loss rate (kHz)', 'xlabel': 'Bias current (uA)'}
pickle.dump(mydict, open('data_final/SM_lossrate_expfit.pkl', 'wb'))


# In[37]:


# export the loss rate parameters
lossratedict = {'ke0 (Hz)': keresult.params['k0'].value*1e3,
                'kealpha (Hz)': keresult.params['alpha'].value*1e3,
                'keix (A)': keresult.params['ix'].value*1e-6,
                'ki0 (Hz)': kiresult.params['k0'].value*1e3,
                'kialpha (Hz)': kiresult.params['alpha'].value*1e3,
                'kiix (A)': kiresult.params['ix'].value*1e-6
                }
lossratedict


# In[38]:


ii = np.linspace(min(Iset), max(Iset), 401)

plt.plot(Iset, kappaint,
         'C0o', markerfacecolor='none', label='kint (Hz)')
plt.plot(ii, lossrate(ii, lossratedict['ki0 (Hz)'], lossratedict['kialpha (Hz)'], lossratedict['kiix (A)']),
         'C0', label='exp fit')
plt.plot(Iset, kappaext,
         'C1o', markerfacecolor='none', label='kext (Hz)')
plt.plot(ii, lossrate(ii, lossratedict['ke0 (Hz)'], lossratedict['kealpha (Hz)'], lossratedict['keix (A)']),
         'C1', label='exp fit')
plt.legend()
plt.show()
plt.close()


# In[39]:


pickle.dump(lossratedict, open(
    'data_processed/sensitivity_estimation/lossrate_expfit.pkl', 'wb'))


# In[1]:


# checking whether the internal loss rate is simply proportional to dw/dI


# In[17]:


f0fitdict = pickle.load(
    open('data_processed/sensitivity_estimation/popt_simulation_JJarrays.pkl', 'rb'))
f0fitdict


# In[14]:


# In[32]:


f0fit = f0model(ii, *f0fitdict)
G1fit = np.gradient(f0fit, ii)
plt.plot(ii, -G1fit/1e7, label='G1')
plt.plot(ii, lossrate(
    ii, lossratedict['ki0 (Hz)'], lossratedict['kialpha (Hz)'], lossratedict['kiix (A)']), label='ki fit')
plt.legend()
plt.show()
plt.close()


# In[35]:


plt.plot(ii, lossrate(
    ii, lossratedict['ki0 (Hz)'], lossratedict['kialpha (Hz)'], lossratedict['kiix (A)'])/(-G1fit))
plt.yscale('log')
plt.xlabel('Bias current (A)')
plt.ylabel('ki/(dwdI)')
plt.show()
plt.close()


# In[36]:


# This is indeed not the case, so it cannot be only current noise
