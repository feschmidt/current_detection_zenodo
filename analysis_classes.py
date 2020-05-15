# %% [markdown]
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# %% [markdown]
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# %% [markdown]
# <h1>Table of Contents<span class="tocSkip"></span></h1>
# <div class="toc"><ul class="toc-item"></ul></div>

# %%
# Libraries

from src.class_filerunner import FileRunner
from src.class_fileplotter import FilePlotter
import glob
import numpy as np

# %%
# Analysis single junction (pos2) - for curr in biascurrents: for pwr in pumppowers

myrunner = FileRunner({'pos': 'pos2', 'dev': 'singleJJ', 'RBWconst': False})
megafiles = sorted(
    glob.glob('data_raw/F32_megameasurement3/F33E*SA_vs_Is_*/*.dat'))

# %%
# Read raw data and extract peaks
myrunner.raw2peaks_current_pwr(megafiles)

# %%
# Read pickled peak data and calculate sensitivity
pklfiles = sorted(glob.glob('data_processed/pos2/F33E*.pkl'))
myrunner.peaks2sens_current(pklfiles, figsave=False, xkey='pwr')

# %%
# Make the plots
pklfiles = sorted(glob.glob('data_sensitivity/pos2/sensitivity_F33E*.pkl'))
myplotter = FilePlotter(
    {'pos': 'pos2', 'id': 'F33', 'ykey': 'Carrier Power (dBm)'})
myplotter.runall_pwr2plots(pklfiles, processlist=[
                           'flip 0,1'], figsave=False, figshow=True, mycmap='coolwarm',unit='Hz')
myplotter.runall_sens2plots(pklfiles, processlist=[
                            'nan_greater 400', 'flip 0,1'], figsave=False, figshow=True, mycmap='coolwarm_r')
myplotter.runall_pwr2plots(pklfiles, processlist=[
                           'flip 0,1'], figsave=False, figshow=True, mycmap='coolwarm',noint=True)

myplotter.runall_noisefloor2plots(pklfiles, processlist=[
                                  'flip 0,1'], figsave=False, figshow=True, mycmap='coolwarm')

# %%
# Analysis small junction (pos2) - for curr in biascurrents: for dF in freqdetuning

myrunner = FileRunner({'pos': 'pos2', 'dev': 'singleJJ'})
megafiles = sorted(
    glob.glob('data_raw/F32_megameasurement3/F33C*SA_vs_Is_*/*.dat'))

# %%
# Read raw data and extract peaks
myrunner.raw2peaks_current_dF(megafiles)

# %%
# Read pickled peak data and calculate sensitivity
pklfiles = sorted(glob.glob('data_processed/pos2/F33C*.pkl'))
myrunner.peaks2sens_current(pklfiles, figsave=False, xkey='dF')

# %%
# Load calculated data and make a single 2D plot for frequency detuning
pklfiles = sorted(glob.glob('data_sensitivity/pos2/sensitivity_F33C*.pkl'))
myplotter = FilePlotter({'pos': 'pos2', 'id': 'F33', 'ykey': 'dF (Hz)'})
myplotter.runall_pwr2plots(pklfiles, processlist=[
                           'flip 0,1'], figsave=False, figshow=True, mycmap='coolwarm',unit='Hz')
myplotter.runall_sens2plots(pklfiles, processlist=[
                            'nan_greater 3000', 'log10', 'flip 0,1'], figsave=False, figshow=True, mycmap='coolwarm_r')
myplotter.runall_pwr2plots(pklfiles, processlist=[
                           'flip 0,1'], figsave=False, figshow=True, mycmap='coolwarm')
myplotter.runall_noisefloor2plots(pklfiles, processlist=[
                                  'flip 0,1'], figsave=False, figshow=True, mycmap='coolwarm')

# %%
