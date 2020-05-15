import stlab
import stlabutils
import numpy as np
import time
from stlab.devices.IVVI import IVVI_DAC
from stlab.devices.BFWrapper import BFWrapper

SA = stlab.adi('TCPIP::192.168.1.114::INSTR')  # FSV-13

f0 = 7.5e9

# Signal analyzer
print('Setting the SA')
SA.SetPoints(2001)
SA.SetVideoBW(5)
SA.SetResolutionBW(5)
SA.SetAveragesType('POW')
SA.SetAverages(21)
SA.SetTraceAverageOn()
SA.SetCenter(f0)
SA.SetSpan(6.5e3)
SA.DisplayOn()
SA.SetReference('EXT')
SA.SetInputAtt(20)

idstring = 'SA_noisefloor_allampson'
prefix = 'F40'
myfilesa = stlabutils.newfile(prefix, idstring, autoindex=False)

datasa = SA.MeasureScreen()
datasa['Input Att (dB)'] = SA.GetInputAtt()

stlabutils.saveframe(myfilesa, datasa)

import matplotlib.pyplot as plt
import os

# save plot with data
datasa.plot(x='Frequency (Hz)', y='Spectrum (dBm)')
plt.title(os.path.split(myfilesa.name)[0])
plt.savefig(os.path.splitext(myfilesa.name)[0] + '.png')
plt.close()

# save plot with other plots
datasa.plot(x='Frequency (Hz)', y='Spectrum (dBm)')
plt.ylim(-130, -85)
plt.title(os.path.split(myfilesa.name)[0])
plt.savefig('plots/' + os.path.split(myfilesa.name)[0] + '.png')
plt.close()