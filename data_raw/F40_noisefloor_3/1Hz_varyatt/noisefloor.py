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
SA.SetResolutionBW(1)
SA.SetAveragesType('POW')
SA.SetAverages(21)
SA.SetTraceAverageOn()
SA.SetCenter(f0)
SA.SetSpan(6.5e3)
SA.DisplayOn()
SA.SetReference('EXT')
SA.SetInputAtt(0)

idstring = 'SA_noisefloor_varyatt_1Hz'
prefix = 'F40'
myfilesa = stlabutils.newfile(prefix, idstring, autoindex=False)

myatts = [0, 10, 20]
for i, myatt in enumerate(myatts):
    SA.SetInputAtt(myatt)
    datasa = SA.MeasureScreen()
    datasa['Input Att (dB)'] = SA.GetInputAtt()

    stlabutils.saveframe(myfilesa, datasa)
    stlabutils.metagen.fromarrays(myfilesa,
                                  datasa['Frequency (Hz)'],
                                  myatts[0:i + 1],
                                  xtitle='Frequency (Hz)',
                                  ytitle='Input Att (dB)',
                                  colnames=list(datasa))