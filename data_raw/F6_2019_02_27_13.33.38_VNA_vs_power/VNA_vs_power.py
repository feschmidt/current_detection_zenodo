# Power dependence of VNA as a reference for RF source

import numpy as np
import time
import stlab
from stlab.devices.BFWrapper import BFWrapper

## Temperature readout
myBF = BFWrapper(addr='172.22.33.89')
T = myBF.GetTemperature(6)

## PNA for RF
PNA = stlab.adi(addr='TCPIP::192.168.1.42::INSTR', reset=False)

# external attenuation: 45 dB + 20dB directional coupler
# external amplification: 2 Miteqs + 20dB directional coupler
rfpowers = np.arange(-30, 31, 1)
PNA.SetIFBW(100.)
PNA.SetPoints(1001)
PNA.SetRange(7.43e9, 7.45e9)
PNA.SetPower(rfpowers[0])
PNA.SetPowerOn()

last_time = time.time()

prefix = 'F'
idstring = 'VNA_vs_power'

for i, rfpow in enumerate(rfpowers):
    print('\n############### Power (dBm):', rfpow)

    PNA.SetPower(rfpow)
    data = PNA.MeasureScreen_pd()
    PNA.AutoScaleAll()
    data['Power (dBm)'] = rfpow
    try:
        T = myBF.GetTemperature(6)
    except:
        T = -1
    data['Temperature (K)'] = T
    if i == 0:
        myfile = stlab.newfile(prefix, idstring, data.keys(), autoindex=True)
    stlab.saveframe(myfile, data)
    stlab.metagen.fromarrays(
        myfile,
        data['Frequency (Hz)'],
        rfpowers[0:i + 1],
        xtitle='Frequency (Hz)',
        ytitle='Power (dBm)',
        colnames=list(data))

caption = 'Powersweep switch position 2 without current bias'
stlab.autoplot(
    myfile,
    'Frequency (Hz)',
    'Power (dBm)',
    'S21dB (dB)',
    title='Power dependence',
    caption=caption,
    cmap='RdBu')
myfile.close()