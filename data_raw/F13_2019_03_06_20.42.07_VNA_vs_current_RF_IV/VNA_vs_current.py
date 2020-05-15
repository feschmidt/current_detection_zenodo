import numpy as np
import time
import stlab
from stlab.devices.IVVI import IVVI_DAC
from stlab.devices.Keithley_2100 import Keithley_2100
from stlab.devices.BFWrapper import BFWrapper

## Temperature readout
myBF = BFWrapper(addr='172.22.33.89')
try:
    T = myBF.GetTemperature(6)
except:
    T = -1

## PNA for RF
PNA = stlab.adi(addr='TCPIP::192.168.1.42::INSTR',reset=False)

# external attenuation: 45 dB + 1 directional coupler
# external amplification: 2 Miteqs + 1 directional coupler
PNA.ClearAll()
PNA.AddTraces('S21')
f0 = 7.4379e9
PNA.SetRange(f0-100e6,f0+20e6) # do not recenter
VNApower_first = -5
PNA.SetPower(VNApower_first)
PNA.SetIFBW(1e3)
PNA.SetPoints(2001)
PNA.SetPowerOn()
_ = PNA.MeasureScreen_pd()
PNA.AutoScaleAll()

## IVVI and Keithley for DC
ivvi = IVVI_DAC(verb=True)
ivvi.RampAllZero(tt=10.)

vmeas = Keithley_2100(verb=True)
vmeas.SetRangeAuto(False)
vmeas.SetRange(10)
#vmeas.write('TRIG:COUN 1')
vmeas.write('VOLT:NPLC 1')
vmeas.write('TRIG:SOUR IMM')
vmeas.write('SENS:GAIN:AUTO OFF')
vmeas.write('SENS:ZERO:AUTO OFF')

vgain = 1000 #V/V gain for vmeas
isgain = 10e-6 #A/V gain for isource (forward bias)
isdac = 3 #dac number for isource

ismax = 8e-6
ismin = -ismax
steps = 81
islist = np.linspace(ismin,ismax,steps)

last_time = time.time()

prefix = 'F'
idstring = 'VNA_vs_current'

ivvi.RampVoltage(isdac,islist[0]/isgain*1e3,tt=10.)
for i_D,curr in enumerate(islist): 
    print('\n############### Iset (A): {}'.format(curr))
    ivvi.SetVoltage(isdac,curr/isgain*1e3) #Take curr in amps, apply gain and convert to millivolts
    vm = float(vmeas.query('READ?'))/vgain
    isset = curr
    current_time = time.time()
    data = PNA.MeasureScreen_pd() #Trigger 2 port measurement and retrieve data in Re,Im format.  Returns OrderedDict
    try:
        T = myBF.GetTemperature(6)
    except:
        T = -1
    data['Power (dBm)'] = VNApower_first
    data['T (K)'] = T
    data['Time (s)'] = current_time
    data['Iset (A)'] = isset
    data['Vmeas (V)'] = vm
    
    if i_D==0:
        myfile = stlab.newfile(prefix,idstring+'_RF_IV', data.keys(), autoindex=True)
    stlab.saveframe(myfile, data)
    stlab.metagen.fromarrays(myfile,data['Frequency (Hz)'],islist[0:i_D+1],xtitle='Frequency (Hz)',ytitle='Iset (A)',colnames=list(data))

ivvi.RampAllZero(tt=10.)
PNA.close()
ivvi.close()
vmeas.close()

caption = 'VNA vs current'
stlab.autoplot(myfile,'Frequency (Hz)','Iset (A)','S21dB (dB)',title='Current dependence',caption=caption,cmap='RdBu')
myfile.close()