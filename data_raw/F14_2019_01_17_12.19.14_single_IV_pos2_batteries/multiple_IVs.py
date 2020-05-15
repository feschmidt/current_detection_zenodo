import stlab
from stlab.devices.IVVI import IVVI_DAC
from stlab.devices.Keithley_2100 import Keithley_2100
import numpy as np
import time
from stlab.devices.BFWrapper import BFWrapper

myBF = BFWrapper(addr='172.22.33.89')

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
isgain = 100e-6 #A/V gain for isource (forward bias)
isdac = 3 #dac number for isource

ismax = 20e-6
ismin = -ismax
steps = 401
islist = np.linspace(ismin,ismax,steps)

prefix = 'F'
idstring = 'single_IV_pos4_batteries'

colnames = ['Time (s)', 'T (K)', 'Iset (A)', 'Vmeas (V)', 'Rmeas (Ohm)']

last_time = time.time()
try:
    T = myBF.GetTemperature(6)
except:
    T = -1

myfile = stlab.newfile(prefix,idstring,colnames,autoindex=True)

for counter in range(1000):
    ivvi.RampVoltage(isdac,islist[0]/isgain*1e3,tt=10.)
    for curr in islist:
        ivvi.SetVoltage(isdac,curr/isgain*1e3) #Take curr in amps, apply gain and convert to millivolts
        im = curr #Do not measure current
        vm = float(vmeas.query('READ?'))/vgain
        isset = curr
        R = vm/im
        current_time = time.time()
        line = [current_time - last_time, T, isset, vm, R]
        stlab.writeline(myfile,line)

    myfile.write('\n')
    stlab.metagen.fromarrays(myfile,islist,range(counter+1),xtitle='Is (A)', ytitle = 'Index ()', colnames = colnames)

ivvi.RampAllZero(tt=30.)
ivvi.close()

vmeas.close()

myfile.close()
