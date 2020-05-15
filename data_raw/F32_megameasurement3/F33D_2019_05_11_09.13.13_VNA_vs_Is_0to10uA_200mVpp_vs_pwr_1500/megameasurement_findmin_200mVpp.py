import stlab
import stlabutils
import numpy as np
import time
from stlab.devices.IVVI import IVVI_DAC
from stlab.devices.BFWrapper import BFWrapper

# Devices
PNA = stlab.adi('TCPIP::192.168.1.42::INSTR')  # PNA5222A
SA = stlab.adi('TCPIP::192.168.1.114::INSTR')  # FSV-13
IVVI = IVVI_DAC(verb=True)  # IVVI S4c
from stlab.devices.Rigol_DG1022 import Rigol_DG1022
LFgen = Rigol_DG1022('USB0::0x1AB1::0x0588::DG1D124204588::INSTR')
# stlab.adi not working for some reason
# Rigol DG1022 merged with IVVI internally
RFgen = stlab.adi('TCPIP::192.168.1.25::INSTR')  # SMB100A
myBF = BFWrapper(addr='172.22.33.89')
T = myBF.GetTemperature(6)
vmeas = stlab.adi('TCPIP::192.168.1.103::INSTR')

# Parameters to be looped over as a function of bias current
VNApowers = np.arange(-20, 31, 1)
deltaFs = np.linspace(1.5e6, -1.5e6, 101)

# PNA (used for searching the resonance)
print('Setting the PNA')
PNA.ClearAll()
PNA.AddTraces('S21')
f0 = 7.4379e9
# PNA.SetCenter(f0)
# PNA.SetSpan(50e6)
PNA.SetRange(f0 - 100e6, f0 + 20e6)  # do not recenter
VNApower_first = 0
PNA.SetPower(VNApower_first)
PNA.SetIFBW(1e3)
PNA.SetPoints(2001)
PNA.SetPowerOn()
_ = PNA.MeasureScreen_pd()
PNA.AutoScaleAll()

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

# Power calibration
print('Power calibration')
from _cal_power import calpower_VNA_RFgen
VNAfile = '_cal_F1_2019_02_28_09.56.05_VNA_vs_power_atten_dircoup_VNA_7to8GHz.dat'  #'..\_cal_F1_2019_02_28_09.56.05_VNA_vs_power_atten_dircoup_VNA_7to8GHz\_cal_F1_2019_02_28_09.56.05_VNA_vs_power_atten_dircoup_VNA_7to8GHz.dat'
RFgenfile = '_cal_F1_2019_02_28_10.25.28_VNA_vs_power_dircoupthrough_VNA_7to8GHz.dat'  #'..\_cal_F1_2019_02_28_10.25.28_VNA_vs_power_dircoupthrough_VNA_7to8GHz\_cal_F1_2019_02_28_10.25.28_VNA_vs_power_dircoupthrough_VNA_7to8GHz.dat'
calibration = calpower_VNA_RFgen(
    VNAfile=VNAfile, RFgenfile=RFgenfile, VNApower=VNApower_first, freqs=f0)

# RF generator
print('Setting the RFgen')
RFgenpower_first = VNApower_first + calibration['dS21dB (dB)']
RFgen.SetPower(RFgenpower_first)
RFgen.SetFrequency(f0)
RFgen.SetReference('EXT')
RFgen.SetPowerOff()

# DC bias with magnetic field
print('Setting the IVVI')
# IVVI.RampAllZero(tt=10.)
isgain = 10e-6  #A/V gain for isource
isdac = 3  #dac number for isource
# imgain = 100e-6
imdac = 1
# Bcurr = 110e-6
ismin = 0
ismax = 10.1e-6
dsteps = 0.1e-6
islist = np.arange(ismin, ismax + dsteps, dsteps)
# islist = np.arange(0, ismax + dsteps, dsteps)
# IVVI.RampVoltage(dac=isdac, mvoltage=islist[0] / isgain * 1e3)
# IVVI.RampVoltage(dac=imdac, mvoltage=Bcurr/imgain*1e3)

# LF generator
print('Setting the LFgen')
LFgen.SetOn()
time.sleep(0.5)
Vpplfgen = 200e-3
while True:
    LFgen.SetVpp(Vpplfgen)
    time.sleep(0.5)
    vpp = LFgen.GetVpp()
    if round(vpp,3) == round(Vpplfgen,3):
        break
lffreq = 1e3
while True:
    LFgen.SetFrequency(lffreq)
    time.sleep(0.5)
    freq = LFgen.GetFrequency()
    if round(freq,3) == round(lffreq,3):
        break
for _ in range(10):
    LFgen.SetOff()
    time.sleep(0.5)

# Keithley for DC measurement
print('Setting the Keithley')
vmeas.FastMeasurementSetup()
vgain = 1e3

print('Creating files')
idstring = 'DC_vs_Is_{}to{}uA'.format(
    int(islist[0] / 1e-6), int(islist[-1] / 1e-6))
prefix = 'F33X'
myfileDC = stlabutils.newfile(prefix, idstring, autoindex=False)
colnamesDC = ['Time (s)', 'T (K)', 'Iset (A)', 'Vmeas (V)', 'Rmeas (Ohm)']
last_time = time.time()

idstring = 'VNA_vs_Is_{}to{}uA_{}mVpp'.format(
    int(islist[0] / 1e-6), int(islist[-1] / 1e-6), int(Vpplfgen / 1e-3))
prefix = 'F33A'
myfilevna = stlabutils.newfile(prefix, idstring, autoindex=False)

idstring = 'SA_vs_Is_{}to{}uA_{}mVpp'.format(
    int(islist[0] / 1e-6), int(islist[-1] / 1e-6), int(Vpplfgen / 1e-3))
prefix = 'F33B'
myfilesa = stlabutils.newfile(prefix, idstring, autoindex=False)

#***********************************************************************
# Instead of fitting the resonance, we look for the lowest value in S21
# and use this value as resonance frequency.
#***********************************************************************

print('Waiting for 10s to give operator a chance to change things...')
time.sleep(10)

for i, curr in enumerate(islist):
    print('\n##### Iset (A):', curr, '#####')
    print('##### Running megameasurement part A #####\n')
    IVVI.RampAllZero()
    # ramp everything back to zero to get back to hysteretic regime
    # IVVI.RampVoltage(dac=imdac,mvoltage=Bcurr/imgain*1e3)
    IVVI.RampVoltage(dac=isdac,mvoltage=curr/isgain*1e3)

    for _ in range(10):
        LFgen.SetOff()
        time.sleep(0.5)
    RFgen.SetPowerOff()
    PNA.SetPowerOff()

    # DC measurement
    im = curr  #Do not measure current
    vm = float(vmeas.query('READ?')) / vgain
    isset = curr
    R = vm / im
    current_time = time.time()
    line = [current_time - last_time, T, isset, vm, R]
    stlabutils.writeline(myfileDC, line)
    myfileDC.write('\n')
    stlabutils.metagen.fromarrays(
        myfileDC, [0],
        islist,
        xtitle='None ()',
        ytitle='Is (A)',
        colnames=colnamesDC)

    PNA.SetPower(VNApower_first)
    PNA.SetPowerOn()
    datavna = PNA.MeasureScreen_pd()
    datavna['Is (A)'] = curr
    try:
        T = myBF.GetTemperature(6)
    except:
        T = -1
        print("ERROR READING TEMPERATURE")
    datavna['T (K)'] = T
    datavna['Power (dBm)'] = VNApower_first
    # datavna['I_mag set (A)'] = Bcurr
    #***************** FIND MIN AND RETURN RESULT ***********************
    z = datavna['S21dB (dB)']
    x = datavna['Frequency (Hz)']
    fnew = x[np.argmin(z)]
    f0 = fnew
    # PNA.SetCenter(f0)
    # PNA.SetSpan(50e6)
    PNA.AutoScaleAll()
    PNA.SetPowerOff()

    print('##### Running megameasurement part B #####\n')
    RFgen.SetFrequency(f0)
    try:
        calibration = calpower_VNA_RFgen(
            VNAfile=VNAfile,
            RFgenfile=RFgenfile,
            VNApower=VNApower_first,
            freqs=f0)
    except ValueError:
        pass  # use old calibration value. This is visible in the measurement by insane f0 or missing sidebands
    RFgen.SetPower(VNApower_first + calibration['dS21dB (dB)'])
    SA.SetCenter(f0)
    RFgen.SetPowerOn()
    for _ in range(10):
        LFgen.SetOn()
        time.sleep(0.5)
    datasa = SA.MeasureScreen()
    datasa['Is (A)'] = curr
    try:
        T = myBF.GetTemperature(6)
    except:
        T = -1
        print("ERROR READING TEMPERATURE")
    datasa['T (K)'] = T
    vpp = LFgen.GetVpp()
    datasa['Vmodpp (V)'] = vpp
    datasa['Imodpp (A)'] = vpp * isgain * 0.01  # IVVI reduction
    datasa['Modfrec (Hz)'] = LFgen.GetFrequency()
    datasa['Carrier Power (dBm)'] = RFgen.GetPower()
    # datasa['I_mag set (A)'] = Bcurr

    stlabutils.saveframe(myfilevna, datavna)
    stlabutils.metagen.fromarrays(
        myfilevna,
        datavna['Frequency (Hz)'],
        islist[:i + 1],
        xtitle='Frequency (Hz)',
        ytitle='Is (A)',
        colnames=list(datavna))
    stlabutils.saveframe(myfilesa, datasa)
    stlabutils.metagen.fromarrays(
        myfilesa,
        datasa['Frequency (Hz)'],
        islist[:i + 1],
        xtitle='Frequency (Hz)',
        ytitle='Is (A)',
        colnames=list(datasa))

    # now for the pump detuning and power sweep

    for idF, deltaF in enumerate(deltaFs):
        print('\n##### dF (Hz):', deltaF, '#####')
        print('##### Iset (A):', curr, '#####')
        print('##### Running megameasurement part C #####\n')
        RFgen.SetFrequency(f0 + deltaF)
        SA.SetCenter(f0 + deltaF)
        RFgen.SetPowerOn()
        [LFgen.SetOn() for _ in range(2)]
        datasa_dF = SA.MeasureScreen()
        datasa_dF['dF (Hz)'] = deltaF
        datasa_dF['Is (A)'] = curr
        try:
            T = myBF.GetTemperature(6)
        except:
            T = -1
            print("ERROR READING TEMPERATURE")
        datasa_dF['T (K)'] = T
        vpp = LFgen.GetVpp()
        datasa_dF['Vmodpp (V)'] = vpp
        datasa_dF['Imodpp (A)'] = vpp * isgain * 0.01  # IVVI reduction
        datasa_dF['Modfrec (Hz)'] = LFgen.GetFrequency()
        datasa_dF['Carrier Power (dBm)'] = RFgen.GetPower()
        # datasa_dF['I_mag set (A)'] = Bcurr

        if idF == 0:
            idstring_dF = 'SA_vs_Is_{}to{}uA_{}mVpp_vs_dF_{}'.format(
                int(islist[0] / 1e-6), int(islist[-1] / 1e-6),
                int(Vpplfgen / 1e-3), int(curr * 1e9))
            prefix_dF = 'F33C'
            myfilesa_dF = stlabutils.newfile(
                prefix_dF, idstring_dF, autoindex=False)

        stlabutils.saveframe(myfilesa_dF, datasa_dF)
        stlabutils.metagen.fromarrays(
            myfilesa_dF,
            datasa_dF['Frequency (Hz)'],
            deltaFs[:idF + 1],
            xtitle='Frequency (Hz)',
            ytitle='dF (Hz)',
            colnames=list(datasa_dF))

    caption = 'Sideband monitor vs pump detuning with resonance tracking via find_min from VNA signal. Frequency scale is sliding so labels are wrong'
    stlabutils.autoplot(
        myfilesa_dF,
        'Frequency (Hz)',
        'dF (Hz)',
        'Spectrum (dBm)',
        title='Spectrum vs pump detuning',
        caption=caption,
        cmap='RdBu_r')

    myfilesa_dF.close()

    for ipwr, rfpower in enumerate(VNApowers):
        IVVI.RampAllZero(tt=10.)
        # ramp everything back to zero to get back to hysteretic regime
        # IVVI.RampVoltage(dac=imdac, mvoltage=Bcurr/imgain*1e3)
        IVVI.RampVoltage(dac=isdac, mvoltage=curr/isgain*1e3)
        print('\n##### Power (dBm):', rfpower, '#####\n')
        print('##### Iset (A):', curr, '#####')
        print('##### Running megameasurement part D #####\n')
        [LFgen.SetOff() for _ in range(2)]
        RFgen.SetPowerOff()
        PNA.SetPower(rfpower)
        PNA.SetPowerOn()
        datavna_pwr = PNA.MeasureScreen_pd()
        datavna_pwr['Is (A)'] = curr
        try:
            T = myBF.GetTemperature(6)
        except:
            T = -1
            print("ERROR READING TEMPERATURE")
        datavna_pwr['T (K)'] = T
        datavna_pwr['Power (dBm)'] = rfpower
        # datavna_pwr['I_mag set (A)'] = Bcurr
        #***************** FIND MIN AND RETURN RESULT ***********************
        z = datavna['S21dB (dB)']
        x = datavna['Frequency (Hz)']
        fnew = x[np.argmin(z)]
        f0 = fnew
        # datavna_pwr['f0 (Hz)'] = f0
        # PNA.SetCenter(f0)
        # PNA.SetSpan(50e6)
        PNA.AutoScaleAll()
        PNA.SetPowerOff()

        RFgen.SetFrequency(f0)
        try:
            calibration = calpower_VNA_RFgen(
                VNAfile=VNAfile,
                RFgenfile=RFgenfile,
                VNApower=rfpower,
                freqs=f0)
        except:
            pass
        RFgen.SetPower(rfpower + calibration['dS21dB (dB)'])
        SA.SetCenter(f0)
        RFgen.SetPowerOn()
        [LFgen.SetOn() for _ in range(2)]
        datasa_pwr = SA.MeasureScreen()
        datasa_pwr['Is (A)'] = curr
        try:
            T = myBF.GetTemperature(6)
        except:
            T = -1
            print("ERROR READING TEMPERATURE")
        datasa_pwr['T (K)'] = T
        vpp = LFgen.GetVpp()
        datasa_pwr['Vmodpp (V)'] = vpp
        datasa_pwr['Imodpp (A)'] = vpp * isgain * 0.01  # IVVI reduction
        datasa_pwr['Modfrec (Hz)'] = LFgen.GetFrequency()
        datasa_pwr['Carrier Power (dBm)'] = RFgen.GetPower()
        # datasa_pwr['I_mag set (A)'] = Bcurr

        if ipwr == 0:
            idstring_vna_pwr = 'VNA_vs_Is_{}to{}uA_{}mVpp_vs_pwr_{}'.format(
                int(islist[0] / 1e-6), int(islist[-1] / 1e-6),
                int(Vpplfgen / 1e-3), int(curr * 1e9))
            prefix_pwr = 'F33D'
            myfilevna_pwr = stlabutils.newfile(
                prefix_pwr, idstring_vna_pwr, autoindex=False)

            idstring_sa_pwr = 'SA_vs_Is_{}to{}uA_{}mVpp_vs_pwr_{}'.format(
                int(islist[0] / 1e-6), int(islist[-1] / 1e-6),
                int(Vpplfgen / 1e-3), int(curr * 1e9))
            prefix_pwr = 'F33E'
            myfilesa_pwr = stlabutils.newfile(
                prefix_pwr, idstring_sa_pwr, autoindex=False)

        stlabutils.saveframe(myfilevna_pwr, datavna_pwr)
        stlabutils.metagen.fromarrays(
            myfilevna_pwr,
            datavna_pwr['Frequency (Hz)'],
            VNApowers[:ipwr + 1],
            xtitle='Frequency (Hz)',
            ytitle='Power (dBm)',
            colnames=list(datavna_pwr))
        stlabutils.saveframe(myfilesa_pwr, datasa_pwr)
        stlabutils.metagen.fromarrays(
            myfilesa_pwr,
            datasa_pwr['Frequency (Hz)'],
            VNApowers[:ipwr + 1],
            xtitle='Frequency (Hz)',
            ytitle='Carrier Power (dBm)',
            colnames=list(datasa_pwr))

    caption = 'VNA vs Power'
    stlabutils.autoplot(
        myfilevna_pwr,
        'Frequency (Hz)',
        'Power (dBm)',
        'S21dB (dB)',
        title='S21 vs VNA power',
        caption=caption,
        cmap='RdBu')
    caption = 'Sideband monitor vs Power with resonance tracking via find_min from VNA signal. Frequency scale is sliding so labels are wrong'
    stlabutils.autoplot(
        myfilesa_pwr,
        'Frequency (Hz)',
        'Carrier Power (dBm)',
        'Spectrum (dBm)',
        title='Spectrum vs carrier power',
        caption=caption,
        cmap='RdBu_r')

    myfilevna_pwr.close()
    myfilesa_pwr.close()

# IVVI.RampAllZero(tt=10.)
IVVI.close()
PNA.SetPowerOn()
PNA.close()
SA.close()
RFgen.SetPowerOff()
RFgen.close()
for _ in range(2):
    LFgen.SetDisplay('ON')
    LFgen.SetOff()
    LFgen.SetLocal()
LFgen.close()
vmeas.close()

caption = 'VNA vs Is with resonance tracking via find_min.'
stlabutils.autoplot(
    myfilevna,
    'Frequency (Hz)',
    'Is (A)',
    'S21dB (dB)',
    title='S21 vs bias current',
    caption=caption,
    cmap='RdBu')
caption = 'Sideband monitor vs Is with resonance tracking via find_min from VNA signal. Frequency scale is sliding so labels are wrong'
stlabutils.autoplot(
    myfilesa,
    'Frequency (Hz)',
    'Is (A)',
    'Spectrum (dBm)',
    title='Spectrum vs bias current',
    caption=caption,
    cmap='RdBu_r')

myfileDC.close()
myfilevna.close()
myfilesa.close()