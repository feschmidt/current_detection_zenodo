import stlabutils
import numpy as np
from scipy.interpolate import interp1d
import copy


def run_routine(datafile, VNApower, freqs):
    """Load and interpolate S21

        Reads in a datafile, looks for the block with matching VNApower and returns the interpolated S21 at a given frequency

        Parameters
        ----------
        datafile : _io.TextIOWrapper or str
            data calibration file handle (can be open or closed) or data file name for reading
        VNApower : float
            power setpoint of VNA
        freqs : float or array of floats
            One or more frequency values for interpolation

        Returns
        -------
        ynew : list
            List of interpolated data

        """

    print('Loading calfile', datafile)
    data0 = stlabutils.readdata.readdat(datafile)
    data = copy.deepcopy(data0)
    powerlist = []
    foundpower = False
    for block in data:
        thepower = block['Power (dBm)'][0]
        powerlist.append(thepower)
        if round(thepower) == round(VNApower):
            foundpower = True
            break
    if not foundpower:
        print('powerlist', powerlist)
        raise ValueError('VNApower outside the range of the cal file!')

    func = interp1d(
        x=block['Frequency (Hz)'], y=block['S21dB (dB)'], kind='cubic')
    ynew = func(np.array(freqs))
    return ynew


def calpower_VNA_RFgen(
        VNAfile='_cal_F1_2019_02_28_09.56.05_VNA_vs_power_atten_dircoup_VNA_7to8GHz\_cal_F1_2019_02_28_09.56.05_VNA_vs_power_atten_dircoup_VNA_7to8GHz.dat',
        RFgenfile='_cal_F1_2019_02_28_10.25.28_VNA_vs_power_dircoupthrough_VNA_7to8GHz\_cal_F1_2019_02_28_10.25.28_VNA_vs_power_dircoupthrough_VNA_7to8GHz.dat',
        VNApower=-10,
        freqs=np.arange(7.3e9, 7.5e9, 1e6)):
    """Calculate power differences
    
    Reads in the two datafiles, looks for the one with matching VNApower and calculates the difference in S21 (dB) between the two setups
    
    Parameters
    ----------
    VNAfile : _io.TextIOWrapper or str
        VNA calibration file handle (can be open or closed) or data file name for reading
    RFgenfile : _io.TextIOWrapper or str
        RFgen calibration file handle (can be open or closed) or data file name for reading
    VNApower : float
        power setpoint of VNA
    freqs : float or array of floats
        One or more frequency values for interpolation
    
    Returns
    -------
    params : dict
        Dict of freqs, interpolated VNA data, interpolated RFgen data, difference between datas
    
    """

    # VNA -> attenuators -> directional coupler -> DUT
    ynew_VNA = run_routine(VNAfile, VNApower, freqs)

    # RFgen -> directional coupler -> DUT
    ynew_RFgen = run_routine(RFgenfile, VNApower, freqs)

    # Calculate the power difference
    delta_VNA_RFgen = ynew_VNA - ynew_RFgen

    params = {
        'Frequency (Hz)': freqs,
        'S21dB VNA (dB)': ynew_VNA,
        'S21dB RFgen (dB)': ynew_RFgen,
        'dS21dB (dB)': delta_VNA_RFgen
    }
    return params