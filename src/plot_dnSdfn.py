# Making overview plots for d^n(S)/df^n and other quantities

import matplotlib.pyplot as plt
import stlabutils
import numpy as np
from scipy.interpolate import interp1d


def Plot_model_inputoutput(Is_meas, Signal_meas, Is_model1, Signal_model1, Is_model2=None, Signal_model2=None, corder=1, porder=1, savefig=False, showfig=True, xlabel=r'Bias current ($\mu$A)', xvariable='Is', ktot=None, note=None, scriptname=None, label1='model ki fix', label2='model ki varies'):

    plt.plot(Is_meas, Signal_meas, 'o',
             label='data', markerfacecolor='none')
    if type(Is_model2) == np.ndarray:
        plt.plot(Is_model1, Signal_model1, label=label1)
        plt.plot(Is_model2, Signal_model2, label=label2)
    else:
        plt.plot(Is_model1, Signal_model1, label='model')
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel('Amplitude (dBm)')
    if ktot:
        [plt.axvline(pf*ktot/2, c='k', ls=':', label='ktot')
         for pf in [+1, -1]]
    if scriptname:
        plt.text(.5, .05, scriptname, ha='center')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    if not note:
        plt.suptitle(
            'Peak number: {}, Correction order: {}'.format(porder, corder))
    else:
        plt.suptitle(
            'Peak number: {}, {}'.format(porder, note))
    if savefig:
        if not note:
            plt.savefig(
                'input-output formalism/model_final_plot_'+xvariable+'_corder{}_porder{}.png'.format(corder, porder))
        else:
            plt.savefig(
                'input-output formalism/model_final_plot_'+xvariable+'_{}_porder{}.png'.format(note, porder))
    if showfig:
        plt.show()
    plt.close()


def Final_plot_inputoutput(Is_dfdI, Sm1, Sp1, Sm1_dBm, meas_curr, meas_signal_p1, corder=1, porder=1, savefig=False, showfig=True, xlabel=r'Bias current ($\mu$A)', xvariable='Is', ktot=None, note=None, scriptname=None):
    # currents can stand for any changing xvariable, for example power or detuning
    _, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
    plt.sca(ax1)
    plt.plot(Is_dfdI, abs(Sm1), '.-', label='$a=-{}$'.format(porder))
    plt.plot(Is_dfdI, abs(Sp1), '.-', label='$a=+{}$'.format(porder))
    plt.ylabel('Amplitude (V)')

    plt.sca(ax2)
    plt.plot(meas_curr, meas_signal_p1, '.-', label='data')
    plt.plot(Is_dfdI, Sm1_dBm, '.-', label='input-output')
    plt.ylabel('Amplitude (dBm)')

    plt.sca(ax3)
    currs = np.linspace(max([min(meas_curr), min(Is_dfdI)]),
                        min([max(meas_curr), max(Is_dfdI)]), 101)
    f1 = interp1d(meas_curr, meas_signal_p1)
    f2 = interp1d(Is_dfdI, Sm1_dBm)
    plt.plot(currs, abs(f1(currs)-f2(currs)),
             label='absolute difference (dB)')
    plt.ylabel('Amplitude (dBm)')

    for ax in [ax1, ax2, ax3]:
        ax.legend()
        ax.set_xlabel(xlabel)
        if ktot:
            [ax.axvline(pf*ktot/1e3, c='k', ls=':', label='ktot')
             for pf in [+1, -1]]
    if scriptname:
        plt.text(.5, .05, scriptname, ha='center')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    if not note:
        plt.suptitle(
            'Peak number: {}, Correction order: {}'.format(porder, corder))
    else:
        plt.suptitle(
            'Peak number: {}, {}'.format(porder, note))
    if savefig:
        if not note:
            plt.savefig(
                'input-output formalism/final_plot_'+xvariable+'_corder{}_porder{}.png'.format(corder, porder))
        else:
            plt.savefig(
                'input-output formalism/final_plot_'+xvariable+'_{}_porder{}.png'.format(note, porder))
    if showfig:
        plt.show()
    plt.close()


def Plot_raw_inputoutput(data, zkey='Spectrum (dBm)', xkey='Frequency (Hz)', ykey='dF (Hz)', cbarlabel='Spectrum (dBm)', rotate=True, xlabel='Detuning (kHz)', ylabel='Frequency (GHz)', xscale=1e3, yscale=1e9, scriptname=None):
    mymtx = stlabutils.framearr_to_mtx(data, zkey, xkey=xkey, ykey=ykey)
    if rotate:
        mymtx.applystep('rotate_ccw')

    wbval = (0.1, 0.1)
    cmap = 'RdBu'
    lims = np.percentile(mymtx.pmtx.values, (wbval[0], 100 - wbval[1]))
    vmin = lims[0]
    vmax = lims[1]
    extents = mymtx.getextents()

    # Full overview
    x = mymtx.pmtx.axes[1]/xscale
    y = mymtx.pmtx.axes[0]/yscale
    X, Y = np.meshgrid(x, y)
    Z = mymtx.pmtx.values
    plt.pcolormesh(X, Y, Z, cmap=cmap, vmin=vmin,
                   vmax=vmax, linewidth=0, rasterized=True)
    # plt.imshow(mymtx.pmtx, aspect='auto', cmap=cmap, extent=(
    #     extents[0]/xscale, extents[1]/xscale, extents[2]/yscale, extents[3]/yscale), vmin=vmin, vmax=vmax)
    cbar = plt.colorbar()
    cbar.set_label(cbarlabel)
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title('Full overview')
    if scriptname:
        plt.text(.5, .05, scriptname, ha='center')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    plt.close()


def plot3_2Dmaps(mydfs, zkey1, zkey2, zkey3,
                 xkey='Frequency (Hz)', ykey='Is (A)', yscale=1e-6, ylabel='Bias current (uA)',
                 applysteps=['rotate_ccw'], cmap='RdBu', supertitle=None, ylim3_2D=None, colorbar=True, scriptname=None):

    mymtx1 = stlabutils.framearr_to_mtx(
        mydfs, zkey1, xkey=xkey, ykey=ykey)
    mymtx2 = stlabutils.framearr_to_mtx(
        mydfs, zkey2, xkey=xkey, ykey=ykey)
    mymtx3 = stlabutils.framearr_to_mtx(
        mydfs, zkey3, xkey=xkey, ykey=ykey)

    if applysteps:
        for mymtx in [mymtx1, mymtx2, mymtx3]:
            for thestep in applysteps:
                mymtx.applystep(thestep)
    else:
        print('No steps to apply. Continuing...')

    wbval = (0.1, 0.1)
    lims1 = np.percentile(mymtx1.pmtx.values, (wbval[0], 100 - wbval[1]))
    vmin1 = lims1[0]
    vmax1 = lims1[1]
    extents1 = mymtx1.getextents()

    lims2 = np.percentile(mymtx2.pmtx.values, (wbval[0], 100 - wbval[1]))
    vmin2 = lims2[0]
    vmax2 = lims2[1]
    extents2 = mymtx2.getextents()

    lims3 = np.percentile(mymtx3.pmtx.values, (wbval[0], 100 - wbval[1]))
    vmin3 = lims3[0]
    vmax3 = lims3[1]
    extents3 = mymtx3.getextents()

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 4))
    plt.sca(ax1)
    x = mymtx1.pmtx.axes[1]/yscale
    y = mymtx1.pmtx.axes[0]/1e9
    X, Y = np.meshgrid(x, y)
    Z = mymtx1.pmtx.values
    plt.pcolormesh(X, Y, Z, cmap=cmap, vmin=vmin1,
                   vmax=vmax1, linewidth=0, rasterized=True)
    # plt.imshow(mymtx1.pmtx, aspect='auto', cmap=cmap, extent=(
    #     extents1[0]/yscale, extents1[1]/yscale, extents1[2]/1e9, extents1[3]/1e9), vmin=vmin1, vmax=vmax1)
    # plt.xlim(0,8)
    # cbar.set_label('S11dB (dB)')

    plt.sca(ax2)
    x = mymtx2.pmtx.axes[1]/yscale
    y = mymtx2.pmtx.axes[0]/1e9
    X, Y = np.meshgrid(x, y)
    Z = mymtx2.pmtx.values
    plt.pcolormesh(X, Y, Z, cmap=cmap, vmin=vmin2,
                   vmax=vmax2, linewidth=0, rasterized=True)
    # plt.imshow(mymtx2.pmtx, aspect='auto', cmap=cmap, extent=(
    #     extents2[0]/yscale, extents2[1]/yscale, extents2[2]/1e9, extents2[3]/1e9), vmin=vmin2, vmax=vmax2)
    # plt.xlim(0,8)
    # cbar.set_label('S11dB (dB)')

    plt.sca(ax3)
    x = mymtx3.pmtx.axes[1]/yscale
    y = mymtx3.pmtx.axes[0]/1e9
    X, Y = np.meshgrid(x, y)
    Z = mymtx3.pmtx.values
    plt.pcolormesh(X, Y, Z, cmap=cmap, vmin=vmin3,
                   vmax=vmax3, linewidth=0, rasterized=True)
    # plt.imshow(mymtx3.pmtx, aspect='auto', cmap=cmap, extent=(
    #     extents3[0]/yscale, extents3[1]/yscale, extents3[2]/1e9, extents3[3]/1e9), vmin=vmin3, vmax=vmax3)
    # plt.xlim(0,8)
    # cbar.set_label('S11dB (dB)')

    for ax, zkey in zip([ax1, ax2, ax3], [zkey1, zkey2, zkey3]):
        ax.set_xlabel(ylabel)
        ax.set_ylabel('Frequency (GHz)')
        ax.set_title(zkey)
        plt.sca(ax)
        if colorbar:
            _ = plt.colorbar()

    if ylim3_2D:
        [ax.set_ylim(ylim3_2D[0], ylim3_2D[1]) for ax in [ax1, ax2, ax3]]
    if scriptname:
        fig.text(.5, .05, scriptname, ha='center')
    if supertitle:
        plt.suptitle(supertitle)
    else:
        print('No supertitle provided. Continuing...')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.show()
    plt.close()


def plot_peakamp_sens_percent(peakno, meas_curr, meas_signal_px, meas_signal_mx, theo_currs, theo_spectrum,
                              meas_sensitivity_px, meas_sensitivity_mx, theo_sensitivity,
                              supertitle='Comparing measurement and theory', savefig=True,
                              ylimsens=(0, 1e3), scriptname=None):

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 4))

    plt.sca(ax1)
    plt.plot(meas_curr/1e-6, meas_signal_px, '.',
             label='measurement +{}f1'.format(peakno))
    plt.plot(meas_curr/1e-6, meas_signal_mx, '.',
             label='measurement -{}f1'.format(peakno))
    plt.plot(theo_currs/1e-6, theo_spectrum, '.', label='theory')
    plt.legend()
    plt.ylim(-90, -48)
    plt.xlabel('Bias current (uA)')
    plt.ylabel('Peak height (dBm)')

    plt.sca(ax2)
    plt.plot(meas_curr/1e-6, meas_sensitivity_px /
             1e-12, '.', label='measurement +{}f1'.format(peakno))
    plt.plot(meas_curr/1e-6, meas_sensitivity_mx /
             1e-12, '.', label='measurement -{}f1'.format(peakno))
    plt.plot(theo_currs/1e-6, theo_sensitivity/1e-12, '.', label='theory')
    plt.legend()
    plt.ylim(ylimsens[0], ylimsens[1])
    plt.xlabel('Bias current (uA)')
    plt.ylabel('Sensitivity (pA/rtHz)')

    plt.sca(ax3)
    meas_sensitivity_int_px = interp1d(meas_curr, meas_sensitivity_px)
    meas_sensitivity_int_mx = interp1d(meas_curr, meas_sensitivity_mx)
    theo_sensitivity_int = interp1d(theo_currs, theo_sensitivity)
    full_currs = np.linspace(max([min(meas_curr), min(theo_currs)]), min(
        [max(meas_curr), max(theo_currs)]), 201)
    plt.plot(full_currs/1e-6, abs(theo_sensitivity_int(full_currs)-meas_sensitivity_int_px(full_currs)
                                  )/theo_sensitivity_int(full_currs)*100, label='(meas-theo)/theo +{}f1'.format(peakno))
    plt.plot(full_currs/1e-6, abs(theo_sensitivity_int(full_currs)-meas_sensitivity_int_mx(full_currs)
                                  )/theo_sensitivity_int(full_currs)*100, label='(meas-theo)/theo -{}f1'.format(peakno))
    plt.legend()
    plt.ylim(0, 50)
    plt.xlabel('Bias current (uA)')
    plt.ylabel('Relative sensitivity difference (%)')

    if scriptname:
        fig.text(.5, .05, scriptname, ha='center')
    if supertitle:
        plt.suptitle('Comparing measurement and theory')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    if savefig:
        plt.savefig('plots/sensitivity_estimation/sensitivity_compare_pm{}.png'.format(peakno),
                    bbox_to_inches='tight')
    plt.show()
    plt.close()


def plot_peakamp_sens_abs(peakno, meas_curr, meas_signal_px, meas_signal_mx, theo_currs, theo_spectrum,
                          meas_sensitivity_px, meas_sensitivity_mx, theo_sensitivity,
                          supertitle='Comparing measurement and theory', savefig=True,
                          ylimsens=(0, 1e3), calcmethod='Taylor-analytical', scriptname=None):

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 4))

    plt.sca(ax1)
    plt.plot(meas_curr/1e-6, meas_signal_px, '.',
             label='measurement +{}f1'.format(peakno))
    plt.plot(meas_curr/1e-6, meas_signal_mx, '.',
             label='measurement -{}f1'.format(peakno))
    plt.plot(theo_currs/1e-6, theo_spectrum, '.', label='theory')
    plt.legend()
    plt.ylim(-90, -48)
    plt.xlabel('Bias current (uA)')
    plt.ylabel('Peak height (dBm)')

    plt.sca(ax2)
    plt.plot(meas_curr/1e-6, meas_sensitivity_px /
             1e-12, '.', label='measurement +{}f1'.format(peakno))
    plt.plot(meas_curr/1e-6, meas_sensitivity_mx /
             1e-12, '.', label='measurement -{}f1'.format(peakno))
    plt.plot(theo_currs/1e-6, theo_sensitivity/1e-12, '.', label='theory')
    plt.legend()
    plt.ylim(ylimsens[0], ylimsens[1])
    plt.xlabel('Bias current (uA)')
    plt.ylabel('Sensitivity (pA/rtHz)')

    plt.sca(ax3)
    meas_sensitivity_int_px = interp1d(meas_curr, meas_signal_px)
    meas_sensitivity_int_mx = interp1d(meas_curr, meas_signal_mx)
    theo_sensitivity_int = interp1d(theo_currs, theo_spectrum)
    full_currs = np.linspace(max([min(meas_curr), min(theo_currs)]), min(
        [max(meas_curr), max(theo_currs)]), 201)
    plt.plot(full_currs/1e-6, abs(theo_sensitivity_int(full_currs)-meas_sensitivity_int_px(full_currs)
                                  ), label='abs(meas-theo) +{}f1'.format(peakno))
    plt.plot(full_currs/1e-6, abs(theo_sensitivity_int(full_currs)-meas_sensitivity_int_mx(full_currs)
                                  ), label='abs(meas-theo) -{}f1'.format(peakno))
    plt.legend()
    plt.xlabel('Bias current (uA)')
    plt.ylabel('Absolute peak difference (dB)')

    if scriptname:
        fig.text(.5, .05, scriptname, ha='center')
    if supertitle:
        plt.suptitle(supertitle+', '+calcmethod)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    if savefig:
        plt.savefig('plots/sensitivity_estimation/sensitivity_compare_'+calcmethod+'_pm{}.png'.format(peakno),
                    bbox_to_inches='tight')
    plt.show()
    plt.close()


def plot_four_overview(data, key='dSdf', unit='(1/Hz)', savefig=False, showfig=True, ii=0, scriptname=None):

    fig, (ax1, ax2) = plt.subplots(2, 2, figsize=(12, 8))
    plt.sca(ax1[0])
    plt.plot(data['Frequency (Hz)'], data
             [key+'abs '+unit], label='data abs')
    plt.plot(data['Frequency (Hz)'], data
             [key+'fitabs '+unit], label='fit abs')
    plt.legend()

    plt.sca(ax1[1])
    plt.plot(data['Frequency (Hz)'], data
             [key+'fitnobgabs '+unit], label='nobg abs')
    plt.legend()

    plt.sca(ax2[0])
    plt.plot(data['Frequency (Hz)'], data
             [key+'fitre '+unit], label='re')
    plt.plot(data['Frequency (Hz)'], data
             [key+'fitim '+unit], label='im')
    plt.legend()

    plt.sca(ax2[1])
    plt.plot(data['Frequency (Hz)'], data
             [key+'fitnobgre '+unit], label='nobg re')
    plt.plot(data['Frequency (Hz)'], data
             [key+'fitnobgim '+unit], label='nobg im')

    plt.legend()
    plt.suptitle(
        key+' '+unit+', I0={}uA'.format(abs(data['Is (A)'][0]/1e-6)))
    if scriptname:
        fig.text(.5, .05, scriptname, ha='center')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    if savefig:
        plt.savefig('plots/sensitivity_estimation/'+key+'_plots/' +
                    key+'_plots_{}.png'.format(ii+1), bbox_to_inches='tight')
    if showfig:
        plt.show()
    plt.close()


def plotall_four_overview(data, key='dSdf', unit='(1/Hz)'):

    doit = input('Do you want to plot and save all of these figures? y/n\n')

    if doit == 'y':
        print('Doing it...')
        for ii, myblock in enumerate(data):
            plot_four_overview(myblock, key, unit,
                               savefig=True, showfig=False, ii=ii)

    else:
        print('Continuing without plotting and saving...')
