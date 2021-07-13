import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from speclite import filters
import astropy.units as u
import seaborn as sns
from scipy import interpolate
import scipy.signal

def filters_viz(z=0, phot_sys='LSST', save=False):
    """
    The function returns a plot of LSST broadband filter response
    curves with composite quasar spectrum at a given redshift.
    """
    
    # Read composite quasar spectrum (Vanden Berk et al. 2001)
    data = pd.read_csv('comp_spec.txt', skiprows=23, header=None, sep=" ", skipinitialspace=True)
    data.columns = ['Wave', 'FluxD', 'e_FluxD']
    
    # Load LSST filters
    if phot_sys=='LSST':
        filt = filters.load_filters('lsst2016-*')
        
    # Load SDSS filters
    if phot_sys=='SDSS':
        filt = filters.load_filters('sdss2010-*')
    
    
    # Plotting
    t1 = data['Wave'].values*u.Angstrom
    t1obs=(z+1)*t1
    filters.plot_filters(filt, wavelength_limits=(700, 11000), legend_loc='upper right')
    plt.plot(t1obs,(0.2*data['FluxD']/data['FluxD'].max()),label=r'$z={}$'.format(z),c='navy',linewidth=1.5)
    plt.legend(loc='upper left', fontsize=13)
    fig = plt.gcf()
    fig.set_size_inches(11,7)
    plt.xlabel(r"$\mathrm{Wavelength \ (\AA)}$",size=16, labelpad=7)
    plt.ylabel("Filter response",size=16, labelpad=7)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    
    CIV_wave = (z+1)*1549
    Ha_wave = (z+1)*6563
    Hb_wave = (z+1)*4861
    MgII_wave = (z+1)*2798
    CIII_wave = (z+1)*1908
    #Lya_wave = (z+1)*1216
    
    if CIV_wave < 11000:
        plt.annotate('C IV', xy =(CIV_wave, 0.1),
                    xytext =(CIV_wave, 0.1), size=13)
    if Ha_wave < 11000:
        plt.annotate(r'H$\alpha$', xy =(Ha_wave, 0.03),
                 xytext =(Ha_wave, 0.04), size=13)
    if Hb_wave < 11000:
        plt.annotate(r'H$\beta$', xy =(Hb_wave, 0.03),
                     xytext =(Hb_wave, 0.03), size=13)
    if MgII_wave < 11000:
        plt.annotate('Mg II', xy =(MgII_wave, 0.04),
                    xytext =(MgII_wave, 0.04), size=13)
    if CIII_wave < 11000:
        plt.annotate('C III', xy =(CIII_wave, 0.06),
                    xytext =(CIII_wave, 0.06), size=13)
    plt.grid(False)
        
    #if Lya_wave < 11000:
    #    plt.annotate(r'Ly$\alpha$ ', xy =(Lya_wave, 0.25),
    #                xytext =(Lya_wave, 0.21), size=13)
    
    if save == True:
        plt.savefig('filters.png',dpi=250)
    
    
def MJD_convert(mjd, time='days'):
    "Convert Modified Julian Date to days or hours, starting from the first data point in the given array"
    if time=='days':
        return mjd - mjd.min() # Converting MJD to days
    elif time=='hours':
        return (mjd - mjd.min())*24
    else:
        print("Error: Keyword argument 'time' is not valid, choose --> {'days', 'hours'})")
            
def mag_to_flux(mag, filt, phot_sys='SDSS'):
    # Staviti UBV, SDSS, LSST konverzije

    if filt=='u':
        return 10**((24.63 - mag)/2.5)
    elif filt =='g':
        return 10**((25.11 - mag)/2.5)
    elif filt == 'r':
        return 10**((24.80 - mag)/2.5)
    elif filt == 'i':
        return 10**((24.36 - mag)/2.5)
    elif filt == 'z':
        return 10**((22.83 - mag)/2.5)
    
def flux_to_mag(flux, filt, phot_sys='SDSS'):
    if filt == 'u':
        return 24.63 - 2.5*np.log10(flux)
    elif filt =='g':
        return 25.11 - 2.5*np.log10(flux)
    elif filt == 'r':
        return 24.80 - 2.5*np.log10(flux)
    elif filt == 'i':
        return 24.36 - 2.5*np.log10(flux)
    elif filt == 'z':
        return 22.83 - 2.5*np.log10(flux)
    
def lc_prep(t, intensity, filt, mag_flux=True, flux_mag=False, mjd_convert=True, time='days',  norm=False, phot_sys='SDSS'):
    if mag_flux == True:
        intensity = mag_to_flux(intensity, filt, phot_sys=phot_sys)
    if flux_mag == True:
        intensity = flux_to_mag(intensity, filt, phot_sys=phot_sys)
    if mjd_convert == True:
        t = MJD_convert(t, time=time)
    if norm == True:
        intensity == (intensity - intensity.min())/(intensity.max() - intensity.min())
    return t, intensity
    
    
def time_lag(acf, ccf):
    "Calculate CCF-ACF"
    # Add the left (inverted) part of the auto-correlation function (acf)
    acf_y = np.append(np.flip(acf['dcf']),acf['dcf'])
    acf_x = np.append(np.flip(acf['tau']*(-1)), acf['tau'])
    
    acf_dim = acf_y.shape[0]
    ccf_dim = ccf.shape[0]
    
    if acf_dim >= ccf_dim:
        # Interpolate the cross-correlation function data so it has the same number of points as acf
        f = interpolate.interp1d(ccf['tau'], ccf['dcf'], kind='nearest',fill_value="extrapolate")
        ccf_acfgrid = f(acf_x)
        # Calculate the differenc between cross-corelation of line & continuum and auto-correlation of the continuum.
        delta = ccf_acfgrid - acf_y
        
    elif ccf_dim > acf_dim:
        # Interpolate the auto-correlation function data so it has the same number of points as ccf
        f = interpolate.interp1d(acf_x, acf_y, kind='nearest',fill_value="extrapolate")
        acf_ccfgrid = f(ccf['tau'])
        # Calculate the differenc between cross-corelation of line & continuum and auto-correlation of the continuum.
        delta = ccf['dcf'] - acf_ccfgrid
    
    return delta

def find_tau(acf,ccf):
    "Find peaks in CCF-ACF"   

    delta = time_lag(acf,ccf)
    acf_x = np.append(np.flip(acf['tau']*(-1)), acf['tau'])
    y = delta
    delta_dim = delta.shape[0]
    acf_dim = acf_x.shape[0]
    ccf_dim = ccf['tau'].shape[0]
    if acf_dim == delta_dim:
        x = acf_x
    elif ccf_dim == delta_dim:
        x = ccf.tau
    indexes, _ = scipy.signal.find_peaks(y)
    peaks = []
    for idx in indexes:
        if x[idx]>0:
            peaks.append(x[idx])
    return peaks
    
    
######################### Plotting functions ###################################

def plot_ccf_acf(acf, ccf, lims=(-8,8), locator=1, label1='ACF', label2='CCF'):
    delta = time_lag(acf,ccf)
    acf_x = np.append(np.flip(acf['tau']*(-1)), acf['tau'])
    acf_dim = acf_x.shape[0]
    ccf_dim = ccf.shape[0]
    
    fig = plt.figure(figsize=(9,6))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=0.95, wspace=0,hspace=0)

    ax1 = fig.add_subplot(211)
    ax1.plot(acf['tau'], acf['dcf'], 'o-r', markersize=3, linewidth=0.5, label=label1)
    ax1.plot(-acf['tau'], acf['dcf'], 'o-r', markersize=3, linewidth=0.5)
    ax1.plot(ccf['tau'], ccf['dcf'], 'o-k', markersize=3, linewidth=0.5, label=label2)
    ax1.set_xlabel("Time (hours)")
    ax1.grid(which='major', axis='x', linestyle='--')
    ax1.xaxis.set_major_locator(plt.MultipleLocator(locator))
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.legend(loc='upper left', fontsize=13)
    l1, l2 = lims
    ax1.set_xlim(l1,l2)

    ax2 = fig.add_subplot(212)
    if acf_dim >= ccf_dim:
        ax2.plot(acf_x, delta, 'o-k', markersize=3, linewidth=0.5, label= label2 + ' - ' + label1)
    elif ccf_dim > acf_dim:
        ax2.plot(ccf['tau'], delta, 'o-k', markersize=3, linewidth=0.5, label= label2 + ' - ' + label1)
    ax2.set_xlabel("Time (hours)", fontsize=17)
    ax2.grid(which='major', axis='x', linestyle='--')
    ax2.xaxis.set_major_locator(plt.MultipleLocator(locator))
    ax2.legend(fontsize=13, loc='lower left')
    ax2.set_xlim(l1,l2)

    fig.text(0.04, 0.5, "Correlation (A.U.)", va='center', rotation='vertical',fontsize=18)
    plt.show()    
    
    
def plot_ccf_acf_grid(acfs, ccfs, lims, label1='ACF', label2='CCF', locator=24, save=False, tau='?'):
    fig = plt.figure(figsize=(15,5))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=0.95, wspace=0.15,hspace=0)
    j = 1
    k = 3
    for acf,ccf,lim in zip(acfs, ccfs, lims):

        delta = time_lag(acf,ccf)
        acf_x = np.append(np.flip(acf['tau']*(-1)), acf['tau'])
        acf_dim = acf_x.shape[0]
        ccf_dim = ccf.shape[0]

        ax = fig.add_subplot(2,2,j) # j=1,2
        ax.plot(acf['tau'], acf['dcf'], 'o-r', markersize=3, linewidth=0.5, label=label1)
        ax.plot(-acf['tau'], acf['dcf'], 'o-r', markersize=3, linewidth=0.5)
        ax.plot(ccf['tau'], ccf['dcf'], 'o-k', markersize=3, linewidth=0.5, label=label2)
        ax.grid(which='major', axis='x', linestyle='--')
        ax.xaxis.set_major_locator(plt.MultipleLocator(1))
        x1,x2 = lim
        ax.set_xlim(x1,x2)
        plt.setp(ax.get_xticklabels(), visible=False)
        if j==1:
            ax.legend(loc='upper left', fontsize=13)
            ax.set_title('Segment 1', fontsize=15)
            ax.xaxis.set_major_locator(plt.MultipleLocator(1))
            ax.text(0.7,0.88, r'$g - \mathrm{continuum}$', transform=ax.transAxes, size=16)
            ax.text(0.7,0.78, r'$r - \mathrm{line \ (H\alpha)}$', transform=ax.transAxes, size=16)
        if j==2:
            ax.set_title('Segment 2', fontsize=15)
            ax.xaxis.set_major_locator(plt.MultipleLocator(locator))

        ax = fig.add_subplot(2,2,k) #k=3,4
        if acf_dim >= ccf_dim:
            ax.plot(acf_x, delta, 'o-k', markersize=3, linewidth=0.5, label=label2+' - '+label1)
        elif ccf_dim > acf_dim:
            ax.plot(ccf['tau'], delta, 'o-k', markersize=3, linewidth=0.5, label=label2+' - '+label1)
        ax.set_xlabel("Time (hours)", fontsize=17)
        ax.grid(which='major', axis='x', linestyle='--')
        if k==4:
            ax.xaxis.set_major_locator(plt.MultipleLocator(locator))
        ax.set_xlim(x1,x2)
        #ax.text(0.04,0.88, 'night '+str(n), transform=ax.transAxes, size=13) 
        if k==3:
            ax.legend(loc='lower left', fontsize=13)
            ax.xaxis.set_major_locator(plt.MultipleLocator(1))
            ax.text(0.8,0.10, r'$\tau = {}$h'.format(tau), transform=ax.transAxes, size=16)
        j=j+1
        k=k+1

    fig.text(0.07, 0.5, "Correlation (A.U.)", va='center', rotation='vertical',fontsize=17)
    if save==True:
        plt.savefig('gr_novo.png', dpi=250)    
    plt.show()    

    
def plot_obs_data(lcs,filters):
    fig = plt.figure(figsize=(13, 8))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=0.95, wspace=0,hspace=0)
    fig.suptitle('Observational data of NGC 4395', size=15)
    index = [1,3,5]
    time='hours'
    for j,lc,lb in zip(index, lcs, filters):
        
        ax = fig.add_subplot(3,2,j)
        #ax.plot(lc['time_{}'.format(time)], lc['flux'], 'ok', markersize = 2)
        ax.errorbar(lc['time_{}'.format(time)], lc['norm_flux'], yerr=lc['norm_flux_err'], markersize = 3, fmt='ok', capsize=2, elinewidth=0.5)
        ax.text(0.04,0.88, lb+' band', transform=ax.transAxes, size=16)
        #ax.set_xlabel('t (hours)', fontsize = 18, labelpad=10)
        #ax.set_ylabel('Magnitude', fontsize = 18, labelpad=10)
        ax.tick_params(direction='in', pad = 5, labelsize=13)
        ax.set_xlim(0,247)
        ax.set_ylim(-0.15,1.15)
        ax.xaxis.set_major_locator(plt.MultipleLocator(40))
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
        if j != 5:
            plt.setp(ax.get_xticklabels(), visible=False)
        
        ax.axvspan(0, 9, alpha=0.4)
        ax.axvspan(70, 130, alpha=0.4)
        #ax.fill_between(x, 0, 1,facecolor='green', alpha=0.5, transform=trans)
        
        ax = fig.add_subplot(3,2,j+1)
        #ax.plot(lc['time_{}'.format(time)], lc['mag'], 'o--k', markersize = 5)
        ax.errorbar(lc['time_{}'.format(time)], lc['norm_flux'], yerr=lc['norm_flux_err'], markersize = 4, fmt='ok', capsize=2, elinewidth=1)
        #ax.set_xlabel('t (hours)', fontsize = 18, labelpad=10)
        #ax.set_ylabel('Magnitude', fontsize = 18, labelpad=10)
        ax.tick_params(direction='in', pad = 5, labelsize=13)
        ax.set_xlim(0.1,8)
        ax.set_ylim(-0.15,1.15)
        plt.setp(ax.get_yticklabels(), visible=False)
        if j+1 != 6:
            plt.setp(ax.get_xticklabels(), visible=False)
        ax.xaxis.set_major_locator(plt.MultipleLocator(1))
        ax.yaxis.set_major_locator(plt.MultipleLocator(0.2))
            
    fig.text(0.5, 0.05, 'Time (hours)', ha='center', fontsize=16)
    fig.text(0.06, 0.5, 'Normalized flux', va='center', rotation='vertical',fontsize=16)
    plt.show()

def plot_cnp_data(lcs_cnn, lcs_cnn8, lcs, filters=['r','i','g']):
    fig = plt.figure(figsize=(17, 8))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=0.95, wspace=0,hspace=0)
    fig.suptitle('Light curve model of NGC 4395 obtained using CNN', size=15)
    
    index = [1,4,7]
    time='hours'
    for j,lc_cnn, lc_cnn8, lc, lb in zip(index, lcs_cnn, lcs_cnn8, lcs, filters):
        
        ax = fig.add_subplot(3,3,j)
        ax.plot(lc['time_{}'.format(time)], lc['flux'], 'ok', markersize = 3)
        ax.plot(lc_cnn['time_{}'.format(time)], lc_cnn['flux'], color='b', alpha=0.7, linewidth=3)
               #,'o-b', markersize = 2, linewidth=2,)
        ax.fill_between(lc_cnn['time_{}'.format(time)], lc_cnn['flux']-lc_cnn['conf'], lc_cnn['flux']+lc_cnn['conf'],alpha=0.3)
        #ax.errorbar(lc['time_{}'.format(time)], lc['flux'], yerr=lc['flux_err'], markersize = 3, fmt='ok', capsize=2, elinewidth=0.5)
        ax.text(0.04,0.88, lb+' band', transform=ax.transAxes, size=17)
        #ax.set_xlabel('t (hours)', fontsize = 18, labelpad=10)
        #ax.set_ylabel('Magnitude', fontsize = 18, labelpad=10)
        ax.tick_params(direction='in', pad = 5, labelsize=13)
        ax.set_xlim(0,247)
        interval=(lc_cnn['flux'].max()-lc_cnn['flux'].min())/4
        ax.set_ylim(lc_cnn['flux'].min()-interval,lc_cnn['flux'].max()+interval)
        ax.xaxis.set_major_locator(plt.MultipleLocator(40))
        ax.yaxis.set_major_locator(plt.MaxNLocator(4))
        if j != 7:
            plt.setp(ax.get_xticklabels(), visible=False)
        
        ax = fig.add_subplot(3,3,j+1)
        #ax.plot(lc['time_{}'.format(time)], lc['flux'], 'ok', markersize = 3)
        ax.plot(lc_cnn['time_{}'.format(time)], lc_cnn['flux'], '-b', markersize = 7,linewidth=3)
        ax.fill_between(lc_cnn['time_{}'.format(time)], lc_cnn['flux']-lc_cnn['conf'], lc_cnn['flux']+lc_cnn['conf'],alpha=0.3)
        ax.errorbar(lc['time_{}'.format(time)], lc['flux'], yerr=lc['flux_err'], markersize = 4, fmt='ok', capsize=2, elinewidth=0.5)
        #ax.set_xlabel('t (hours)', fontsize = 18, labelpad=10)
        #ax.set_ylabel('Magnitude', fontsize = 18, labelpad=10)
        ax.tick_params(direction='in', pad = 5, labelsize=13)
        ax.set_xlim(0.01,8)
        ax.set_ylim(lc_cnn['flux'].min()-interval,lc_cnn['flux'].max()+interval)
        plt.setp(ax.get_yticklabels(), visible=False)
        if j+1 != 8:
            plt.setp(ax.get_xticklabels(), visible=False)
        ax.xaxis.set_major_locator(plt.MultipleLocator(1))
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        
        ax = fig.add_subplot(3,3,j+2)
        #ax.plot(lc['time_{}'.format(time)], lc['flux'], 'ok', markersize = 5)
        ax.errorbar(lc['time_{}'.format(time)], lc['flux'], yerr=lc['flux_err'], markersize = 4, fmt='ok', capsize=2, elinewidth=0.5)
        ax.plot(lc_cnn8['time_{}'.format(time)], lc_cnn8['flux'], color='b', alpha=0.7, linewidth=3)
        ax.fill_between(lc_cnn8['time_{}'.format(time)], lc_cnn8['flux']-lc_cnn8['conf'], lc_cnn8['flux']+lc_cnn8['conf'],alpha=0.3)
        ax.tick_params(direction='in', pad = 5, labelsize=13)
        ax.set_xlim(0.01,8)
        ax.set_ylim(lc_cnn['flux'].min()-interval,lc_cnn['flux'].max()+interval)
        plt.setp(ax.get_yticklabels(), visible=False)
        if j+2 != 9:
            plt.setp(ax.get_xticklabels(), visible=False)
        ax.xaxis.set_major_locator(plt.MultipleLocator(1))
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
            
    fig.text(0.5, 0.05, 'Time (hours)', ha='center', fontsize=16)
    fig.text(0.07, 0.5, 'Flux (A.U.)', va='center', rotation='vertical',fontsize=16)
    plt.show()
    
def plot_ccf_acf2(acfs, ccfs, lims=(-8,8), locator=1, labels1=['ACF','ACF'], labels2=['CCF','CCF'], tau=[2.7,2.3]):
    fig = plt.figure(figsize=(15,5))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=0.95, wspace=0.1,hspace=0)
    for i,acf,ccf,label1,label2,t in zip(range(1,3),acfs,ccfs, labels1,labels2,tau):
    
        delta = time_lag(acf,ccf)
        acf_x = np.append(np.flip(acf['tau']*(-1)), acf['tau'])
        acf_dim = acf_x.shape[0]
        ccf_dim = ccf.shape[0]
    
        ax1 = fig.add_subplot(2,2,i)
        ax1.plot(acf['tau'], acf['dcf'], 'o-r', markersize=3, linewidth=0.5, label=label1)
        ax1.plot(-acf['tau'], acf['dcf'], 'o-r', markersize=3, linewidth=0.5)
        ax1.plot(ccf['tau'], ccf['dcf'], 'o-k', markersize=3, linewidth=0.5, label=label2)
        ax1.set_xlabel("Time (hours)")
        ax1.grid(which='major', axis='x', linestyle='--')
        ax1.xaxis.set_major_locator(plt.MultipleLocator(locator))
        plt.setp(ax1.get_xticklabels(), visible=False)
        ax1.legend(loc='upper left', fontsize=13)
        l1, l2 = lims
        ax1.set_xlim(l1,l2)
        if i==1:
            ax1.text(0.7,0.88, r'$i - \mathrm{continuum}$', transform=ax1.transAxes, size=16)
            ax1.text(0.7,0.78, r'$r - \mathrm{line \ (H\alpha)}$', transform=ax1.transAxes, size=16)
        elif i==2:
            ax1.text(0.7,0.88, r'$i - \mathrm{continuum}$', transform=ax1.transAxes, size=16)
            ax1.text(0.7,0.78, r'$g - \mathrm{line \ (H\beta)}$', transform=ax1.transAxes, size=16)
        
        
        
        ax2 = fig.add_subplot(2,2,i+2)
        if acf_dim >= ccf_dim:
            ax2.plot(acf_x, delta, 'o-k', markersize=3, linewidth=0.5, label= label2 + ' - ' + label1)
        elif ccf_dim > acf_dim:
            ax2.plot(ccf['tau'], delta, 'o-k', markersize=3, linewidth=0.5, label= label2 + ' - ' + label1)
        ax2.set_xlabel("Time (hours)", fontsize=17)
        ax2.grid(which='major', axis='x', linestyle='--')
        ax2.xaxis.set_major_locator(plt.MultipleLocator(locator))
        ax2.legend(fontsize=13, loc='lower left')
        ax2.set_xlim(l1,l2)
        ax2.text(0.8,0.10, r'$\tau = {}$h'.format(str(t)), transform=ax2.transAxes, size=16)

    fig.text(0.07, 0.5, "Correlation (A.U.)", va='center', rotation='vertical',fontsize=18)
    plt.show()

