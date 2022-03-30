"""
Authors: Isidora Jankov & dr Andjelka Kovačević
-----------------------------------------------

Description:
------------
A set of functions for:

- generating artificial AGN light curves using a Damped Random Walk process (see Kovačević et al. 2021).
- generating artificial AGN photometric light curves for photometric reverberation mapping (see Jankov et al. 2022)
- utility functions for preparing these light curves for ZDCF and PLIKE programs (Alexander 1997)

Attribution
-----------
If you use these functions for scientific work leading to a publication,
please cite:

- Kovačević et al. (2021, MNRAS, 505, 5012-5028)
- Jankov et al. (2022, Astronomische Nachrichten, 343, e210090)

where the method and code for generating artificial light curves is described, 
tested, and applied.
"""

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from scipy import interpolate
from speclite import filters
import astropy.units as u
from add_asym import add_asym

def LC_conti(T, deltatc=1, oscillations=True, A=0.14, noise=0.00005, z=0, method='default', rblr_flag=False, lag=40):
    """ 
    Function originally written by dr Andjelka Kovačević, adapted and improved by Isidora Jankov.
    
    Generate one artificial light curve using a stochastic model based on the Damped random walk (DRW)
    proccess. Parameters describing the model are characteristic amplitude ("logsig2" in code) and time 
    scale of the exponentially-decaying variability ("tau" in code), both infered from physical 
    quantities such are supermassive black hole mass and/or luminosity of the AGN. For further details
    regarding the model see Kovačević et al. (2021) and references therein.
 
    
    Parameters:
    -----------
    T: int
        Total time span of the light curve. It is recommended to generate light curves to be at least 
        10 times longer than their characteristic timescale (Kozłowski 2017). 
    deltatc: int, default=1
        Cadence (or sampling rate) - time interval between two consecutive samplings of the light 
        curve in days.
    oscillations: bool, default=True
        If True, light curve simulation will take an oscillatory signal into account.
    A: float, default=0.14
        Amplitude of the oscillatory signal in magnitudes (used only if oscillations=True).
    noise: float, default=0.00005
        Amount of noise to include in the light curve simulation.
    z: float, default=0
        Redshift.
    method: {'default', 'Kelly'}, default='default'
        Method for calculating DRW model parameters.
    rblr_flag: bool, default=False
        If True, the time-lag is user defined and keyword parameter 'lag' indicates the value. 
        Otherwise, the time-lag is randomly generated.
    lag: float, default=40
        The value of the time-lag. Considered only if the rblr_flag is True.
    
    Returns:
    --------
    tt: np.array
        Days when the light curve was sampled.
    yy: np.array
        Magnitudes of the simulated light curve.
    err: np.array
        Light curve error (Ivezić et al. 2019)
    rblr: float
        BLR radius of the simulated light curve in days (i.e., time-lag).
        
    References:
    -----------
    Ivezić, Ž., et al. 2019, ApJ, 873, 111 (https://iopscience.iop.org/article/10.3847/1538-4357/ab042c)
    Kelly, B.C., Bechtold, J., & Siemiginowska, A. 2009, ApJ, 698, 895 (https://iopscience.iop.org/article/10.1088/0004-637X/698/1/895)
    Kovačević, A., et al. 2021, submitted to MNRAS (https://github.com/LSST-sersag/white_paper/blob/main/data/paper.pdf)
    Kozłowski, S. 2017, A&A, 597, A128 (https://www.aanda.org/articles/aa/full_html/2017/01/aa29890-16/aa29890-16.html)
    """
    
    # Constants
    const1 = 0.455*1.25*1e38
    const2 = np.sqrt(1e09)
    meanmag = 23.
    
    # Generating survey days 
    tt = np.arange(1, T, int(deltatc))
    times = tt.shape[0]
    
    # Generating L_bol
    if rblr_flag == True:
        lumbol = ((lag)/(10**(1.527)))**(1/0.533)*1e44
    else:
        loglumbol = np.random.uniform(42.2,45.5,1)
        lumbol = np.power(10,loglumbol)
    print("log (L_bol) = {:.2f}".format(np.log10(float(lumbol))))

    # Calculate M_{SMBH}
    msmbh=np.power((lumbol*const2/const1),2/3.)
    print("M_bh = {:2e}".format(float(msmbh)))
    
    # Calculate DRW model parameters: 
    # damping time scale (tau) & amplitude of correlation decay (sig)
    if method == 'default':
        tau = 80.4*np.power(lumbol/1e45,-0.42)*np.power(msmbh/1e08,1.03)
        tau = tau*(1+z)
        logsig2 = -3.83-0.09*np.log(lumbol/1e45)-0.25*np.log(msmbh/1e08)
        sig = np.sqrt(np.power(10,logsig2))/np.sqrt(1+z)
        
    elif method == 'Kelly': 
        logtau = -8.13+0.24*np.log10(lumbol)+0.34*np.log10(1+z) # Eq 22, Kelly et al. 2009
        tau = np.power(10,logtau)*(1+z)
        logsig2 = 8-0.27*np.log10(lumbol)+0.47*np.log10(1+z) # Eq 17, Kelly et al. 2009
        sig = np.sqrt(np.power(10,logsig2))/np.sqrt(1+z)
    
    # Calculate the broad line region radius
    logrblr = 1.527 + 0.533*np.log10(lumbol/1e44)
    rblr = np.power(10,logrblr)
    rblr=rblr.item()
    print("Time-lag = {:.2f} days".format(rblr))
    
    # Calculating light curve magnitudes
    ss = np.zeros(times)
    ss[0] = meanmag # light curve is initialized
    SFCONST2=sig*sig
    ratio = -deltatc/tau

    for i in range(1, times):
        ss[i] = np.random.normal(ss[i-1]*np.exp(ratio) + meanmag*(1-np.exp(ratio)),
                                     np.sqrt(10*0.5*tau*SFCONST2*((1-np.exp(2*ratio)))),1)
        
    # Calculating light curve error (Ivezic et al. 2019)    
    err = lc_err(ss)
    
    # Final light curve with oscillations
    if oscillations == True:
        # Calculate underlying periodicity
        conver=173.145 # convert from LightDays to AU
        lightdays=10
        P = np.sqrt(((lightdays*conver)**3)/(msmbh))
        # Calculating and adding oscillatory signal
        sinus=A*np.sin(2*np.pi*tt/(P*365))
        ss = ss + sinus
        yy = np.zeros(times)
        for i in range(times):
            # Adding noise to each magnitude value
            yy[i] = ss[i] + np.random.normal(0,((noise*ss[i])),1)
    
        return tt, yy, err, rblr
    
    # Final light curve without oscillations
    elif oscillations == False:
        yy = np.zeros(times)
        for i in range(0,times):
            # Adding noise to each magnitude value
            yy[i] = ss[i] + np.random.normal(0,((noise*ss[i])),1)
    
        return tt, yy, err, rblr



def LC_two_bands(tc, fxc, errxc, rblr, wl=0.2, wc=0.8):
    """
    Function originally written by dr Andjelka Kovačević, adapted by Isidora Jankov.
    
    Using the continuum light curve generated by LC_conti(), simulate a light curve in a hypothetical
    band covering emission line together with the provided continuum. The method is described in Section 2
    of Chelouche & Daniel (2012).
    
    Parameters:
    -----------
    tc: np.array
    	Time dimension of the continuum light curve.
    fxc: np.array
    	Flux dimension of the continuum light curve.
    errxc: np.array
        Error in flux of the continuum light curve.
    wl: float, default=0.2
        Amount of contribution of the emission line in the total integrated flux.
    wc: float, default=0.8
        Amount of contribution of the continuum in the total integrated flux.
    
    Returns:
    --------
    x_band: pd.DataFrame
        Time, flux and error of the continuum light curve packed into a DataFrame.
    y_band: pd.DataFrame
        Time, flux and error of the continuum + line light curve packed into a DataFrame.
    line_response: np.array
        Time, flux and error of the emission line response function packed into a DataFrame.
    """

    fxc_norm = fxc/fxc.max() # normalize flux
    
    # 1. Y band continuum
    # --> Convolution between continuum flux (X band) and some transfer function (Chelouche & Daniel 2012)
    # --> We approximate this Y continuum to be the same as the X band continuum
    
    fy_c = fxc_norm
    tcem = np.arange(len(fy_c))
    
    
    # 2. Y band emission line (response)
    # --> Convolution between continuum flux (X band) and a Gaussian transfer function
    
    mu = rblr
    std = np.sqrt(mu)/2
    
    transfer_l = lambda t: 23*np.exp(-np.power(t - mu, 2.) / (np.power(std, 2.)))
    
    # Plot the Gaussian kernel
    plt.plot(transfer_l(np.arange(150)), label=r'$\psi \ (Y_l)$')
    plt.legend(fontsize=15)
    plt.title("Gaussian kernel")
    plt.xlabel("Time (days)")
    
    y_norm = np.convolve(np.ones_like(fxc_norm), transfer_l(tc), mode='full')
    valid_indices = (y_norm > 0.)
    y_norm = y_norm[valid_indices]
    response = np.convolve(fxc_norm, transfer_l(tc), 'full')[valid_indices]/y_norm
    
    tem = np.arange(len(response))
    
    
    # 3. Final Y band light curve
    # --> Light curve covering emission line and surrounding continuum
    # --> Weighted sum between Y band continuum and Y band emission line
    
    merged = wc*fxc_norm[:len(tc)] + wl*response[:len(tc)]
    
    # 4. Calculate Y band light curve error
    
    err_merg_e = lc_err(merged)
    errc_e = lc_err(fxc_norm)
    err_em_e = lc_err(wl*response)
    
    err_merg = np.sqrt(np.max(err_merg_e/merged)*merged)
    errc = np.sqrt(np.max(errc_e/fxc_norm)*fxc_norm)
    err_em = np.sqrt(np.max(err_em_e/wl*response)*wl*response)
           
    # Create DataFrames for continuum LC (X band), continuum + line LC (Y band)
    # and for emission line response
    
    data_cont = {'t':tc,
            'flux':fxc_norm[:len(tc)],
            'err':errc[:len(tc)]}
    
    data_response = {'t':tc,
            'flux':wl*response[:len(tc)],
            'err':err_em[:len(tc)]}
    
    data_merged = {'t':tc,
            'flux':merged[:len(tc)],
            'err':err_merg[:len(tc)]}
    
    x_band = pd.DataFrame(data_cont)
    y_band = pd.DataFrame(data_merged)
    line_response = pd.DataFrame(data_response)
        
    return x_band, y_band, line_response


def lc_err(lc):
    # Calculating light curve error (Ivezic et al. 2019)
    gamma=0.039
    m5=24.7
    x=np.zeros(lc.shape)
    x=np.power(10, 0.4*(lc-m5))

    err = (0.005*0.005) + (0.04-gamma)*x + gamma*x*x
    return err

def filters_viz(z=0, phot_sys='LSST', save=False):
    """
    Function originally written by dr Andjelka Kovačević, adapted by Isidora Jankov to include ZTF filters.
    
    The function returns a plot of LSST/SDSS/ZTF broadband filter response
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
    
    # Load ZTF filters
    if phot_sys=='ZTF':
        directory_name = './ZTF_data/ZTF_filters'
        fg_name = os.path.join(directory_name, 'ztf-g.ecsv')
        fr_name = os.path.join(directory_name, 'ztf-r.ecsv')
        fi_name = os.path.join(directory_name, 'ztf-i.ecsv')
        filt = filters.load_filters(fg_name,fr_name,fi_name)
    
    
    # Plotting
    t1 = data['Wave'].values*u.Angstrom
    t1obs=(z+1)*t1
    filters.plot_filters(filt, wavelength_limits=(700, 11000), legend_loc='upper right')
    plt.plot(t1obs,(0.2*data['FluxD']/data['FluxD'].max()),label=r'$z={}$'.format(z),c='navy',linewidth=1.5)
    plt.legend(loc='upper left', fontsize=14)
    fig = plt.gcf()
    fig.set_size_inches(11,7)
    plt.xlabel(r"$\mathrm{Wavelength \ (\AA)}$",size=17, labelpad=7)
    plt.ylabel("Filter response",size=18, labelpad=9)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    
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
        plt.savefig('filters.pdf',dpi=10)
    
    
def LC_plot(tt, yy, T):
    """
    Simple plotting function.
    
    Parameters:
    -----------
    tt: np.array
        Days when the light curve was sampled.
    yy: np.array
        Light curve magnitudes.
    T: int
        Total time span of the light curve. 
    """
    
    fig = plt.figure(figsize=(15,5))
    
    ax = fig.add_subplot(111)
    ax.plot(tt, yy, 'ko', markersize = 1, label='1 day cadence')

    custom_xlim = (0, T)
    custom_ylim = (yy.min()-0.1, yy.max()+0.1)
    ax.set_xlabel('t [days]', fontsize = 18, labelpad=10)
    ax.set_ylabel('magnitude', fontsize = 18, labelpad=10)
    ax.tick_params(direction='in', pad = 5, labelsize=13)
    plt.setp(ax, xlim=custom_xlim, ylim=custom_ylim)
    ax.legend(fontsize=15)
    ax.grid(True)


#---------------------------------------------------------------------------------------#
# Functions for preparaing light curves for ZDCF and PLIKE programs (by Isidora Jankov) #
#---------------------------------------------------------------------------------------#

def add_inverted_acf(acf):
    new_acf = pd.DataFrame(columns=['tau', '-sig(tau)', '+sig(tau)', 'dcf', '-err(dcf)', '+err(dcf)', '#bin'])
    new_acf['dcf'] = np.append(np.flip(acf['dcf']),acf['dcf'])
    new_acf['tau'] = np.append(np.flip(acf['tau']*(-1)), acf['tau'])
    new_acf['+sig(tau)'] = np.append(np.flip(acf['+sig(tau)']*(-1)), acf['+sig(tau)'])
    new_acf['-sig(tau)'] = np.append(np.flip(acf['-sig(tau)']*(-1)), acf['-sig(tau)'])
    new_acf['+err(dcf)'] = np.append(np.flip(acf['+err(dcf)']),acf['+err(dcf)'])
    new_acf['-err(dcf)'] = np.append(np.flip(acf['-err(dcf)']),acf['-err(dcf)'])
    new_acf['#bin'] = np.append(np.flip(acf['#bin']),acf['#bin'])
    return new_acf

def interp(a,b):
    # a: df with common grid of tau values
    # b: df which we want to force upon that grid
    
    new_b = pd.DataFrame(columns=['tau', '-sig(tau)', '+sig(tau)', 'dcf', '-err(dcf)', '+err(dcf)', '#bin'])
    
    f0 = interpolate.interp1d(b['tau'], b['dcf'], kind='quadratic',fill_value="extrapolate")
    new_b['dcf'] = f0(a['tau'])
    
    
    for y in ['dcf','-sig(tau)','+sig(tau)','-err(dcf)','+err(dcf)','#bin']:
        if y == 'dcf':
            kind = 'quadratic'
        else:
            kind = 'nearest'
        f = interpolate.interp1d(b['tau'], b[y], kind=kind, fill_value="extrapolate")
        new_b[y] = f(a['tau'])
    
    new_b['tau'] = a['tau']
        
    return new_b


def delta_ccf(acf, ccf):
    """
    Subtract ACF from CCF and use error propagation (Laursen et al. 2019) to estimate asymmetric errors
    in the resulting function.
    """
    
    delta = pd.DataFrame(columns=['tau', '-sig(tau)', '+sig(tau)', 'dcf', '-err(dcf)', '+err(dcf)', '#bin'])
    
    delta['tau'] = acf['tau']
    delta['#bin'] = acf['#bin']

    for i in ccf.index:
        x0 = [ccf.loc[i,'dcf'],-acf.loc[i,'dcf']]
        s1 = [ccf.loc[i,'-err(dcf)'], acf.loc[i,'-err(dcf)']]
        s2 = [ccf.loc[i,'+err(dcf)'], acf.loc[i,'+err(dcf)']]
        delta_val, err1, err2 = add_asym(x0,s1,s2,order=1)
        delta.loc[i,'dcf'] = delta_val
        delta.loc[i,'-err(dcf)'] = err1
        delta.loc[i,'+err(dcf)'] = err2
        
    for i in ccf.index:
        x0 = [ccf.loc[i,'tau'],-acf.loc[i,'tau']]
        s1 = [ccf.loc[i,'-sig(tau)'], acf.loc[i,'-sig(tau)']]
        s2 = [ccf.loc[i,'+sig(tau)'], acf.loc[i,'+sig(tau)']]
        delta_val, err1, err2 = add_asym(x0,s1,s2,order=1)
        delta.loc[i,'-sig(tau)'] = err1
        delta.loc[i,'+sig(tau)'] = err2
    
    # Change data types (must be done for PLIKE to work)
    delta['-sig(tau)'] = delta['-sig(tau)'].astype('int64')
    delta['+sig(tau)'] = delta['+sig(tau)'].astype('int64')
    delta['dcf'] = delta['dcf'].astype('float64')
    delta['-err(dcf)'] = delta['-err(dcf)'].astype('float64')
    delta['+err(dcf)'] = delta['+err(dcf)'].astype('float64')
    delta['#bin'] = delta['#bin'].astype('int64')
        
    return delta


def plot_ccf_acf(delta, ccf, acf, locator=10, save=False, peak=False, tau=0, err_low=0, err_high=0, x1=-20, x2=80, y1=-0.25, y2=1.25):
    """
    Plot CCF, ACF and their difference. Optionally, you can add peak location (tau) and associated errors.
    
    Parameters:
    -----------
    delta: pd.DataFrame
    	Table returned by delta_ccf() function. Contains information about CCF-ACF.
    ccf: pd.DataFrame
    	Table containing information about CCF.
    acf: pd.DataFrame
    	Table containing information about ACF.
    locator: int, default=10
        Parameter indicating a step for plotting x-axis ticks.
    save: bool, default=False
        If True, the resulting plot will be saved as ccf-acf.pdf
    peak: bool, default=False
        If True, plot will show peak location (but you also need to specify 'tau' keyword parameter)
    tau: float, default=0
    	Time-lag value.
    err_low: float, default=0
    	Time-lag (-) error. Enter with minus sign.
    err_high: float, default=0
    	Time-lag (+) error.
    x1: int
    	Lower limit on x-axis of both plots.
    x2: int
    	Upper limit on x-axis of both plots.
    y1: int
    	Lower limit on y-axis of both plots.
    y2: int
    	Upper limit on y-axis of both plots.
    """

    
    fig = plt.figure(figsize=(9,6))
    plt.subplots_adjust(left=None, bottom=None, right=None, top=0.95, wspace=0,hspace=0)
    
    
    ax1 = fig.add_subplot(211)
    
    ax1.plot(acf['tau'], acf['dcf'], 'o-r', markersize=3, linewidth=1, label='ACF')
    ax1.plot(ccf['tau'], ccf['dcf'], 'o--k', markersize=3, linewidth=1, label='CCF')
    ax1.set_xlabel("Time")
    ax1.grid(which='major', axis='x', linestyle='--')
    ax1.xaxis.set_major_locator(plt.MultipleLocator(locator))
    ax1.yaxis.set_major_locator(plt.MultipleLocator(0.2))
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.legend(loc='lower left', fontsize=13)
    #ax1.set_xlim(lim1,lim2)
    ax1.set_xlim(x1,x2)
    ax1.set_ylim(y1,y2)
    
    
    ax2 = fig.add_subplot(212)
    
    lower_error_x =  delta['-sig(tau)'] 
    upper_error_x =  delta['+sig(tau)']
    asymmetric_error_x = np.array(list(zip(lower_error_x, upper_error_x))).T
    
    lower_error_y =  delta['-err(dcf)']
    upper_error_y =  delta['+err(dcf)']
    asymmetric_error_y = np.array(list(zip(lower_error_y, upper_error_y))).T
    
    ax2.errorbar(delta['tau'], delta['dcf'], xerr=asymmetric_error_x, fmt='o-k',
                 markersize=2, linewidth=0.5, label= 'CCF - ACF', ecolor='gray', capsize=2)
    if peak == True:
        nearest_idx = np.abs(delta['tau'] - tau).argmin()
        ax2.vlines(tau, ymin=-1, ymax= delta['dcf'][nearest_idx],linestyles='dashed', colors='royalblue')
        ax2.axvspan(tau+err_low, tau+err_high, alpha=0.2)
    
    #ax2.plot(delta_ccf_x, delta_ccf_y , 'o-k', markersize=2, linewidth=0.5, label='CCF - ACF')
    
    ax2.set_xlabel("Time (days)", fontsize=17, labelpad=8)
    ax2.grid(which='major', axis='x', linestyle='--')
    ax2.xaxis.set_major_locator(plt.MultipleLocator(locator))
    ax2.yaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax2.legend(fontsize=13, loc='lower left')
    ax2.set_ylim(-0.25,0.25)
    ax2.set_xlim(x1,x2)
    if peak==True:
        ax2.text(0.72,0.05, r'$\tau$ = {:.1f} ({:.1f},+{:.1f}) d'.format(tau,err_low,err_high), transform=ax2.transAxes, size=12.5)
    
    fig.text(0.03, 0.5, "Correlation (A.U.)", va='center', rotation='vertical',fontsize=18)
    if save==True:
        plt.savefig('ccf-acf.pdf',dpi=800)
    plt.show()
    
