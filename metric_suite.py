import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import multiprocessing
import matplotlib.pyplot as plt
import scipy.stats as spstats
warnings.filterwarnings("ignore")

from glob import glob
from astropy.coordinates import SkyCoord, Angle
from astropy.visualization import PercentileInterval

"""
Functions for calculative light curve statistical
metrics and querying periodicity information from
the scope_features catalogs on Kowalski machines

Author: Zach Vanderbosch (Caltech)
"""

# Minimum number of valid data points required per source
NLIMIT = 20

# Directory containing zero-point calibration info
CAL_DIR = '/data/ZVAR-Forced-Phot/cals'

FUNC_ORDER = {
    'stdev':0,
    'pp_scat':1,
    'ms_successive_diff':2,
    'skew':3,
    'skewP':4,
    'kurtosis':5,
    'norm_ev':6,
    'von_neumann':7,
    'interq_range':8,
    'welch_i':9,
    'stetson_J':10,
    'stetson_K':11
}

# Standard Deviation Function
def stdev(vals,errs,N):
    mean = np.nanmean(vals, axis=0)
    sd = (np.nansum((vals-mean)**2, axis=0)/N)**0.5
    sd[N<NLIMIT] = np.nan
    return sd,'stdev'

# Point-2-Point Scatter Function
def pp_scat(vals,errs,N):
    v1 = vals[0:-1,:]
    v2 = vals[1:,:]
    vdiff = v2-v1
    pp_avg = (np.nansum(vdiff**2, axis=0)/(N-1))**0.5
    pp_avg[N<NLIMIT] = np.nan
    return pp_avg,'pp_scat'

# Mean Square Successive Difference function
def ms_successive_diff(vals,errs,N):
    v1 = vals[0:-1,:]
    v2 = vals[1:,:]
    vdiff = v2-v1
    msdiff = np.nansum(vdiff**2, axis=0)/(N-1)
    msdiff[N<NLIMIT] = np.nan
    return msdiff,'ms_successive_diff'

# Function to Calculate Skew
def skew(vals,errs,N):
    mean = np.nanmean(vals, axis=0)
    nfactor = abs((N/(N-1)/(N-2)))
    sk = nfactor * np.nansum((vals-mean)**3/(errs**3), axis=0)
    sk[N<NLIMIT] = np.nan
    return sk,'skew'

# Function to Calculate Pearson Skew
def skewP(vals,errs,N):
    mean = np.nanmean(vals, axis=0)
    median = np.nanmedian(vals, axis=0)
    sd,_ = stdev(vals,errs,N)
    sp = (mean-median)/sd
    sp[N<NLIMIT] = np.nan
    return sp,'skewP'

# Function to Calculate kurtosis
def kurtosis(vals,errs,N):
    mean = np.nanmean(vals, axis=0)
    term1 = N*(N+1)/(N-1)/(N-2)/(N-3)
    term2 = np.nansum((mean-vals)**4/(errs**4), axis=0)
    term3 = 3.0*(N-1)**2/(N-2)/(N-3)
    k = (term1 * term2) - term3
    k[N<NLIMIT] = np.nan
    return k,'kurtosis'

# Function to Normalized Excess Variance
def norm_ev(vals,errs,N):
    mean = np.nanmean(vals, axis=0)
    term1 = 1./(N*mean**2)
    term2 = np.nansum((mean-vals)**2 - errs**2, axis=0)
    nev = term1 * term2
    nev[N<NLIMIT] = np.nan
    return nev,'norm_ev'

# Function to Normalized Excess Variance
def von_neumann(vals,errs,N):
    mssd,_ = ms_successive_diff(vals,errs,N)
    sd,_ = stdev(vals,errs,N)
    vari = sd**2
    vn = mssd/vari
    vn[N<NLIMIT] = np.nan
    return vn,'von_neumann'

# Functin to calculate the inter-quartile range
def interq_range(vals,errs,N):
    iqr = spstats.iqr(vals, axis=0, nan_policy='omit')
    iqr[N<NLIMIT] = np.nan
    return iqr,'interq_range'


# Function to calculate the Welch/Stetson I metric
def welch_i(vals,errs,N):
    mean = np.nanmean(vals, axis=0)
    nfactor = N/(N-1)
    term1 = (vals-mean)/errs
    term2 = (vals[::-1,:]-mean)/errs[::-1,:]
    wi = nfactor * np.nansum(term1*term2, axis=0)
    wi[N<NLIMIT] = np.nan
    return wi,'welch_i'

# Function to calculate Stetson J metric
def stetson_J(vals,errs,N):
    mean = np.nanmean(vals, axis=0)
    nfactor = np.sqrt(N/(N-1))
    vdiff1 = vals - mean
    vdiff2 = vals[::-1,:] - mean
    del1 = vdiff1/errs
    del2 = vdiff2/errs[::-1,:]
    vsign = np.sign(vdiff1*vdiff2)
    vroot = np.sqrt(abs(del1*del2))
    J = nfactor * np.nansum(vsign * vroot, axis=0)
    J[N<NLIMIT] = np.nan
    return J,'stetson_J'

# Function to calculate Stetson K metric
def stetson_K(vals,errs,N):
    mean = np.nanmean(vals, axis=0)
    nfactor = np.sqrt(1./N)
    vdiff1 = vals - mean
    vdiff2 = vals[::-1,:] - mean
    del1 = vdiff1/errs
    del2 = vdiff2/errs[::-1,:]
    numer = np.nansum(abs(del1*del2), axis=0)
    denom = np.nansum(del1*(del2**2), axis=0)
    K = nfactor * numer / np.sqrt(denom)
    K[N<NLIMIT] = np.nan
    return K,'stetson_K'

def mag_corr (mag, field, ccd, quad, filtr):
    """
    Parameters
    ----------
    mag : numpy array
        Array of magnitudes
    field : int
        ZTF field number
    ccd : int
        CCD chip number (1-16)
    quad : int
        CCD quadrant number (1-4)
    filtr : string
        Filter name {'g', 'r', 'i'}
    
    Returns
    -------
    calmag : float
        Magnitude calibrated to ZTF system
    """

    caldata_dtype = {
        'names': (
            'field', 'ccd', 'quad', 
            'zp_corr_g', 'n_calib_g', 
            'zp_corr_r', 'n_calib_r', 
            'zp_corr_i', 'n_calib_i'
        ),
        'formats': (
            'i4','i4','i4','f4','i4','f4','i4','f4','i4'
        )
    }

    if field < 1000:
        calfile = f"{CAL_DIR}/zero_point_corrections_primary.txt"
    else:
        calfile = f"{CAL_DIR}/zero_point_corrections_secondary.txt"
    caldata = np.loadtxt(calfile,dtype=caldata_dtype)
        
    cal_index = (
        (caldata['field'] == field) &
        (caldata['ccd'] == ccd) &
        (caldata['quad'] == quad)
    )

    if not any(cal_index):
        return mag
        
    if filtr == 'g':
        magcorr = caldata['zp_corr_g'][cal_index][0]
    elif filtr == 'r':
        magcorr = caldata['zp_corr_r'][cal_index][0]
    elif filtr == 'i':
        magcorr = 0.0

    calmag = mag + magcorr

    return calmag



def flux2mag (diff_flux, diff_flux_unc, magzp, 
    field, ccd, quad, filtr, output):
    """
    Parameters
    ----------
    phot : pandas DataFrame
        DataFrame containing flux/flux_err photometry of target
    magzp : float
        Zero-point magnitude of the reference image
    field : int
        ZTF field number
    ccd : int
        CCD chip number (0-15)
    quad : int
        CCD quadrant number (0-3)
    filtr : string
        Filter name {'g', 'r', 'i'}
    output : string
        DoPHOT output file path+name
    
    Returns
    -------
    abs_mag : float
        Absolute Magnitude (not calibrated to ZTF system)
    abs_mag_unc : float
        Absolute Magnitude error

    """

    """Get photometry results from ref image"""

    Nstars = diff_flux.shape[1]
    ref_mag = output[:Nstars,4]
    ref_mag_unc = output[:Nstars,5]

    ref_flux = 10**(-0.4*ref_mag)
    ref_flux_unc = ref_mag_unc/1.0857 * ref_flux
    total_flux = diff_flux + ref_flux
    total_flux_unc = np.sqrt(diff_flux_unc**2 + ref_flux_unc**2)
    snr = total_flux/total_flux_unc

    """ Convert total flux values to magnitudes """

    SNT = 0.0 # S/N threshold for total flux
    q = snr > SNT
    abs_mag = np.zeros(diff_flux.shape)*np.nan
    abs_mag_unc = np.zeros(diff_flux.shape)*np.nan
    abs_mag[q] = magzp - 2.5*np.log10(total_flux[q])
    abs_mag_unc[q] = (2.5/np.log(10.))*(1./snr[q])

    """ Apply magnitude correction """

    abs_mag = mag_corr(abs_mag, field, ccd, quad, filtr)

    return abs_mag, abs_mag_unc


def get_scope_features(gloria, params, field, ccd, quad, qlimit=10000):
    """
    Parameters
    ----------
    gloria : penquins Kowalski Instance
        A gloria instance
    params : dict
        Saved parameters with filenames
    field : int
        ZTF field number
    ccd : int
        CCD chip number (0-15)
    quad : int
        CCD quadrant number (0-3)
    qlimit : int
        Maximum number of returned entries per Gloria query
    
    Returns
    -------
    scope_data : float
        Array containing cross-matched SCOPE data
    """

    # Load in source coordinates
    source_id, ra, dec = np.loadtxt(
        params['starlist'],
        unpack=True,
        usecols=(0,1,2),
        dtype='S64,float,float',
        encoding='utf-8'
    )
    starlist_coords = SkyCoord(
        ra=ra,
        dec=dec,
        unit='deg',
        frame='icrs'
    )
    Nstars = len(ra)

    # Initialize array
    scope_data = np.zeros((Nstars,3))*np.nan

    # First query to figure out how many 
    # entries are available in quadrant
    doc_query = {
        'query_type': 'count_documents',
        'query' : {
            'catalog': 'ZTF_source_features_DR5',
            'filter': {
                'field': int(field),
                'ccd': int(ccd),
                'quad': int(quad),
            }
        }
    }
    doc_response = gloria.query(query=doc_query)
    doc_number = doc_response['data']
    if doc_number == 0:
        print('No SCOPE data in this quadrant.')
        return scope_data
    print(f'   {Nstars} sources in quadrant')

    # Figure out how many queries are needed
    num_queries = int(np.ceil(doc_number / qlimit))
    skip_values = np.arange(0,num_queries,1,dtype=int)*qlimit

    data = []
    for i in range(num_queries):

        # Define the SCOPE data query
        scope_query = {
            'query_type': 'find',
            'query' : {
                'catalog': 'ZTF_source_features_DR5',
                'filter': {
                    'field': int(field),
                    'ccd': int(ccd),
                    'quad': int(quad),
                },
                'projection': {
                    'ra': 1,
                    'dec': 1,
                    'n': 1,
                    'period': 1,
                    'significance': 1,
                    'median': 1,
                }
             },
            'kwargs': {
                'skip': int(skip_values[i]),
                'limit': int(qlimit)
            }
        }

        # Perform the query
        r = gloria.query(query=scope_query)
        d = r.get("data")
        data += d

        # Be a good samaritan
        time.sleep(1.0)

    Nper = len(data)
    if Nper == 0:
        print('No SCOPE data in this quadrant.')
        return scope_data
           
    # Extract relevant info     
    all_ras = np.array([x.get("ra") for x in data])
    all_decs = np.array([x.get("dec") for x in data])
    all_nump = np.array([x.get("n") for x in data])
    all_periods = np.array([x.get("period") for x in data])
    all_psigs = np.array([x.get("significance") for x in data])
    all_medmags = np.array([x.get("median") for x in data])

    scope_coords = SkyCoord(
        ra=all_ras,
        dec=all_decs,
        unit="deg",
        frame='icrs'
    )

    # Crossmatch starlist coordinates with scope entries
    # match_idx, match_sep2d,_ = starlist_coords.match_to_catalog_sky(scope_coords)
    match_idx_starlist, match_idx_scope,_,_ = scope_coords.search_around_sky(
        starlist_coords,
        seplimit=Angle(1.0, unit='arcsec')
    )

    uidx = np.unique(match_idx_starlist)
    Nmatch = len(uidx)

    # There may be multiple SCOPE entries per source. To simplify
    # things, just choose entry based on the largest number of
    # data points "n"
    scope_match_idx = []
    for ux in uidx:
        idx1 = np.where(match_idx_starlist==ux) 
        idx2 = match_idx_scope[idx1]
        n = all_nump[idx2]
        scope_idx = idx2[n == max(n)][0]
        scope_match_idx.append(scope_idx)
    scope_match_idx = np.array(scope_match_idx,dtype=int)

    print(f'   {Nper} SCOPE entries found')
    print(f'   {Nmatch} sources matched with SCOPE entries')

    # Populate array with matched period & significance values
    scope_data[uidx,0] = all_medmags[scope_match_idx]
    scope_data[uidx,1] = all_periods[scope_match_idx]
    scope_data[uidx,2] = all_psigs[scope_match_idx]

    return scope_data


# Define a Function to calculate all of the ZTF
# variability metrics desired for a single object
def calc_lc_metrics(res_dir, scope_data, params, 
    field, ccd, quad, filtr, percentile=80.,threads=1):

    # List of stat functions
    lc_stat_funcs = [
        stdev, 
        pp_scat, 
        ms_successive_diff, 
        skew, 
        skewP,
        kurtosis, 
        norm_ev, 
        von_neumann,
        interq_range,
        welch_i,
        stetson_J, 
        stetson_K
    ]
    Nstats = len(lc_stat_funcs)


    """ Load in Photometry Data """

    # Get required filenames
    epochs_file = f"{res_dir}epochs.txt"
    starlist_file = params['starlist']
    dophot_file = params['output']
    phot_files = sorted(glob(f"{res_dir}ztf*psf.cat"))

    # Load in epoch and starlist files
    epochs = pd.read_csv(epochs_file, delim_whitespace=True,header=None)
    starlist = np.loadtxt(starlist_file)
    refphot = np.loadtxt(dophot_file)
    Nstars = len(starlist)

    # Raise error if mismatching number of epochs and psf.cat files
    if len(epochs) != len(phot_files):
        print('\nERROR: Number of entries in epochs.txt file ({})'.format(len(epochs)))
        print('       does not match number of psf.cat files ({})\n'.format(len(phot_files)))
        return

    # Initialize array storing SCOPE data and metrics
    scope_cols = scope_data.shape[1]
    lc_stats = np.zeros((Nstars,Nstats+scope_cols+1))
    lc_stats[:,0] = starlist[:,0]
    lc_stats[:,1:scope_cols+1] = scope_data

    # Generate some basic quality flags
    bjd = epochs.iloc[:,1].values
    infobits = epochs.iloc[:,9].values
    seeing = epochs.iloc[:,10].values
    magzp = epochs.iloc[:,11].values
    zprms = epochs.iloc[:,13].values
    limmag = epochs.iloc[:,16].values
    airmass = epochs.iloc[:,17].values
    ncal = epochs.iloc[:,18].values

    q = (
        (infobits == 0) & 
        (seeing < 3.0) &
        (zprms < 0.06) &
        (airmass < 2.5) &
        (ncal > 80)
    )
    Nepochs = len(q[q])
    
    if Nepochs > 5:

        # Initialize matrices
        flux = np.zeros((Nepochs, Nstars))
        eflux = np.zeros((Nepochs, Nstars))
        flags = np.zeros((Nepochs, Nstars))

        # Filter out poor quality data
        phot_files = [f for qval,f in zip(q,phot_files) if qval]

        # Load photometry data
        for i,f in enumerate(phot_files):

            epoch_data = np.fromfile(f, dtype='<f')
            epoch_data = epoch_data.reshape((2,Nstars))
            epoch_flux = epoch_data[0,:]
            epoch_eflux = epoch_data[1,:]
            epoch_nanidx = (
                np.isnan(epoch_flux) | 
                np.isnan(epoch_eflux)
            )

            # Save data to matrices
            flux[i,:] = epoch_flux
            eflux[i,:] = epoch_eflux
            flags[i,:] = epoch_nanidx*1


        # Apply percentile range to data
        percentiles = np.nanpercentile(flux, q=(10.,90.), axis=0)
        flux[flux<percentiles[0,:]] = np.nan
        flux[flux>percentiles[1,:]] = np.nan
        eflux[flux<percentiles[0,:]] = np.nan
        eflux[flux>percentiles[1,:]] = np.nan
        flags[flux<percentiles[0,:]] = 1
        flags[flux>percentiles[1,:]] = 1

        # Get number of non-NaN observations per source
        N = np.sum((1-flags), axis=0)

        # Calculate light curve stats
        for i,func in enumerate(lc_stat_funcs):
            stat_vals, func_name = func(flux, eflux, N)
            stat_index = FUNC_ORDER[func_name]
            lc_stats[:,stat_index+scope_cols+1] = stat_vals
        #with multiprocessing.Pool(threads) as pool:
        #    pool_args = (flux, eflux, N)
        #    res = [
        #        pool.apply_async(func, args=pool_args) 
        #        for func in lc_stat_funcs
        #    ]
        #    for r in res:
        #        stat_vals, func_name = r.get(timeout=120)
        #        stat_index = FUNC_ORDER[func_name]
        #        lc_stats[:,stat_index+scope_cols+1] = stat_vals


    # Save stats to file
    fname = f'{res_dir}lc_metrics.txt'
    fmt = "%-20d %6.3f %10.6f %10.6f %12.6f %12.6f %15.6f %13.6f %9.6f %14.6f %16.6f %9.6f %10.6f %10.6f %10.6f %10.6f\n"
    with open(fname,'w') as fp:
        for i in range(Nstars):
            fp.write(fmt %tuple(lc_stats[i,:]))

    return
