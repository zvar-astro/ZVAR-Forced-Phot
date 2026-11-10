import os
import sys
import json
import time
import pickle
import logging
import argparse
import requests
import subprocess
import numpy as np
import pandas as pd

from glob import glob
from functools import partial
from multiprocessing import Pool
from multiprocessing.pool import ThreadPool
from multiprocessing.context import TimeoutError
from astropy.io import fits, ascii
from astropy import wcs
from astropy.time import Time
from penquins import Kowalski

# Set some environment variables
BASE_DIR = '/data'
os.environ['NUMEXPR_MAX_THREADS'] = '32'
sys.path.append(f'{BASE_DIR}/ZVAR-Forced-Phot/psf')

# Import from local scripts/modules
from PSF_Phot_PM import get_phot
from psf_dophot import get_phot_ref_image
from build_matchfile import save_to_file
from metric_suite import calc_lc_metrics, get_scope_features

"""
PSF forced photometry on subtracted images

Example:
python3 psf-phot.py 760 14 3 r --jd_start 2458200 --jd_end 2458217.9
python3 psf-phot.py 760 14 3 r --jd_start 2458200 --jd_end 2458600

P. Mroz @ Caltech, 28 Nov 2020

Modified by Z. Vanderbosch @ Caltech, Apr 2022 - Aug 2024
"""

ZTF_BASEURL = 'https://irsa.ipac.caltech.edu/ibe/data/ztf/products'
DATA_DIR = f'{BASE_DIR}/data'
RES_DIR = f'/matchfiles/data/res'
LOG_DIR = f'{BASE_DIR}/ZVAR-Forced-Phot/logs/log_files'
BUILD_DIR = '{BASE_DIR}/ZVAR-Forced-Phot/matchfiles'
LOGIN_URL = "https://irsa.ipac.caltech.edu/account/signon/login.do"
EPOCH_NCOL = 21 # Number of columns to expect in epochs.txt files

PUBLIC_DATE = 2460071.5 # public data collected before 2023-05-07 (DR18)
PRIVATE_DATE = 2460157.5 # private data collected before 2023-08-01

# Fields for which light curve metrics will be calculated.
# Currently none, but as an example below are the original
# 20 SCoPE fields.
METRIC_FIELDS = []
#METRIC_FIELDS = [
#    296,297,423,424,487,488,562,563,682,683,
#    699,700,717,718,777,778,841,842,852,853
#]

def get_image_list (mel, field_id, ccd_id, quad_id, jd_start=2450000.0, 
    jd_end=2470000.0, filters=[1,2,3], verbose=False):
    """ 
    Query the Kowalski DB to get a list of images of a given field.
    
    Parameters
    ----------
    mel : Kowalski object
        a Melman instance
    field_id : int
        ZTF field number
    ccd_id : int
        CCD chip number (1-16)
    quad_id : int
        CCD quadrant number (1-4)
    jd_start : float, optional
        minimal Julian Date (default 2450000.0)
    jd_end : float, optional
        maximal Julian Date (default 2460000.0)
    filters : array_like
        a list of ZTF filters (1 - g, 2 - r, 3 - i), default: [1,2,3]
    verbose : bool
        if True, prints names of found images to the std out
    
    Returns
    -------
    frames, jd : lists
        lists of found images and their JDs
    """

    q = {
        'query_type': 'find',
        'query' : {
            'catalog': 'ZTF_ops',
            'filter': {
                'field': field_id,
                'jd_start': {'$gte': jd_start},
                'jd_end': {'$lte': jd_end},
                'filter': {'$in': filters}
            },
            'projection': {'_id': 0}
        }
    }

    r = mel.query(query=q)
    
    if r['status'] != 'success':
        print('Error executing the Kowalski query.')
        return [], []

    frames, jd = [], []

    if verbose: print('The following images were found:')
    if verbose: print('Field filter pid fileroot jd_start')
    for t in r['data']:
        if str(t['fileroot']) == 'nan':
            continue
        if t['pid'] == 0: 
            continue
        if (t['pid'] == 1) & (t['jd_start'] >= PUBLIC_DATE): 
            continue
        if verbose: 
            print('%4d %1d %s %1s %13.5f'%(t['field'],t['filter'],t['pid'],\
                t['fileroot'],t['jd_start']))
        frames.append(t['fileroot'])
        jd.append(t['jd_start']+0.5*t['exp']/86400.)
        
    return frames, jd

def get_cookie(username, password):
    """Get a cookie from the IPAC login service

    Parameters
    ----------
    username: `str`
        The IPAC account username
    password: `str`
        The IPAC account password
    """
    cfile = os.path.join(os.path.expanduser('~'),'.ztf_cookie')
    if not os.path.exists(cfile):
        url = '%s?josso_cmd=login&josso_username=%s&josso_password=%s' \
            %(LOGIN_URL,username,password)
        response = requests.get(url)
        cookies = response.cookies
        with open(cfile,'wb') as out:
            pickle.dump(cookies, out)
    else:
        file = open(cfile,'rb')
        cookies = pickle.load(file)
    return cookies

def load_file(url, localdir='/tmp', auth=('username','password')):
    """Load a file from the specified URL and save it locally.

    Parameters
    ----------
    url : `str`
        The URL from which the file is to be downloaded
    localdir : `str`
        The local directory to which the file is to be saved
    auth : tuple
        A tuple of (username, password) to access the ZTF archive
    """
    cookies = get_cookie(auth[0],auth[1])
    response = requests.get(url,stream=True,cookies=cookies)
    if response.status_code == 200:
        file = os.path.join(localdir,url[url.rindex('/') + 1:])
        with open(file,'wb') as handle:
            for block in response.iter_content(1024):
                handle.write(block)
        return os.stat(file).st_size
    else:
        return 0
        
def download_data(urls, data_dir, auth=('username','password'), threads=16):
    """
    Download images from the archive using multiple processes.
    
    Parameters
    ----------
    urls : array_like
        List of image URLs
    data_dir : string
        Directory to store the data
    auth : tuple
        A tuple of (username, password) to access the ZTF archive
    threads : int, optional
        Number of parallel threads
    """
    
    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)

    # Check if any images already exist
    new_urls = []
    for u in urls:
        imname = u.split("/")[-1]
        if not os.path.exists(f"{data_dir}/{imname}"):
            new_urls.append(u)
    if len(new_urls) == 0:
        return
        
    # Perform download
    with Pool(threads) as p:
        for job in p.imap_unordered(partial(load_file,localdir=data_dir,auth=auth),new_urls):
            pass   


def download_ZTF_images (kow, field_id, ccd_id, quad_id, data_dir, 
    jd_start=2450000.0, jd_end=2460000.0, filters=[1,2,3], verbose=False,
    threads=16, irsa_auth=('username','password'), existing_phot=[]):
    """
    Download ZTF images of a given field / ccd / quadrant
    
    Parameters
    ----------
    kow : Kowalski object
        a Kowalski instance
    field_id : int
        ZTF field number
    ccd_id : int
        CCD chip number (1-16)
    quad_id : int
        CCD quadrant number (1-4)
    data_dir : string
        Directory to store the data
    jd_start : float, optional
        minimal Julian Date (default 2450000.0)
    jd_end : float, optional
        maximal Julian Date (default 2460000.0)
    filters : array_like
        a list of ZTF filters (1 - g, 2 - r, 3 - i), default: [1,2,3]
    verbose : bool
        if True, prints names of found images to the std out
    threads : int, optional
        Number of parallel threads
    irsa_auth : tuple
        A tuple of (username, password) to access the ZTF archive
    existing_phot : list
        List of frames that have already been processed
    
    Returns
    -------
    frames, jd : lists
        lists of found images and their JDs
    """
    
    # Query Kowalski and get names of the images
    frames, jd = get_image_list(kow,field_id,ccd_id,quad_id,
        jd_start=jd_start,jd_end=jd_end,filters=filters,verbose=verbose)

    # Update frames and JDs based on existing photometry
    jd = [y for x,y in zip(frames,jd) if x not in existing_phot]
    frames = [x for x in frames if x not in existing_phot]

    # Download subtracted images
    urls = ['%s/sci/%s/%s/%s/%s_c%02d_o_q%1d_scimrefdiffimg.fits.fz'%(ZTF_BASEURL,\
        x[4:8],x[8:12],x[12:18],x,ccd_id,quad_id) for x in frames]
    print('Downloading %d Difference images...'%(len(urls)))
    
    download_data(urls,data_dir,threads=threads,auth=irsa_auth)
    
    # Download PSF images
    urls = ['%s/sci/%s/%s/%s/%s_c%02d_o_q%1d_sciimgdaopsfcent.fits'%(ZTF_BASEURL,\
        x[4:8],x[8:12],x[12:18],x,ccd_id,quad_id) for x in frames]
    print('Downloading %d PSF Files..'%(len(urls)))
    
    download_data(urls,data_dir,threads=threads,auth=irsa_auth)
        
    # Download science images
    urls = ['%s/sci/%s/%s/%s/%s_c%02d_o_q%1d_sciimg.fits'%(ZTF_BASEURL,\
        x[4:8],x[8:12],x[12:18],x,ccd_id,quad_id) for x in frames]
    print('Downloading %d Science images...'%(len(urls)))
    
    download_data(urls,data_dir,threads=threads,auth=irsa_auth)
        
    return frames, jd


def get_image_info (hdr):
    """
    Get image stats. This function is based on Kevin's code.
    
    Parameters
    ----------
    hdr : Header instance
        FITS file header
    
    Returns
    -------
    image_info: list
        List containing various image metadata parameters
    """
    
    pid = hdr['PROGRMID']
    seeing = hdr['SEEING']
    magzp = hdr['MAGZP']
    magzpunc = hdr['MAGZPUNC']
    magzprms = hdr['MAGZPRMS']
    clrcoeff = hdr['CLRCOEFF']
    clrcounc = hdr['CLRCOUNC']
    airmass = hdr['AIRMASS']
    nps1matches = hdr['NMATCHES']
    imstat = hdr['STATUS']
    infobits = hdr['INFOBITS']
    limmag = hdr['MAGLIM']
    exptime = hdr['EXPOSURE']
    field = hdr['FIELDID']
    ccd = hdr['CCDID']
    quad = hdr['QID']
    qcomm = hdr['QCOMMENT']
    qcomm = qcomm.replace(" ","_")
            
    image_info = [
        exptime, pid, field, ccd, quad, imstat, infobits, 
        seeing, magzp, magzpunc, magzprms, clrcoeff, clrcounc, 
        limmag, airmass, nps1matches, qcomm
    ]

    return image_info
    
    
def process_single_image_psf (result_file, diff_image, sci_image, 
    psf_image, ra, dec, magzp_fid=26.275, stamp=3, gain=6.2):
    """
    Calculate forced PSF photometry on a single ZTF image.
    
    Parameters
    ----------
    result_file : string
        Name of the result file
    diff_image : string
        Subtracted image name
    sci_image : string
        Science image name
    psf_image : string
        PSF image name
    ra : array_like
        Array of right ascension (deg)
    dec : array_like
        Array of declination (deg)
    magzp_fid : float, optional
        Fiducial magnitude zero-point, default 26.275
    stamp : int, optional
        Stamp size for the PSF photometry, default 3
    gain : float, optional
        Gain on the science image, default 6.2
        
    Returns
    -------
    jd - JD of observation
    result_file - name of the result file
    """

    #print('Processing file %s'%diff_image)
    
    """
    Open the subtracted image and performed PSF photometry with fixed
    centroids.
    """

    try:
        with fits.open(diff_image) as hdul:
            
            """Open image header"""
            
            hdr = hdul[1].header
            fil = hdr['FILTER'].replace(' ','_')
            
            """Transform equatorial coordinates to pixel coordinates
            and select stars located within the image (otherwise
            BasicPSFPhotometry will raise NoOverlapErrror)."""
            
            image_wcs = wcs.WCS(hdr)
            xy = image_wcs.all_world2pix(np.vstack((ra,dec)).T,0)
            x = xy[:,0]
            y = xy[:,1]
            f = (x > 0) & (x < hdr['NAXIS1']) & (y > 0) & (y < hdr['NAXIS2'])
            
            if sum(f) == 0:
                #print('Object Not In the Image',diff_image)
                return None

            image_info = get_image_info(hdr)

            """ Skip over poor-quality images """
            ZPthres = 25.8
            if fil == 'ZTF_r':
                if any([
                    image_info[6] != 0, # infobits
                    image_info[5] < 1, # status
                    image_info[7] > 3.0, # seeing
                    image_info[0] > 31.0, # exptime
                    image_info[8] > 26.65-0.15*image_info[14], # zpmag+airmass
                    image_info[8] < 25.80-0.15*image_info[14], # zpmag+airmass
                    image_info[10] > 0.05, # magzprms
                    image_info[15] < 120]): #nps1matches
                    return None
            if fil == 'ZTF_g':
                if any([
                    image_info[6] != 0,
                    image_info[5] < 1,
                    image_info[7] > 3.0,
                    image_info[0] > 31.0,
                    image_info[8] > 26.70-0.20*image_info[14],
                    image_info[8] < 25.80-0.15*image_info[14], 
                    image_info[10] > 0.06,
                    image_info[15] < 80]): 
                    return None
            if fil == 'ZTF_i':
                if any([
                    image_info[6] != 0,
                    image_info[5] < 1,
                    image_info[7] > 3.0,
                    image_info[0] > 31.0,
                    image_info[8] > 26.00-0.07*image_info[14],
                    image_info[8] < 25.80-0.07*image_info[14], 
                    image_info[10] > 0.06, 
                    image_info[15] < 100]): 
                    return None
                
            jd = hdr['OBSJD'] + 0.5*hdr['EXPTIME']/86400.0
            n_stars = len(ra)
            
            # Perform forced photometry
            flux, ferr = get_phot(
                diff_image,
                sci_image,
                psf_image,
                x+0.5,
                y+0.5,
                gain,
                n_stars,
                stamp
            )
            
            magzp = image_info[8]
            flux *= pow(10.0,0.4*(magzp_fid-magzp))
            ferr *= pow(10.0,0.4*(magzp_fid-magzp))
            
            res = np.vstack((flux,ferr))
            
            res.astype('<f').tofile(result_file)

            return [jd] + image_info + [result_file]
            
    except (FileNotFoundError,OSError):
        #print('File Not Found',diff_image)
        pass
        
    return None


def process_images_psf (results, diff_images, sci_images, psf_images, 
    ra, dec, magzp_fid=26.275, stamp=5, threads=16, gain=6.2):
    """
    Parameters
    ----------
    results : list of string
        Names of result files
    diff_images : list of strings
        Subtracted image names
    sci_images : list of strings
        Science image names
    psf_images : list of strings
        PSF image names
    ra : array_like
        Array of right ascension (deg)
    dec : array_like
        Array of declination (deg)
    magzp_fid : float, optional
        Fiducial magnitude zero-point, default 26.275
    stamp : int, optional
        Stamp size for the PSF photometry, default 3
    gain : float, optional
        Gain on the science image, default 6.2
    threads : int, optional
        Number of parallel threads, default 16
        
    Returns
    -------
    jd - array of JDs of observations
    outfiles - array of result files
    """

    N = len(diff_images)

    all_image_info = []
    with Pool(threads) as pool:
        
        res = [
            pool.apply_async(
                process_single_image_psf,
                args=(
                    results[i],
                    diff_images[i],
                    sci_images[i],
                    psf_images[i],
                    ra,
                    dec,
                    magzp_fid,
                    stamp,
                    gain
                )
            ) for i in range(N)
        ]
        
        for result in res:
            try:
                image_info = result.get(timeout=120)
            except TimeoutError:
                continue
            if image_info is None: 
                continue
            all_image_info.append(image_info)

    field_names = ['jd','exptime','pid','field','ccd','quad',
        'imstat','infobits','seeing','magzp','magzpunc',
        'magzprms','clrcoeff','clrcounc','limmag','airmass',
        'nps1matches','qcomment','result_file']

    all_image_info = pd.DataFrame(all_image_info,columns=field_names)
            
    return all_image_info



def get_image_info_preprocessed (mel, field, ccd, quad, fid, 
    proc_frames, res_dir, data_dir, threads=16, 
    irsa_auth=('username','password')):
    """
    Get image metadata for images that have already 
    been processed but don't have the correct metadata
    already in the epochs.txt file. Uses Kowalski and
    IRSA tables when possible, and downloads the diff
    images if not in tables.

    """

    # ZTF_ops query, no date constraints this time
    ops_query = {
        'query_type': 'find',
        'query' : {
            'catalog': 'ZTF_ops',
            'filter': {
                'field': int(field),
                'filter': int(fid),
            },
            'projection': {'_id': 0}
        }
    }

    # ZTF_exposures query for image metadata
    exp_query = {
        'query_type': 'find',
        'query' : {
            'catalog': 'ZTF_exposures_20230309',
            'filter': {
                'field': int(field),
                'ccd': int(ccd),
                'quad': int(quad),
                'filter': int(fid),
            },
            'projection': {
                '_id': 0
            }
        }
    }

    ops_response = mel.query(query=ops_query)
    ops_data = ops_response.get("data")
    ops_frames = [x['fileroot'] for x in ops_data]

    exp_response = mel.query(query=exp_query)
    exp_data = exp_response.get("data")
    exp_fracdays = [x['fracday'] for x in exp_data]

    # Use IRSA API to get seeing data
    base_search_url = 'https://irsa.ipac.caltech.edu/ibe/search/ztf/products/sci?'
    where_clause = 'WHERE=field={}%20AND%20ccdid={}%20AND%20qid={}%20AND%20fid={}'.format(
        field,ccd,quad,fid
    )
    pars = '&COLUMNS=filefracday,seeing&ct=csv'
    irsa_scisearch_url = base_search_url + where_clause + pars
    response = requests.get(irsa_scisearch_url)
    irsa_table = ascii.read(response.text).to_pandas()
    irsa_fracdays = list(irsa_table.filefracday)


    preproc_info = []
    missing_frames = []

    for f in proc_frames:

        ops_index = ops_frames.index(f)
        ops_entry = ops_data[ops_index]
        fracday = int(f.split("_")[1])

        if (fracday in exp_fracdays) & (fracday in irsa_fracdays):

            exp_index = exp_fracdays.index(fracday)
            exp_entry = exp_data[exp_index]
            irsa_entry = irsa_table[irsa_table.filefracday == fracday]

            # Get metadata
            imstat = 1
            pid = exp_entry['programid']
            seeing = irsa_entry['seeing'].iloc[0]
            magzp = exp_entry['magzp']
            magzpunc = exp_entry['magzpunc']
            magzprms = exp_entry['magzprms']
            clrcoeff = exp_entry['clrcoeff']
            clrcounc = exp_entry['clrcounc']
            airmass = exp_entry['airmass']
            nps1matches = exp_entry['nps1matches']
            infobits = exp_entry['infobits']
            limmag = exp_entry['limitmag']
            exptime = exp_entry['exptime']
            field = exp_entry['field']
            ccd = exp_entry['ccd']
            quad = exp_entry['quad']
            jd = ops_entry['jd_start']+0.5*exptime/86400.
            qcomm = ops_entry['qcomment']
            qcomm = qcomm.replace(" ","_")
            res_name = "{:s}{:s}_c{:02d}_o_q{:1d}_psf.cat".format(
                res_dir,f,ccd,quad
            )

            image_info = [
                jd, exptime, pid, field, ccd, quad, imstat, 
                infobits, seeing, magzp, magzpunc, magzprms,
                clrcoeff, clrcounc, limmag, airmass, nps1matches,
                qcomm, res_name
            ]
            preproc_info.append(image_info)
        else:
            missing_frames.append(f)

    Nmissed = len(missing_frames)
    if Nmissed > 0:

        print('Downloading %d missing difference images...'%(Nmissed))

        # Download subtracted images
        missing_urls = ['%s/sci/%s/%s/%s/%s_c%02d_o_q%1d_scimrefdiffimg.fits.fz'%(
            ZTF_BASEURL,x[4:8],x[8:12],x[12:18],x,ccd,quad) 
            for x in missing_frames]

        download_data(missing_urls,data_dir,threads=threads,auth=irsa_auth)

        missing_diff_names = np.array(
            [data_dir+'/%s_c%02d_o_q%1d_scimrefdiffimg.fits.fz'%(x,ccd,quad) 
            for x in missing_frames]
        )
        missing_res_names = np.array(
            [res_dir+'%s_c%02d_o_q%1d_psf.cat'%(x,ccd,quad) 
            for x in missing_frames]
        )

        all_missing_info = []
        for i,mdn in enumerate(missing_diff_names):

            # If image was not downloaded, delete the 
            # psf.cat file for this epoch
            if os.path.exists(mdn):
                with fits.open(mdn) as mdiff:

                    hdr = mdiff[1].header
                    image_info = get_image_info(hdr)
                    jd = hdr['OBSJD'] + 0.5*hdr['EXPTIME']/86400.0

                mdiff_info = [jd] + image_info + [missing_res_names[i]]
                all_missing_info.append(mdiff_info)
            else:
                os.remove(missing_res_names[i])
                with open('./missing_images.txt','a') as mif:
                    mif.write(f"{mdn}\n")

    if Nmissed > 0:
        preproc_info += all_missing_info

    all_preproc_image_info = pd.DataFrame(
        preproc_info,
        columns=['jd','exptime','pid','field','ccd','quad',
            'imstat','infobits','seeing','magzp','magzpunc',
            'magzprms','clrcoeff','clrcounc','limmag','airmass',
            'nps1matches','qcomment','result_file'
        ]
    )

    return all_preproc_image_info


def download_reference_images (field, ccd, quad, filtr, data_dir, threads=16, \
    irsa_auth=('username','password')):
    """
    Download reference images from IPAC
    
    Parameters
    ----------
    field : int
        ZTF field name
    ccd: : int
        ZTF CCD number
    quad : int
        ZTF quadrant number
    pflags : list
        List of three bools e.g. [True, False, False] saying
        which filters (g,r,i) have already been processed
    data_dir : string
        folder to store the data
    threads : int, optional
        Number of parallel threads, default 16
    irsa_auth : tuple
        A tuple of (username, password) to access the ZTF archive
    """
    
    print('Downloading Reference Image...')
    refim_name = 'ztf_%06d_z%s_c%02d_q%1d_refimg.fits'%(field,filtr,ccd,quad)
    refim_path = f'{data_dir}/{refim_name}'
    if field < 1000: # Primary Fields
        urls = ['%s/ref/000/field%06d/z%s/ccd%02d/q%1d/%s'%\
            (ZTF_BASEURL,field,filtr,ccd,quad,refim_name)]
    else:            # Secondary Fields
        urls = ['%s/ref/001/field%06d/z%s/ccd%02d/q%1d/%s'%\
            (ZTF_BASEURL,field,filtr,ccd,quad,refim_name)]

    
    download_data(urls,data_dir,threads=threads,auth=irsa_auth)
    
    return refim_path


        
def get_phot_ref_image_wrapper (args):
    return get_phot_ref_image(*args)
        

def process_ref_images (kow, img_name, data_dir, res_dir, threads=16):        
    """
    Process reference images using DoPhot
    
    Parameters
    ----------
    kow : Kowalski object
        Kowalski instance
    field : int
        ZTF field number
    filtr : string ('r','g','i')
        ZTF filter name
    data_dir : string
        Folder to store the data
    res_dir : string
        Folder to store the results
    threads : int, optional
        Number of parallel threads, default 16
    """

    ref_res_dir = f'{res_dir}/' + img_name.split("/")[-1].strip('_refimg.fits')
    
    args = [(kow, img_name, data_dir, ref_res_dir, passwd)]
    
    """NOTE: Use ThreadPool instead of Pool because python's opensssl
    library (which is used by Kowalski) has some issues with the latter
    See: https://github.com/dmitryduev/ztf-variable-marshal/blob/master/nb/label.ipynb
    """
    with ThreadPool(processes=threads) as pool:
        for job in pool.imap_unordered(get_phot_ref_image_wrapper,args):
            pass

def check_if_reference_images_are_processed (field, ccd, quad, filtr, res_dir, verbose=True):
    
    root = 'ztf_%06d_z%s_c%02d_q%1d'%(field,filtr,ccd,quad)
    path = '%s/%s/'%(res_dir,root)
    
    file1 = path + root + '_starlist.txt'
    file2 = path + root + '_dophot_input.txt'
    file3 = path + 'params.json'
    file4 = path + root + '_dophot_output.txt'
    
    flag1 = os.path.isfile(file1)
    flag2 = os.path.isfile(file2)
    flag3 = os.path.isfile(file3)
    
    if flag1 and flag2 and flag3:
        if verbose:
            print('Reference Images already processed.')
        return True, path, file1, file4
    else:
        if verbose:
            print('Reference Image not yet processed.')
        return False, path, None, None


# New function that checks for existing difference image photometry
def check_if_diffim_photometry_exists (quad_res_dir):
    """
    Check for any pre-existing photometry of ZTF difference images

    Parameters
    ----------
    quad_res_dir : directory storing photometry results for a quadrant

    Returns
    -------
    flag: Boolean
        True if phot files and epoch data are present
    phot_frames: list
        List of ZTF frames corresponding to phot files
    epoch_data: DataFrame
        Epoch data (JD, BJD, pid, phot_filenames)
    """

    phot_files = sorted(glob(quad_res_dir + "*psf.cat"))
    Nphot = len(phot_files)

    if Nphot == 0:
        return False, [], []
    else:
        # Check for an epochs.txt file
        f_epoch = f"{quad_res_dir}epochs.txt"
        if os.path.exists(f_epoch):

            # Figure out how many columns in epochs.txt
            # Old format = 4, new_format = 21
            with open(f_epoch,'r') as ef:
                first_line = ef.readline().strip('\n')
                ecol = len(first_line.split())

            if ecol == 4:
                epoch_data = pd.read_csv(
                    f_epoch,
                    header=None,
                    delim_whitespace=True,
                    names=['jd','bjd','pid','result_file']
                )
            elif ecol == EPOCH_NCOL:
                epoch_data = pd.read_csv(
                    f_epoch,
                    header=None,
                    delim_whitespace=True,
                    names=[
                        'jd','bjd','fid','exptime','pid','field','ccd','quad',
                        'imstat','infobits','seeing','magzp','magzpunc',
                        'magzprms','clrcoeff','clrcounc','limmag','airmass',
                        'nps1matches','qcomment','result_file'
                    ]
                )
            
            epoch_files = epoch_data['result_file'].values
            epoch_files = [quad_res_dir+x.split("/")[-1] for x in epoch_files]
            Nepoch = len(epoch_files)

            # Make sure epoch data and phot_files match
            if (Nepoch == Nphot) & (all([x in epoch_files for x in phot_files])):
                print(f'{Nepoch} Difference Images Previously Processed.')
                phot_frames = [x.split("/")[-1].split("c")[0][:-1] for x in phot_files]
                epoch_data['result_file'] = epoch_files
                return True, phot_frames, epoch_data
            else:

                # Try to reconcile epochs.txt and psf.cat files
                drop_idx = []
                for i,erf in enumerate(epoch_files):
                    if erf not in phot_files:
                        drop_idx.append(i)
                        #print(f'{erf} MISSING psf.cat file!')
                if len(drop_idx) > 0:
                    epoch_data.drop(drop_idx,inplace=True)

                for prf in phot_files:
                    if prf not in epoch_files:
                        os.remove(prf)
                        #print(f'{prf} MISSING in epochs.txt!')
                phot_files = sorted(glob(quad_res_dir + "*psf.cat"))
                phot_frames = [x.split("/")[-1].split("c")[0][:-1] for x in phot_files]
                Nphot = len(phot_frames)

                if Nphot > 0:
                    print(f'{Nphot} Difference Images Previously Processed.')
                    return True, phot_frames, epoch_data
                else:
                    return False, [], []

                # log_msg = f"f{args.field:06d}_c{args.ccd:02d}_q{args.q:1d}_z{args.filtr:1s} FAILED"
                # err_msg = 'ERROR: Number of entries in epochs.txt file ({})\n'.format(Nepoch) + \
                #           '       does not match number of psf.cat files ({})'.format(Nphot)
                # logger.error(f"{log_msg}\n{err_msg}")
                # exit(f"\n{err_msg}")
                # return False, [], [] 

        else:
            return False, [], []


def process_quad (kow, glo, mel, field, ccd, quad, filtr, data_dir, res_dir, \
    jd_start=2450000.0, jd_end=2470000.0, stamp=5, threads=16, \
    irsa_auth=('username','password')):
    """
    Perform forced PSF photometry for an entire quadrant
    
    Parameters
    ----------
    kow : Kowalski object
        a Kowalski instance
    glo : Kowalski object
        a Gloria instance
    mel : Kowalski object
        a Melman instance
    field : int
        ZTF field number
    ccd : int
        CCD chip number (0-15)
    quad : int
        CCD quadrant number (0-3)
    filtr : string
        Filter name {'g', 'r', 'i'}
    data_dir : string
        Directory to save the images
    res_dir : string
        Directory to save the results
    jd_start : float, optional
        minimal Julian Date (default 2450000.0)
    jd_end : float, optional
        maximal Julian Date (default 2460000.0)
    stamp : int, optional
        Stamp size for the PSF photometry, default 3
    threads : int, optional
        Number of parallel threads, default 16
    irsa_auth : tuple
        A tuple of (username, password) to access the ZTF archive
    """
    
    if filtr == 'g':
        filter_name = 'zg'
        filter_id = 1
    elif filtr == 'r':
        filter_name = 'zr'
        filter_id = 2
    elif filtr == 'i':
        filter_name = 'zi'
        filter_id = 3
    else:
        exit('Wrong filter name')


    """Check for Existing Photometry"""

    phot_flag, phot_frames, epoch_data = check_if_diffim_photometry_exists(res_dir)
    n_preproc = len(phot_frames)
        
    """Download images from the ZTF archive"""
    
    time1 = time.time()
    
    frames,_ = download_ZTF_images(kow,field,ccd,quad,\
        data_dir,jd_start=jd_start,jd_end=jd_end,filters=[filter_id],\
        threads=threads,irsa_auth=irsa_auth,existing_phot=phot_frames)

    n_imgs = len(frames)
    
    time2 = time.time()
    
    if n_imgs == 0:
        if not phot_flag:
            exit('No new images or previous photometry to process.')
        else:

            """ Incorporate Previous Photometry """

            if epoch_data.shape[1] == EPOCH_NCOL:
                epoch_data.drop(labels=['bjd','fid'],inplace=True,axis=1)
            else:
                print('Getting preprocessed metadata...')
                epoch_data = get_image_info_preprocessed(
                    mel, 
                    field, 
                    ccd, 
                    quad, 
                    filter_id, 
                    phot_frames, 
                    res_dir,
                    data_dir,
                    threads=threads,
                    irsa_auth=irsa_auth
                )

            jd = epoch_data['jd'].values
            pid = epoch_data['pid'].values
            output_files = epoch_data['result_file'].values

            # Exclude epochs outside the date limits
            date_limits = (
                (jd > jd_start) &
                (((pid == 1) & (jd < PUBLIC_DATE)) |
                 ((pid >  1) & (jd < jd_end)))
            )
            image_metadata = epoch_data[date_limits].copy(
                deep=True).reset_index(drop=True)

            # Sort metadata by date
            image_metadata.sort_values(by='jd',inplace=True)

            # Remove any photometry files outside the date limits
            # to keep the epochs.txt file and phot files consistent
            for of in output_files[~date_limits]:
                os.remove(of)

            # Load in params.json file
            fname = res_dir + "params.json"
            with open(fname) as fp:
                params = json.load(fp)

            # Load in starlist file
            source_id, ra, dec = np.loadtxt(params['starlist'],unpack=True,
                usecols=(0,1,2),dtype='S64,float,float',encoding='utf-8')
            
            # Populate time params
            time3 = time.time()
            time4 = time.time()

    else:
    
        print('Processing new images...')

        """List images"""
        
        diff_image_names = np.array(
            [data_dir+'/%s_c%02d_o_q%1d_scimrefdiffimg.fits.fz'%(x,ccd,quad) 
            for x in frames]
        )
        psf_image_names = np.array(
            [data_dir+'/%s_c%02d_o_q%1d_sciimgdaopsfcent.fits'%(x,ccd,quad) 
            for x in frames]
        )
        sci_image_names = np.array(
            [data_dir+'/%s_c%02d_o_q%1d_sciimg.fits'%(x,ccd,quad) 
            for x in frames]
        )
        result_names = np.array(
            [res_dir+'%s_c%02d_o_q%1d_psf.cat'%(x,ccd,quad) 
            for x in frames]
        )
        
        # Load params.json file
        fname = res_dir + 'params.json'
        with open(fname) as fp:
            params = json.load(fp)
        
        # Load starlist file
        source_id, ra, dec = np.loadtxt(params['starlist'],unpack=True,
            usecols=(0,1,2),dtype='S64,float,float',encoding='utf-8')
        
        time3 = time.time()
        
        """Perform PSF photometry with fixed centroids"""
        
        image_metadata = process_images_psf(
            result_names,
            diff_image_names,
            sci_image_names,
            psf_image_names,
            ra, dec,
            magzp_fid=params['magzp'],
            stamp=stamp,
            threads=threads,
            gain=params['gain']
        )

        n_imgs_processed = len(image_metadata)
        print(f'{n_imgs_processed} new images passed quality checks.')
        
        time4 = time.time()

        """Incorporate previous photometry if available"""

        if phot_flag:

            if epoch_data.shape[1] == EPOCH_NCOL:
                epoch_data.drop(labels=['bjd','fid'],inplace=True,axis=1)
            else:
                print('Getting preprocessed metadata...')
                epoch_data = get_image_info_preprocessed(
                    mel, 
                    field, 
                    ccd, 
                    quad, 
                    filter_id, 
                    phot_frames, 
                    res_dir,
                    data_dir,
                    threads=threads,
                    irsa_auth=irsa_auth
                )

            jd_preproc = epoch_data['jd'].values
            pid_preproc = epoch_data['pid'].values
            output_files_preproc = epoch_data['result_file'].values

            # Exclude epochs outside the date limits
            date_limits = (
                (jd_preproc > jd_start) &
                (((pid_preproc == 1) & (jd_preproc < PUBLIC_DATE)) |
                 ((pid_preproc >  1) & (jd_preproc < jd_end)))
            )
            epoch_data = epoch_data[date_limits].copy(
                deep=True).reset_index(drop=True)

            # Combine new and preprocessed metadata
            image_metadata = pd.concat([image_metadata,epoch_data])

            # Sort metadata by date
            image_metadata.sort_values(by='jd',inplace=True)


            # Remove any photometry files outside the date limits
            # to keep the epochs.txt file and phot files consistent
            for ofp in output_files_preproc[~date_limits]:
                os.remove(ofp)


    
    params['field'] = field
    params['ccd'] = ccd
    params['quad'] = quad
    params['filter'] = filter_name
    params['res_dir'] = res_dir

    num_images = len(image_metadata)
    image_metadata.insert(1,'filterid',[filter_id]*num_images)

    matchfile_name = save_to_file(image_metadata,params)

    time5 = time.time()

    """ Get SCOPE Features """
    
    if field in METRIC_FIELDS:
        print('Querying SCoPE Periods/Metrics...')
        if os.path.isfile(params['res_dir']+'epochs.txt'):
            scope_data = get_scope_features(
                glo,
                params,
                field,
                ccd,
                quad,
                qlimit=10000
            )

    """ Calculate light curve metrics """

    time6 = time.time()

    if field in METRIC_FIELDS:
        print('Calculating LC metrics...')
        if os.path.isfile(params['res_dir']+'epochs.txt'):
            try:
                lc_metrics = calc_lc_metrics(
                    res_dir,
                    scope_data,
                    params,
                    field,
                    ccd,
                    quad,
                    filtr,
                    percentile=80.,
                    threads=threads
                )
            except Exception as e:
                log_msg = f"f{args.field:06d}_c{args.ccd:02d}_q{args.q:1d}_z{args.filtr:1s} MetricError"
                logger.exception(log_msg)


    """ Save data to an HDF5 file """

    time7 = time.time()
    if os.path.isfile(params['res_dir']+'epochs.txt'):
        subprocess.run(
            [f'{BUILD_DIR}/build_matchfile',
            params['starlist'],
            params['output'],
            params['res_dir']+'lc_metrics.txt',
            params['res_dir']+'epochs.txt',
            str(params['magzp']),
            matchfile_name]
        )

    time8 = time.time()
    
    print('# Downloading data      : %.2f sec'%(time2-time1))
    print('# Doing photometry      : %.2f sec'%(time4-time3))
    print('# SCOPE Query/Crossmatch: %.2f sec'%(time6-time5))
    print('# Calculating metrics   : %.2f sec'%(time7-time6))
    print('# Saving matchfile      : %.2f sec'%(time8-time7))
    print('# JD Start              : %.1f'%(jd_start))
    print('# JD End (Private)      : %.1f'%(jd_end))
    print('# JD End (Public)       : %.1f'%(PUBLIC_DATE))
    print('# Number of stars: %d'%(len(ra)))
    print('# Number of images available: %d'%(n_imgs + n_preproc))
    print('# Number of images processed: %d'%(num_images))
    log_msg = '\n   # Downloading data      : %.2f sec\n'%(time2-time1) + \
              '   # Doing photometry      : %.2f sec\n'%(time4-time3) + \
              '   # SCOPE Query/Crossmatch: %.2f sec\n'%(time6-time5) + \
              '   # Calculating metrics   : %.2f sec\n'%(time7-time6) + \
              '   # Saving matchfile      : %.2f sec\n'%(time8-time7) + \
              '   # Number of stars: %d\n'%(len(ra)) + \
              '   # Number of images available: %d\n'%(n_imgs + n_preproc) + \
              '   # Number of images processed: %d'%(num_images)
    return log_msg
    


# Command line arguments
parser = argparse.ArgumentParser()
parser.add_argument('field',type=int)
parser.add_argument('ccd',type=int)
parser.add_argument('q',type=int)
parser.add_argument('filtr',type=str)
parser.add_argument('--jd_start',type=float,default=2450000.0)
parser.add_argument('--jd_end',type=float,default=PRIVATE_DATE)
parser.add_argument('--logfile',type=str,default=None)
args = parser.parse_args()

# Set up the logger
logger = logging.getLogger(__name__)
if args.logfile is None:
    logfilename = f'{LOG_DIR}/pipeline_{args.field:04d}.log'
else:
    logfilename = args.logfile

logging.basicConfig(
    level=logging.INFO,
    filename=logfilename,
    format='%(asctime)s:%(levelname)s:%(message)s'
)


# Load username/password info
with open(f'{BASE_DIR}/passwd.json','r') as f:
    passwd = json.load(f)

# Connect to Kowalski
k = Kowalski(
    username=passwd['kowalski']['username'],
    password=passwd['kowalski']['password'],
    protocol='https', 
    host='kowalski.caltech.edu', 
    port=443, 
    timeout=240.
)
if not k.ping():
    log_msg = f"f{args.field:06d}_c{args.ccd:02d}_q{args.q:1d}_z{args.filtr:1s} FAILED"
    logger.info(log_msg)
    logger.error('Could not connect to Kowalski!')
    exit('Error while connecting to Kowalski.')

# Connect to Gloria
g = Kowalski(
    token=passwd['gloria']['token'],
    protocol='https',
    host='gloria.caltech.edu',
    port=443,
    timeout=240.
)
if not g.ping():
    log_msg = f"f{args.field:06d}_c{args.ccd:02d}_q{args.q:1d}_z{args.filtr:1s} FAILED"
    logger.info(log_msg)
    logger.error('Could not connect to Gloria!')
    exit('Error while connecting to Gloria.')

# Connect to Melman
m = Kowalski(
    token=passwd['gloria']['token'],
    protocol='https',
    host='melman.caltech.edu',
    port=443,
    timeout=240.
)
if not m.ping():
    log_msg = f"f{args.field:06d}_c{args.ccd:02d}_q{args.q:1d}_z{args.filtr:1s} FAILED"
    logger.info(log_msg)
    logger.error('Could not connect to Melman!')
    exit('Error while connecting to Melman.')

if not os.path.isdir(DATA_DIR):
    os.mkdir(DATA_DIR)
if not os.path.isdir(RES_DIR):
    os.mkdir(RES_DIR)


""" Check whether reference images have already been processed """

flag, quad_res_dir, starlist_fname, output_fname = check_if_reference_images_are_processed(
    args.field,
    args.ccd,
    args.q,
    args.filtr,
    RES_DIR,
    verbose=True)


if not flag:
    
    start_time = time.time()

    # Download quadrant reference image
    ref_path = download_reference_images(
        args.field,
        args.ccd,
        args.q,
        args.filtr,
        DATA_DIR,
        threads=16,
        irsa_auth=(
            passwd['irsa']['username'],
            passwd['irsa']['password']
        )
    )
    if not os.path.exists(ref_path):
        exit('Reference image not available at IRSA.\n')
    
    # Process reference image
    process_ref_images(
        k, 
        ref_path, 
        DATA_DIR, 
        RES_DIR, 
        threads=16
    )
    
    end_time = time.time()
    print('# Processing reference images: %.2f sec'%(end_time-start_time))

    # Update processing flag and starlist/output filenames
    flag,_,starlist_fname, output_fname = check_if_reference_images_are_processed(
        args.field,
        args.ccd,
        args.q,
        args.filtr,
        RES_DIR,
        verbose=False
    )

if flag:
    
    # Check that the params file has correct filenames
    fname = quad_res_dir + 'params.json'
    with open(fname) as fp:
        params = json.load(fp)
    if params['starlist'] != starlist_fname:
        params['starlist'] = starlist_fname
        with open(fname,'w') as fp:
            json.dump(params,fp)
    if params['output'] != output_fname:
        params['output'] = output_fname
        with open(fname,'w') as fp:
            json.dump(params,fp)
    
    start_time = time.time()

    # Perform photometry on quadrant sources
    try:
        pmsg = process_quad(
            k,g,m,
            args.field,
            args.ccd,
            args.q,
            args.filtr,
            DATA_DIR,
            quad_res_dir,
            jd_start=args.jd_start,
            jd_end=args.jd_end,
            threads=16,
            stamp=3,
            irsa_auth=(
                passwd['irsa']['username'],
                passwd['irsa']['password']
            )
        )
        log_msg = f"f{args.field:06d}_c{args.ccd:02d}_q{args.q:1d}_z{args.filtr:1s} SUCCESS"
        logger.info(f"{log_msg}{pmsg}")
    except Exception as e:
        log_msg = f"f{args.field:06d}_c{args.ccd:02d}_q{args.q:1d}_z{args.filtr:1s} FAILED"
        logger.exception(log_msg)


    end_time = time.time()
    print('# Execution time: %.2f sec\n'%(end_time-start_time))
