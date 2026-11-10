import numpy as np
import os
import subprocess
import json
from astropy.io import fits
from astropy import wcs
import pandas as pd

"""
Routines for calling dophot on reference images
"""

DOPHOT_DIR = "/data/ZVAR-Forced-Phot/dophot"

def get_info_from_ref_image (filename):
    """
    Get info from a reference image
    
    Parameters
    ----------
    filename : string
        Reference image name
    
    Returns
    -------
    ra_cen, dec_cen - equatorial coordinates of the center of the image
    magzp - magnitude zero-point
    ref_wcs - WCS of the reference image
    """
    
    params = {}
    
    with fits.open(filename) as hdul:
        hdr = hdul[0].header
        ref_img_wcs = wcs.WCS(hdr)
        params['magzp'] = hdr['MAGZP']
        x_cen = hdr['NAXIS1']/2
        y_cen = hdr['NAXIS2']/2
        xx = ref_img_wcs.all_pix2world(np.reshape([x_cen,y_cen],(1,2)),0)
        params['ra_cen'] = xx[0][0]
        params['dec_cen'] = xx[0][1]
        params['nx'] = hdr['NAXIS1']
        params['ny'] = hdr['NAXIS2']
        params['saturate'] = hdr['SATURATE']
        params['fwhm'] = hdr['MEDFWHM']
        params['background'] = hdr['GLOBMED']
    
    return params, ref_img_wcs

def radec_to_xy (ra, dec, image_wcs, origin=0):
    """
    Transform equatorial coordinates to pixel coordinates.
    
    Parameters
    ----------
    ra : array
        Array of RA (deg)
    dec : array
        Array of Decl. (deg)
    image_wcs : astropy WCS object
        Image WCS
    origin : int
        coordinate in the upper left corner image, default 0
        
    Returns
    -------
    x, y: pixel coordinates
    """
    
    xy = image_wcs.all_world2pix(np.vstack((ra,dec)).T,origin)
    x = xy[:,0]
    y = xy[:,1]
    
    return x, y
    
def get_Gaia_stars (kow, ra, dec, radius=0.64):
    """
    Query the Kowalski DB to get a list of Gaia sources.
    
    Parameters
    ----------
    kow : Kowalski object
        a Kowalski instance
    ra : float
        Right ascention (deg)
    dec : float
        Declination (deg)
    radius : float, optional
        Search radius (deg), default: 0.64 deg
    
    Returns
    -------
    source_id, ra, dec : numpy arrays containing Gaia ID, RA, and Decl
    """    
    
    q = {'query_type': 'cone_search',
        'query': {
            'object_coordinates': {
               'radec': '[(%.6f,%.6f)]'%(ra,dec), 
                'cone_search_radius': '%.4f'%radius,
                'cone_search_unit': 'deg'
            },
            'catalogs': {
                'Gaia_EDR3': {
                    'filter': {},
                    'projection': {
                        'ra': 1, 
                        'dec': 1,
                        'source_id': 1
                    }
                }
            }
        }
        }
    
    r = kow.query(query=q)
    
    if r['status'] != 'success':
        print('Error executing the Kowalski query.')
        return np.array([]), np.array([]), np.array([])
    
    source_id, ra_all, dec_all = [], [], []
    key = list(r['data']['Gaia_EDR3'].keys())[0]
    for x in r['data']['Gaia_EDR3'][key]:
        source_id.append(str(x['source_id']))
        ra_all.append(x['ra'])
        dec_all.append(x['dec'])
    
    return np.array(source_id), np.array(ra_all), np.array(dec_all)

def get_PS1_stars (kow, ra, dec, radius=0.64):
    """
    Query the Kowalski DB to get a list of Pan-STARRS sources.
    
    Parameters
    ----------
    kow : Kowalski object
        a Kowalski instance
    ra : float
        Right ascention (deg)
    dec : float
        Declination (deg)
    radius : float, optional
        Search radius (deg), default: 0.64 deg
    
    Returns
    -------
    source_id, ra, dec : numpy arrays containing PS ID, RA, and Decl
    """    
    
    MIN_NUM_PS_DETECTIONS = 5
    
    q = {'query_type': 'cone_search',
        'query': {
            'object_coordinates': {
               'radec': '[(%.6f,%.6f)]'%(ra,dec), 
                'cone_search_radius': '%.4f'%radius,
                'cone_search_unit': 'deg'
            },
            'catalogs': {
                'PS1_DR1': {
                    'filter': {},
                    'projection': {
                        'raMean': 1, 
                        'decMean': 1,
                        '_id': 1,
                        'nDetections': 1,
                        'gMeanPSFMag': 1,
                        'rMeanPSFMag': 1
                    }
                }
            }
        }
        }
    
    r = kow.query(query=q)
    
    if r['status'] != 'success':
        print('Error executing the Kowalski query.')
        return np.array([]), np.array([]), np.array([])
    
    source_id, ra_all, dec_all, psg_all, psr_all = [], [], [], [], []
    key = list(r['data']['PS1_DR1'].keys())[0]
    for x in r['data']['PS1_DR1'][key]:
        if x['nDetections'] <= MIN_NUM_PS_DETECTIONS: continue
        source_id.append(str(x['_id']))
        ra_all.append(x['raMean'])
        dec_all.append(x['decMean'])
        try:
            psg_all.append(x['gMeanPSFMag'])
        except:
            psg_all.append(np.nan)
        try:
            psr_all.append(x['rMeanPSFMag'])
        except:
            psr_all.append(np.nan)
    
    return np.array(source_id), np.array(ra_all), np.array(dec_all), np.array(psg_all), np.array(psr_all)

def create_dophot_parameter_input_file (basename, params, data_dir, result_dir):
    
    fname = os.path.join(result_dir,basename + 'dophot_params.txt')
        
    with open(fname,'w') as fp:
        fp.write('PSFTYPE = \'PGAUSS\'\n')
        fp.write('FWHM = %.2f\n'%params['fwhm'])
        fp.write('SKY = %.1f\n'%params['background'])
        fp.write('EPERDN = %.1f\n'%params['gain'])
        fp.write('RDNOISE = %.1f\n'%params['rdnoise'])
        fp.write('PARAMS_DEFAULT = \'%s/param_default_c\'\n'%(DOPHOT_DIR))
        fp.write('PARAMS_OUT = \'{}/param_out\'\n'.format(result_dir))
        fp.write('IMAGE_IN = \'{}\'\n'.format(data_dir+'/'+basename+'refimg.fits'))
        fp.write('IMAGE_OUT = \'{}\'\n'.format(result_dir+'/'+basename+'refimg_out.fits'))
        fp.write('OBJECTS_IN = \'{}\'\n'.format(result_dir+'/'+basename+'dophot_input.txt'))
        fp.write('OBJTYPE_IN = \'INTERNAL\'\n')
        fp.write('OBJECTS_OUT = \'{}\'\n'.format(result_dir+'/'+basename+'dophot_output.txt'))
        fp.write('OBJTYPE_OUT = \'COMPLETE\'\n')
        fp.write('SHADOWFILE_OUT = \'%s\'\n'%(result_dir+'/shadowfile_out'))
        fp.write('ERRORS_OUT = \'%s\'\n'%(result_dir+'/error_out'))
        fp.write('LOGFILE = \'{}\'\n'.format(result_dir+'/log.txt'))
        fp.write('LOGVERBOSITY = 1\n')
        fp.write('FIXPOS = \'YES\'\n')
        fp.write('THRESHMIN = 20.0\n')
        fp.write('THRESHMAX = 20000.0\n')
        fp.write('CTPERSAT = 6.0e4\n')
        fp.write('ITOP = %.1f\n'%(params['saturate']))
        fp.write('CENTINTMAX = %.1f\n'%(params['saturate']))
        fp.write('END\n')


def get_phot_ref_image (kow, imagename, data_dir, result_dir, passwd):
    
    if not os.path.isdir(result_dir):
        os.mkdir(result_dir)
    
    basename = imagename.split('/')[-1].replace('refimg.fits','')
    
    """Grab info from the image header: coordinates of the center 
    of the field, magnitude zero-point, WCS
    """
    
    params, ref_img_wcs = get_info_from_ref_image(imagename)

    """Get a list of objects from Gaia and transform to image coordinates
    """

    print('Performing Pan-STARRS1 Query...')
    source_id, ra, dec, mg, mr = get_PS1_stars(kow,params['ra_cen'],\
        params['dec_cen'],radius=0.64)

    x, y = radec_to_xy(ra,dec,ref_img_wcs)

    """Select stars located within the reference image footprint"""

    f = (x > 0) & (x < params['nx']) & (y > 0) & (y < params['ny'])
    
    source_id = source_id[f]
    ra = ra[f]
    dec = dec[f]
    x = x[f]
    y = y[f]
    mg = mg[f]
    mr = mr[f]
    
    nstars = sum(f)
    print(f'{nstars} PS1 Sources on Reference Image.')


    ps_output = {
        'source_id': source_id,
        'ra':ra,
        'dec':dec,
        'gmag':mg,
        'rmag':mr
    } 
    ps_df = pd.DataFrame(ps_output)
    ps_df.to_csv(f"{result_dir}/ps_sources.csv",index=False)
    
    """Save the starlist to a text file"""
    
    fname = os.path.join(result_dir,basename + 'starlist.txt')
    params['starlist'] = fname
    fmt = '%s %10.6f %10.6f %7.2f %7.2f'
    aux = np.zeros(nstars,dtype=[('var1','U32'),('var2',float),\
        ('var3',float),('var4',float),('var5',float)])
    aux['var1'] = source_id
    aux['var2'] = ra
    aux['var3'] = dec
    aux['var4'] = x + 1.0
    aux['var5'] = y + 1.0
    np.savetxt(fname,aux.T,fmt=fmt)
    
    """Create dophot input files"""
    
    fname = os.path.join(result_dir,basename + 'dophot_input.txt')
    sigma = (params['fwhm']/2.35)**2
    idx = np.arange(1,nstars+1)
    fmt = '%%d 1 %5.1f %%7.2f %%7.2f 1000 %4.2f 0.00 %4.2f '%\
        (params['background'],sigma,sigma)
    fmt += '0 '*10
    np.savetxt(fname,np.array((idx,x,y)).T,fmt=fmt)
    
    """According to info from headers of science images, gain and 
    readnoise are the identical for all detectors.
    This is a bit suspicious
    """
    
    params['gain'] = 6.2
    params['rdnoise'] = 10.0
    
    create_dophot_parameter_input_file(basename,params,data_dir,result_dir)
    
    """Run DoPhot"""
    
    fname = os.path.join(result_dir,basename + 'dophot_params.txt')
    subprocess.run([f'{DOPHOT_DIR}/dophot',fname])
    
    fname = os.path.join(result_dir,basename + 'dophot_output.txt')
    params['output'] = fname

    #return source_id, ra, dec, params

    fname = os.path.join(result_dir,'params.json')
    with open(fname,'w') as fp:
        json.dump(params,fp)
