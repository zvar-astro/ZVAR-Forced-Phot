import numpy as np
import os

from astropy.coordinates import SkyCoord, EarthLocation
from astropy.time import Time

PALOMAR = EarthLocation(lat='33:21:26.3',lon='-116:51:36',height=1712.)

MATCHFILE_DIR = "/data/matchfiles"

def JD_to_BJD (jd, ra, dec):
    """
    Transform JD to BJD
    
    Parameters
    ----------
    jd : array like
        Array of JDs
    ra : float
        Right ascention (deg)
    dec : float
        Declination (deg)
        
    Returns
    -------
    BJD - barycentric julian date
    """
    
    t = Time(jd,format='jd',scale='utc')
    t2 = t.tdb
    coord = SkyCoord(ra,dec,unit='deg')
    delta_t = t2.light_travel_time(coord,kind='barycentric',location=PALOMAR)
    BJD_TDB = t2 + delta_t
    return BJD_TDB.value

def save_to_file (epoch_image_info, params):

    n_imgs = len(epoch_image_info)

    if n_imgs == 0:
        # Generate the matchfile path+name
        fname = '%s/%04d/data_%04d_%02d_%1d_%2s.h5'%(MATCHFILE_DIR,params['field'],\
                     params['field'],params['ccd'],params['quad'],params['filter'])
        return fname

    jd = epoch_image_info['jd'].values
    bjd = JD_to_BJD(jd,params['ra_cen'],params['dec_cen'])
    epoch_image_info.insert(1,'bjd',bjd)
    epoch_image_info.sort_values(by='jd',inplace=True)
    field_names = list(epoch_image_info.columns)
    
    fname = params['res_dir'] + 'epochs.txt'
    with open(fname,'w') as fp:
        for i in range(n_imgs):
            row_data = [epoch_image_info[x].iloc[i] for x in field_names]
            fp.write('%14.6f %14.6f %1d %5.2f %1d %4d %2d %1d %1d %9d %6.3f %9.6f %12.10f %8.6f %9.6f %12.10f %5.2f %5.3f %6d %20s %s\n' %tuple(row_data))
    
    path = '%s/%04d/'%(MATCHFILE_DIR,params['field'])
    if not os.path.isdir(path):
        os.mkdir(path)

    # Generate the matchfile path+name
    fname = '%s/%04d/data_%04d_%02d_%1d_%2s.h5'%(MATCHFILE_DIR,params['field'],\
                 params['field'],params['ccd'],params['quad'],params['filter'])

    return fname

