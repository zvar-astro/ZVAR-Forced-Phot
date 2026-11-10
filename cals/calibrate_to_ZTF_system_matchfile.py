import numpy as np
import argparse
import json
import matplotlib.pyplot as plt
import os
import h5py

from penquins import Kowalski

DATABASE = 'ZTF_sources_20201201'
MATCH_RADIUS = 1.0

def get_ZTF_stars (kow, ra, dec, radius=0.01, filtr=2):
    """
    Query the Kowalski DB to get a list of ZTF sources.
    
    Parameters
    ----------
    kow : Kowalski object
        a Kowalski instance
    ra : float
        Right ascention (deg)
    dec : float
        Declination (deg)
    radius : float, optional
        Search radius (arcsec), default: 0.01 deg
    
    Returns
    -------
    ra, dec, mag : numpy arrays containing RA, Dec, and median magnitude
    """    
        
    q = {'query_type': 'cone_search',
        'query': {
            'object_coordinates': {
               'radec': '[(%.6f,%.6f)]'%(ra,dec), 
                'cone_search_radius': '%.4f'%radius,
                'cone_search_unit': 'deg'
            },
            'catalogs': {
                DATABASE: {
                    'filter': {},
                    'projection': {
                        'ra' : 1,
                        'dec' : 1,
                        'filter' : 1,
                        'data.mag' : 1
                    }
                }
            }
        }
        }
    
    r = kow.query(query=q)
        
    if r['status'] != 'success':
        print('Error executing the Kowalski query.')
        return np.array([]), np.array([]), np.array([])
    
    ra, dec, mag = [], [], []
    
    key = list(r.get('data').get(DATABASE).keys())[0]

    for lc in r.get('data').get(DATABASE).get(key):
        if lc['filter'] != filtr:
            continue
        ra.append(lc['ra'])
        dec.append(lc['dec'])
        x = np.array([x['mag'] for x in lc['data']])
        mag.append(np.median(x))
    
    return np.array(ra), np.array(dec), np.array(mag)
    
def match_stars (ra_dophot, dec_dophot, mag_dophot, type_dophot, 
    ra_ztf, dec_ztf, mag_ztf):
        
    N_stars = len(ra_ztf)
    mag_diff = []

    for i in range(N_stars):
        adist = ((ra_dophot-ra_ztf[i])*np.cos(dec_ztf[i]*np.pi/180.0))**2
        adist += (dec_dophot-dec_ztf[i])**2
        adist = np.sqrt(adist)
        idx = np.argmin(adist)
        if adist[idx]*3600.0 <= MATCH_RADIUS and type_dophot[idx] == 1:
            mag_diff.append(mag_ztf[i]-mag_dophot[idx])

    zp_corr = np.median(mag_diff)
    N = len(mag_diff)

    #plt.hist(mag_diff,range=(-0.5,0.5),bins=100)
    #plt.show()

    return zp_corr, N


parser = argparse.ArgumentParser()
parser.add_argument('matchfile',type=str)
parser.add_argument('filter',type=str)
parser.add_argument('ra',type=float)
parser.add_argument('dec',type=float)
parser.add_argument('--radius',default=0.01,type=float)
args = parser.parse_args()

if args.filter == 'g':
    filtr = 1
elif args.filter == 'r':
    filtr = 2
elif args.filter == 'i':
    filtr = 3
else:
    exit('Wrong filter name')

with h5py.File(args.matchfile,'r') as f:
    sources = f['data']['sources']
    ra = sources['ra']
    dec = sources['decl']
    mag_raw = sources['mag_ref']
    typ = sources['objtype']
    mag_zp = f['data'].attrs['magzp']
    mag_raw += mag_zp
    
#-----Initialize the Kowalski connection-----#

with open('passwd.json','r') as f:
    passwd = json.load(f)

protocol, host, port = 'https', 'kowalski.caltech.edu', 443

k = Kowalski(username=passwd['kowalski']['username'],
    password=passwd['kowalski']['password'],
    protocol=protocol, host=host, port=port)

if not k.ping(): exit('Error while connecting to Kowalski.')

#-----Get median magnitudes from Kowalski-----#

ra_ztf, dec_ztf, mag_ztf = get_ZTF_stars(k,args.ra,args.dec,filtr=filtr,
    radius=args.radius)

#-----Find the zero-point correction-----#

zp, n_match = match_stars(ra,dec,mag_raw,typ,ra_ztf,dec_ztf,mag_ztf)

print('Zero-point correction = {:.3f}'.format(zp))
print('Number of common stars = {}'.format(n_match))
