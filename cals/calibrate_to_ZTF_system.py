import numpy as np
import argparse
import json
import matplotlib.pyplot as plt
import os

from penquins import Kowalski

DATABASE = 'ZTF_sources_20230821'
RESULTS_DIR = '/matchfiles/data/res/'
MATCH_RADIUS = 1.0

def get_ZTF_stars (kow, ra, dec, field, radius=0.01, filtr=2):
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
                        'field' : 1,
                        'filter' : 1,
                        'data.mag' : 1
                    }
                }
            }
        }
        }
    
    r = kow.query(query=q)
        
    if r['default']['status'] != 'success':
        print('Error executing the Kowalski query.')
        return np.array([]), np.array([]), np.array([])
    
    ra, dec, mag = [], [], []
    
    key = list(r['default'].get('data').get(DATABASE).keys())[0]

    for lc in r['default'].get('data').get(DATABASE).get(key):
        if lc['field'] != field:
            continue
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
parser.add_argument('field',type=int)
parser.add_argument('ccd',type=int)
parser.add_argument('quad',type=int)
parser.add_argument('filter',type=str)
parser.add_argument('--radius',default=0.1,type=float)
args = parser.parse_args()

if args.filter == 'g':
    filtr = 1
elif args.filter == 'r':
    filtr = 2
elif args.filter == 'i':
    filtr = 3
else:
    exit('Wrong filter name')

base_name = 'ztf_{:06d}_z{}_c{:02d}_q{:1d}'.format(args.field,
    args.filter,args.ccd,args.quad)
    
starlist = '{}{}/{}_starlist.txt'.format(RESULTS_DIR,base_name,base_name)
dophot_output = '{}{}/{}_dophot_output.txt'.format(RESULTS_DIR,
    base_name,base_name)
param_file = '{}{}/params.json'.format(RESULTS_DIR,base_name)

if not os.path.exists(starlist):
    exit('Cannot open file {}'.format(starlist))
if not os.path.exists(dophot_output):
    exit('Cannot open file {}'.format(dophot_output))
if not os.path.exists(param_file):
    exit('Cannot open file {}'.format(param_file))

#-----Read data from text files----#

ra, dec = np.loadtxt(starlist,unpack=True,usecols=(1,2),dtype=float)
typ, mag_raw = np.loadtxt(dophot_output,unpack=True,usecols=(1,4),
    dtype='int,float')
with open(param_file,'r') as f:
    params = json.load(f)

N = len(ra)
typ = typ[:N]
mag_raw = mag_raw[:N] + params['magzp']

#-----Initialize the Kowalski connection-----#

with open('passwd.json','r') as f:
    passwd = json.load(f)

protocol, host, port = 'https', 'gloria.caltech.edu', 443

k = Kowalski(
    token=passwd['gloria']['token'],
    protocol=protocol, 
    host=host, 
    port=port, 
    timeout=240.0
)

if not k.ping(): exit('Error while connecting to Kowalski.')

#-----Get median magnitudes from Kowalski-----#

ra_ztf, dec_ztf, mag_ztf = get_ZTF_stars(k,params['ra_cen'],
    params['dec_cen'],args.field,filtr=filtr,radius=args.radius)

#-----Find the zero-point correction-----#

zp, n_match = match_stars(ra,dec,mag_raw,typ,ra_ztf,dec_ztf,mag_ztf)

#print('Zero-point correction = {:.3f}'.format(zp))
#print('Number of common stars = {}'.format(n_match))
print('{:3d} {:2d} {:1d} {} {:6.3f} {:5d}'.format(args.field,args.ccd,args.quad,args.filter,zp,n_match))
