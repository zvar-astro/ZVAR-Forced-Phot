import argparse
import json

from penquins import Kowalski
from utils import *

"""
This programs returns the field, ccd, and quadrant number of a given 
point, as well as the PanSTARRS DR1 identifier of the star. 

These calculations are approximate and may fail if the object is located 
close to the edge of the field of view / edge of the reference image.

The program queries Kowalski and it attempts to read a JSON file  
passwd.json with your Kowalski credentials:
cat passwd.json
{"kowalski": {"username": "YOUR_KOWALSKI_USERNAME", 
"password": "YOUR_KOWALSKI_PASSWORD"}}

Usage: radec2field.py [--radius RADIUS] ra dec

P. Mroz @ Caltech, 26 Jan 2021
"""

parser = argparse.ArgumentParser()
parser.add_argument('ra',type=float,help='Right ascension (deg)')
parser.add_argument('dec',type=float,help='Declination (deg)')
parser.add_argument('--radius',default=2.0,type=float,
    help='Search radius (arcsec), default: 2.0 arcsec')
args = parser.parse_args()

#-----Initialize the Kowalski connection-----#

with open('passwd.json','r') as f:
    passwd = json.load(f)

protocol, host, port = 'https', 'kowalski.caltech.edu', 443

k = Kowalski(username=passwd['kowalski']['username'],
    password=passwd['kowalski']['password'],
    protocol=protocol, host=host, port=port)

if not k.ping(): exit('Error while connecting to Kowalski.')

#-----Read the field list and find the field/ccd/quadrant-----#

fieldno, ra, dec = read_fields()

res = get_field_id(args.ra,args.dec,fieldno,ra,dec)

#-----Get PS1 identifiers-----#

ident_ps, ra_ps, dec_ps, adist_ps = get_PS1_stars(k,args.ra,args.dec,radius=args.radius)

print('# matchfile PS1_ID field ccd quad ang_dist')
for x in res:
    field, ccd, quad = x
    for i,ident in enumerate(ident_ps):
        fname = 'matchfiles/%04d/data_%04d_%02d_%1d_zr.h5'%(field,field,ccd,quad)
        print('%s %d %3d %2d %1d %5.2f'%(fname,ident,field,ccd,quad,adist_ps[i]))
