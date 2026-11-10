import argparse
import h5py
import os
import numpy as np

"""
This program pulls out the light curve of a given object. The program
prints BJD, differential flux and its uncertainty, measurement flag,
and program ID to the standard output.

Usage:
python3 getobj.py MATCHFILE OBJECT_ID

P. Mroz @ Caltech, 26 Jan 2021
"""

def getobj (objid, matchfile):
    if not os.path.exists(matchfile): 
        tmp = np.empty([])
        return tmp, tmp, tmp, tmp, tmp
    with h5py.File(matchfile,'r') as fp:
        source_id = fp['data']['sources']['gaia_id']
        _idx = np.where(source_id == objid)
        if len(_idx) == 0: 
            tmp = np.empty([])
            return tmp, tmp, tmp, tmp, tmp
        bjd = fp['data']['exposures']['bjd']
        pid = fp['data']['exposures']['pid']
        n_exp = len(bjd)
        idx = _idx[0][0]
        tmp = idx*n_exp + np.arange(n_exp)
        flux = fp['data']['sourcedata']['flux'][tmp]
        ferr = fp['data']['sourcedata']['flux_err'][tmp]
        flag = fp['data']['sourcedata']['flag'][tmp]
        return bjd, flux, ferr, flag, pid


parser = argparse.ArgumentParser()
parser.add_argument('matchfile',type=str,help='Matchfile name')
parser.add_argument('objid',type=int,help='PanSTARRS objects ID')
args = parser.parse_args()

bjd, flux, ferr, flag, pid = getobj(args.objid,args.matchfile)

npts = len(bjd)

for i in range(npts):
    if np.isnan(flux[i]) or np.isnan(ferr[i]): continue
    print('%13.5f %10.3f %10.3f %1d %1d'%(bjd[i],flux[i],ferr[i],flag[i],pid[i]))

