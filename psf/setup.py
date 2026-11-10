from distutils.core import setup, Extension
import numpy as np
import sys

# python3 setup.py build_ext --inplace

module = Extension(
           name='PSF_Phot_PM', 
           sources=['psf_wrapper.c','psf.c'], 
           include_dirs=[np.get_include(),'/home/pmroz/software/cfitsio-3.49/include',],
           library_dirs=['/usr/local/lib','/home/pmroz/software/cfitsio-3.49/lib',],
           libraries=['cfitsio'],
           language='c')

setup(
    name = 'PSF_Phot_PM', 
    ext_modules = [module],
)
