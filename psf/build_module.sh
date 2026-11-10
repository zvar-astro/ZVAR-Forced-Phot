#!/bin/bash

# Simple script to build and rename the PSM_Phot_PM module
python3 setup.py build_ext --inplace
mv PSF_Phot_PM.cpython-39-x86_64-linux-gnu.so PSF_Phot_PM.so
