#!/bin/bash


# Set field range
field_start=1241
field_end=1897


# Set base directory for pipeline routines
base_dir="/data/ZVAR-Forced-Phot"


# Iterate over fields/quadrants
for field in $(eval echo {$field_start..$field_end})
do

for ccd in {1..16}
do

    for quad in {1..4}
    do
   
    resdirg=$(printf "/matchfiles/data/res/ztf_%06d_zg_c%02d_q%1d" "$field" "$ccd" "$quad")
    rm -rf /home/zvanderbosch/.ztf_cookie
    echo $field $ccd $quad g 
    /home/pmroz/software/env/bin/python3 $base_dir/psf-phot-pm.py $field $ccd $quad g
    rm -rf /data/data/*
    rm -rf $resdirg/*_out.fits
    rm -rf $resdirg/shadowfile_out
    rm -rf $resdirg/log.txt
    rm -rf $resdirg/error_out
    
    resdirr=$(printf "/matchfiles/data/res/ztf_%06d_zr_c%02d_q%1d" "$field" "$ccd" "$quad")
    rm -rf /home/zvanderbosch/.ztf_cookie
    echo $field $ccd $quad r
    /home/pmroz/software/env/bin/python3 $base_dir/psf-phot-pm.py $field $ccd $quad r
    rm -rf /data/data/*
    rm -rf $resdirr/*_out.fits
    rm -rf $resdirr/shadowfile_out
    rm -rf $resdirr/log.txt
    rm -rf $resdirr/error_out
    
    resdiri=$(printf "/matchfiles/data/res/ztf_%06d_zi_c%02d_q%1d" "$field" "$ccd" "$quad")
    rm -rf /home/zvanderbosch/.ztf_cookie
    echo $field $ccd $quad i
    /home/pmroz/software/env/bin/python3 $base_dir/psf-phot-pm.py $field $ccd $quad i
    rm -rf /data/data/*
    rm -rf $resdiri/*_out.fits
    rm -rf $resdiri/shadowfile_out
    rm -rf $resdiri/log.txt
    rm -rf $resdiri/error_out
    
    done

done

done
