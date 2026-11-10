#!/bin/bash


# Set base directory for pipeline routines
base_dir="/data/ZVAR-Forced-Phot"


# Set filenames for input quadrant list and output log
logfile=$(printf "%s/logs/log_files/pipeline_rerun.log" "$base_dir")
input=$(printf "%s/logs/failed_quads.txt" "$base_dir")


# Remove previous log file if it exists
if [ -f $logfile ]
then
    rm $logfile
fi


# Iterate over failed quadrants
while read -r line
do
    infoarray=($line)
    resdir=$(printf "/matchfiles/data/res/ztf_%06d_z%s_c%02d_q%1d" "${infoarray[0]}" "${infoarray[3]}" "${infoarray[1]}" "${infoarray[2]}")
    rm -rf /home/zvanderbosch/.ztf_cookie
    echo ${infoarray[0]} ${infoarray[1]} ${infoarray[2]} ${infoarray[3]}
    /home/pmroz/software/env/bin/python3 psf-phot-pm.py --logfile $logfile ${infoarray[0]} ${infoarray[1]} ${infoarray[2]} ${infoarray[3]}
    rm -rf /data/data/*
    rm -rf $resdir/*_out.fits
    rm -rf $resdir/shadowfile_out
    rm -rf $resdir/log.txt
    rm -rf $resdir/error_out
  
done < "$input"
