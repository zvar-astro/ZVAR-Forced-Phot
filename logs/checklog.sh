#!/bin/bash


# Set base directory for log files
base_dir="/data/ZVAR-Forced-Phot/logs"


# Remove existing failed_quads file if it exists
if [ -f $base_dir/failed_quads.txt ] 
then
    rm $base_dir/failed_quads.txt
fi

# Set field range for logs to check
field_start=245
field_end=255


# Iterate over fields
for field in $(eval echo {$field_start..$field_end})
do

  log_file=$(printf "%s/log_files/pipeline_%04d.log" "$base_dir" "$field")
  python checklog.py --logfile=$log_file 

done
