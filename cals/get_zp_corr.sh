#!/bin/bash

# Use {245..881} for primary grid or {1241..1896} for secondary grid
for field in {245..881}
do

for ccd in {1..16}
do

    for quad in {1..4}
    do

    echo $field $ccd $quad

    gval=$(python /data/ZVAR-Forced-Phot/cals/calibrate_to_ZTF_system.py --radius=0.1 $field $ccd $quad g)
    rval=$(python /data/ZVAR-Forced-Phot/cals/calibrate_to_ZTF_system.py --radius=0.1 $field $ccd $quad r)
    ival=$(python /data/ZVAR-Forces-Phot/cals/calibrate_to_ZTF_system.py --radius=0.1 $field $ccd $quad i)

    garray=($gval)
    rarray=($rval)
    iarray=($ival)
    gmagcorr=${garray[4]}
    rmagcorr=${rarray[4]}
    imagcorr=${iarray[4]}
    gstars=${garray[5]}
    rstars=${rarray[5]}
    istars=${iarray[5]} 

    if [ $field -lt 1000 ]
    then
        printf "%04d %02d %1d %6.3f %5d %6.3f %5d %6.3f %5d\n" "$field" "$ccd" "$quad" "$gmagcorr" "$gstars" "$rmagcorr" "$rstars" "$imagcorr" "$istars" >> /data/ZVAR-Forced-Phot/cals/zero_point_corrections_primary.txt
    else
        printf "%04d %02d %1d %6.3f %5d %6.3f %5d %6.3f %5d\n" "$field" "$ccd" "$quad" "$gmagcorr" "$gstars" "$rmagcorr" "$rstars" "$imagcorr" "$istars" >> /data/ZVAR-Forced-Phot/cals/zero_point_corrections_secondary.txt
    fi

    done

done

done
