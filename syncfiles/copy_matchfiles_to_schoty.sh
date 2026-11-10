#!/bin/bash

# This script uses rsync over ssh, which allows for
# automatic authentication using a public ssh key that
# has already been added to ~/.ssh/authorized_keys on
# the destination machine (schoty in this case). You
# can follow steps here to add an ssh key to schoty
# for your own account:
#
# https://www.digitalocean.com/community/tutorials/how-to-copy-files-with-rsync-over-ssh
#
# Author: Zach Vanderbosch (Caltech)

username="zvanderbosch"

for field in {245..1876}
do
    
    localdir=$(printf "/data/matchfiles/%04d" "$field")
    remotedir=$(printf "%s@schoty.caltech.edu:/gdata/Data/Matchfiles/matchfiles/%04d/" "$username" "$field") 
    if [ -d $localdir ]; then
        echo $field
        rsync -avh --update --log-file=rsync.log -e ssh $localdir/*.h5 $remotedir
    fi
done
