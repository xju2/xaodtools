#!/bin/bash

dir_=$1
file_=$2

source /afs/cern.ch/atlas/software/dist/AtlasSetup/scripts/asetup.sh AtlasProduction,17.2.13.7,here

if [ -f ${dir_}/$file_ ]; then

    cp ${dir_}/$file_ .
    out_name="merged_${file_}.root"

    acmd.py merge-files `cat ${file_}` -o $out_name

    scp ${out_name} $dir_/
else
    echo "no input $file_"
fi
