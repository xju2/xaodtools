#!/bin/bash

dir_=$1
file_=$2

. /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.16-x86_64-slc6-gcc49-opt/bin/thisroot.sh
. /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc493_x86_64_slc6/slc6/gcc49/setup.sh

cp ${dir_}/$file_ .
out_name="merged_${file_}.root"
hadd ${out_name} `cat ${file_}`

#scp ${out_name} $dir_/
out_dir=`echo ${dir_} | awk -F/ '{printf("%s/%s",$9,$10)}'`
echo "${GROUPEOSDIR}bphys/merged/$out_dir/${out_name}"

xrdcp -f ${out_name} ${GROUPEOSDIR}bphys/merged/$out_dir/${out_name}
