#!/bin/bash

#setup ROOTCORE
input=$1
output=$2

echo "---Start the Job---"
echo ${ROOTCOREBIN}
echo ${ROOTCOREDIR}
if [ "x${ROOTCOREBIN}" == "x" ]; then
    shift $#
    source /afs/cern.ch/user/x/xju/work/upsilon/code/rcSetup.sh
else
. /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc493_x86_64_slc6/slc6/gcc49/setup.sh /afs/cern.ch/sw/lcg/contrib
fi

which root
which gcc

echo "${ROOTCOREBIN}/bin/x86_64-slc6-gcc49-opt/analysis $input ${output}"

#analysis upsilon $input -1 useBphy1=1 >& log
analysis upsilon $input -1 useBphy1=0 >& log

#mv reduced_ntup.root $output
scp reduced_ntup.root $output
scp log ${output}.log
