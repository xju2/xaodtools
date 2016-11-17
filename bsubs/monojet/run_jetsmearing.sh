#!/bin/bash

input_=$1
output_=$2

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

script="/afs/cern.ch/user/x/xju/work/upsilon/code/MyXAODTools/scripts/jetsmearing.py"

python $script $input_ test.root --get_hist
cp test.root $output_
