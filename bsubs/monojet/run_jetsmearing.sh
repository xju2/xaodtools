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

script="/afs/cern.ch/user/x/xju/work/upsilon/code/MyXAODTools/scripts/monojet/jetsmearing.py"
echo $script
echo $
#new_input=`echo ${input_} | awk -F/ '{printf("%s/%s/%s",$9,$10,$11)}'`
#new_input=${GROUPEOSDIR}bphys/merged/$new_input
#echo "${new_input}"

#python $script make $new_input $output_ --oldPtCut
python $script $input_ $output_
