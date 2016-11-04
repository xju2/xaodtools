#!/bin/bash

input_=$1
output_=$2

#. /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.16-x86_64-slc6-gcc49-opt/bin/thisroot.sh
#. /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc493_x86_64_slc6/slc6/gcc49/setup.sh
#source /afs/cern.ch/user/x/xju/work/h4l/h4lcode/workspaces/setup.sh
source /afs/cern.ch/user/x/xju/setup.sh

which root
which python
script="/afs/cern.ch/user/x/xju/work/upsilon/code/MyXAODTools/scripts/draw_upsilon.py"
#python /afs/cern.ch/user/x/xju/work/upsilon/run/data12_v3/new_upsilon.py  make $input_ $output_

new_input=`echo ${input_} | awk -F/ '{printf("%s/%s/%s",$9,$10,$11)}'`
new_input=${GROUPEOSDIR}bphys/merged/$new_input
echo "${new_input}"

#python $script make $new_input $output_ --oldPtCut
python $script make $new_input $output_
