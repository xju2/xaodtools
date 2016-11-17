#!/bin/bash

input_=$1
output_=$2

echo "---Start the Job---"

# setup root
#. /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/root/6.04.16-x86_64-slc6-gcc49-opt/bin/thisroot.sh
#. /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/x86_64/Gcc/gcc493_x86_64_slc6/slc6/gcc49/setup.sh
source /afs/cern.ch/user/x/xju/setup.sh

which root
which gcc
which python

script="/afs/cern.ch/user/x/xju/work/upsilon/code/MyXAODTools/scripts/read_monojet_minitree.py"

echo $script

python ${script} ${input_} smeared test.root -s

cp test.root ${output_}
