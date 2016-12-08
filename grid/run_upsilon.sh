#!/bin/sh
is_data=$1

tar -xf hello.tar
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:RootCore/lib:

echo "------------"
echo "|          |"
echo "---start----"
echo "------------"
#echo $LD_LIBRARY_PATH
which  gcc
which root

export ROOTCOREDIR="RootCore/"
export ROOTCOREBIN="RootCore/"
export CALIBPATH="/cvmfs/atlas.cern.ch/repo/sw/database/GroupData:${ROOTCOREBIN}data:http//atlas.web.cern.ch/Atlas/GROUPS/DATABASE/GroupData"
export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:"$ROOTCOREBIN/include"
echo $CALIBPATH
#cat toberun.txt
cat toberun.txt | sed 's/,/\n/g' > jj.txt
mv jj.txt toberun.txt
cat toberun.txt

#isData=0 isAtlfast=0 noSys=1 doSmear=0 doPhoton=0 debug=0
#./bin/analysis upsilon toberun.txt -1 isData=${is_data} useBphy1=1
./bin/analysis upsilon toberun.txt -1 isData=${is_data} useBphy1=0

echo "---------------"
echo "The job is done"
