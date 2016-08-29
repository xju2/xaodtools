#!/bin/bash

#update libararies
if [ "x${ROOTCOREBIN}" == "x" ]; then
  echo "ROOTCOREDIR is not defined!"
  exit
fi
griddir=`pwd`
##copy RootCore
mkdir RootCore
mkdir RootCore/lib    
mkdir RootCore/data  
mkdir RootCore/python 
mkdir RootCore/include
cd ${ROOTCOREBIN}/lib/x86_64-slc6-gcc49-opt; cp -H * ${griddir}/RootCore/lib/
cd ${ROOTCOREBIN}/data; cp -rH * ${griddir}/RootCore/data/
cd ${ROOTCOREBIN}/python; cp -rH * ${griddir}/RootCore/python/
cd ${ROOTCOREBIN}/include; cp -rH * ${griddir}/RootCore/include/

cd ${griddir}
###executable directory
mkdir bin; 

cp ${ROOTCOREBIN}/bin/x86_64-slc6-gcc49-opt/upsilon bin/
cp ${ROOTCOREBIN}/bin/x86_64-slc6-gcc49-opt/analysis bin/

tar -czf hello.tar bin/ RootCore/ 
rm -rf  bin/ RootCore/ 
