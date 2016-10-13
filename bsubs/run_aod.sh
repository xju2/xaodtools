#!/bin/bash

dir_=$1
file_=$2

current_dir=$PWD

test_area="/afs/cern.ch/user/x/xju/work/upsilon/athena/"
cd  $test_area
source /afs/cern.ch/atlas/software/dist/AtlasSetup/scripts/asetup.sh AtlasProduction,17.2.13.7,here

echo $current_dir
cd $current_dir

ls

jo="/afs/cern.ch/user/x/xju/work/upsilon/athena/PhysicsAnalysis/AnalysisCommon/AnalysisExamples/share/VFitZmmOnAOD_jobOptions_AutoConfig_template.py"
if [ -f ${dir_}/$file_ ]; then

    cp $jo .
    cp ${dir_}/$file_ .

    file_name=""
    for line in `cat $file_`
    do
        file_name="'$line',$file_name"
    done
    echo "jp.AthenaCommonFlags.FilesInput = [$file_name]" > input.py
    echo "my input: `cat input.py`"

    athena.py VFitZmmOnAOD_jobOptions_AutoConfig_template.py

    out_name="VFitZmmOnAOD_${file_}.root"
    mv VFitZmmOnAOD.root $out_name
    scp  $out_name $dir_/
else
    echo "no input $file_"
fi
