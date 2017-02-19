#!/bin/bash

input_file=$1
out_name=$2
Seed=$3
firstEvt=$4

# copy to local machine
cp $input_file .
cp /afs/cern.ch/user/x/xju/work/h4l/highmass/showering_interference/MC15.344784.gg2vvPy8EG_ggH300_5SMW_ZZ_4l.py . 
file_name=`basename ${input_file}`
JO_NAME="MC15.344784.gg2vvPy8EG_ggH300_5SMW_ZZ_4l.py"

# count lines and decided how many jobs needed
n_lines=`wc -l ${file_name} | awk '{print $1}'`
n_jobs=`expr ${n_lines} / 11 / 5000`
echo "total jobs: ${n_jobs}"

# change input file
LHE_Name="group.phys-gener.gg2VV.344784.ggH300_5SMW_ZZ_4l_13TeV.TXT.mc15_v1._00001.events"
tar -cf $LHE_Name $file_name

# setup athena and run

  export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
  source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
  export JOBOPTSEARCHPATH=/cvmfs/atlas.cern.ch/repo/sw/Generators/MC15JobOptions/latest/common:$JOBOPTSEARCHPATH

  for (( ijob=0; ijob<${n_jobs}; ijob++ ))
  do
    echo "running job ${ijob}"
    asetup --unset
    asetup 19.2.5.14,AtlasProduction,here
    rm test.evgen.pool.root
    new_event=`expr ${firstEvt} + 5000 \* ${ijob}`
    new_seed=`expr ${Seed} + ${ijob}`
    Generate_tf.py --ecmEnergy=13000. --maxEvents=5000 --runNumber=344784 --firstEvent=${new_event} --randomSeed=${new_seed} --outputEVNTFile=test.evgen.pool.root --jobConfig=${JO_NAME} --inputGeneratorFile=$LHE_Name  >& log.gen

    ## run the truth xAOD
    if [ ! -d truth ];then
        mkdir truth
    fi
    cd truth
    asetup --unset
    asetup 20.1.5.7,AtlasDerivation,gcc48,here
    Reco_tf.py --inputEVNTFile ../test.evgen.pool.root --outputDAODFile test.TruthDAOD.root --reductionConf TRUTH1  >& log.truth

    ## copy files
    #cp DAOD_TRUTH4.test.TruthDAOD_2.root ${out_name//.root/_2.root}
    new_out_name=`echo ${out_name} | sed "s/.root/_${ijob}.root/"`
    xrdcp DAOD_TRUTH4.test.TruthDAOD.root root://eosatlas//eos/atlas/unpledged/group-wisc/users/xju/H4l/showering_interference_truth1/${new_out_name}
    rm -f *

    cd ../
  done
