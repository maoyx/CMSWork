#!/bin/bash
#$1 number of jobs
cd /home/maoy/working/GeneratorStudy/CMSSW_5_3_20/src
eval `scramv1 runtime -sh`

source /home/maoy/working/GeneratorStudy/CMSSW_5_3_20/src/setup.sh

cd /home/maoy/working/GeneratorStudy/CMSSW_5_3_20/src/POWHEG/POWHEG-BOX/Dijet/testrun-lhc


pwd=`pwd`

nJobs=399
i=0
while [ $i -le $nJobs ];
do
  jobname=powheg_Job_${i}
  seed=1000${i}
  export LS_JOBNAME=${jobname}
  export TheSeed=${seed}
  echo "$TheSeed"
   SUBID=$(sbatch /home/maoy/working/GeneratorStudy/CMSSW_5_3_20/src/powheg.slurm ${TheSeed},${LS_JOBNAME})
   echo "$SUBID seed=$TheSeed" >> batchSubmission.$(date +%s).log
  let "i++"

done

echo "done for powheg lhe generation submission!"

