#!/bin/bash
#$1 number of jobs
cd /home/maoy/working/GeneratorStudy/RIVET/CMSSW_5_3_20/src
eval `scramv1 runtime -sh`


cd /home/maoy/working/GeneratorStudy/RIVET/CMSSW_5_3_20/src


   for file in `cat filelistAll.dat`; do
#   for file in `cat try.dat`; do
    dataset+=("${file[@]}")
   done

export Z2ENERGY=7000.

for ((i=0; i<${#dataset[@]}; i++)); do
  export INPUT="${dataset[${i}]}"
  export OUTPUT="TuneZ2_POWHEG_CT10NLO_7TeVR5_supfactor250_ktmin5_file${i}.aida"
  export JOBNAME=RIVET_Job_${i}
  echo "input = ${INPUT} with file ${i}"
   SUBID=$(sbatch rivet.slurm ${INPUT},${OUTPUT},${Z2ENERGY},${JOBNAME})
   echo "$SUBID input=$INPUT" >> batchSubmission.$(date +%s).log
done

echo "submit all the RIVET analysis for LHE files !"

