#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=4G
#SBATCH --mem-per-cpu=4G
#SBATCH --time=1:00:00
#SBATCH --account=cms_stage3
#SBATCH -o job_rivet_%A.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=yaxian.mao@vanderbilt.edu

export SCRAM_ARCH=slc6_amd64_gcc472

cd /home/maoy/working/GeneratorStudy/RIVET/CMSSW_5_3_20/src
eval `scramv1 runtime -sh`

cd /home/maoy/working/GeneratorStudy/RIVET/CMSSW_5_3_20/src


cmsRun rivet_cfg.py inputFiles="$INPUT" output="$OUTPUT" sqrtS="$Z2ENERGY" &> run-${JOBNAME}.log

echo "cmsRun is finished, now converting the output to root"
aida2root ${OUTPUT}


echo  "converting to root is done!"

