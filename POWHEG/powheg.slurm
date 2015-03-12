#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --mem-per-cpu=8G
#SBATCH --time=10:00:00
#SBATCH --account=cms_stage3
#SBATCH -o slurm_%A.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=yaxian.mao@vanderbilt.edu

export SCRAM_ARCH=slc6_amd64_gcc472

cd /home/maoy/working/GeneratorStudy/CMSSW_5_3_20/src
eval `scramv1 runtime -sh`

source /home/maoy/working/GeneratorStudy/CMSSW_5_3_20/src/setup.sh

cd /home/maoy/working/GeneratorStudy/CMSSW_5_3_20/src/POWHEG/POWHEG-BOX/Dijet/testrun-lhc

#name=${jobname}
#TheSeed=${seed}
echo "$TheSeed"
export Templatefile=powheg.input.bak
export Targetfile=powheg-new.input
export WorkDir=/scratch/maoy/${LS_JOBNAME}
mkdir -p ${WorkDir}
cd ${WorkDir}
cp /home/maoy/working/GeneratorStudy/CMSSW_5_3_20/src/POWHEG/POWHEG-BOX/Dijet/testrun-lhc/${Templatefile} ${Targetfile}
#cp powheg.input.bak ${Templatefile}
sed -i s/MYSEED/${TheSeed}/ ${Targetfile}
cp ${Targetfile} powheg.input

# Run the powheg generator
/home/maoy/working/GeneratorStudy/CMSSW_5_3_20/src/POWHEG/POWHEG-BOX/Dijet/pwhg_main > ${LS_JOBNAME}.out

mv *.out /scratch/maoy/logs
mv *.lhe ${LS_JOBNAME}.lhe
#sed '$d' < ${LS_JOBNAME}.lhe > temp_${LS_JOBNAME}.lhe
#mv temp_${LS_JOBNAME}.lhe ${LS_JOBNAME}.lhe
cp *.lhe /scratch/maoy/logs
cp ${LS_JOBNAME}.lhe /cms/store/user/ymao/pA5TEV/Generator/POWHEG/SevenTevMinkt20
cd /scratch/maoy
rm -rf ${WorkDir}
