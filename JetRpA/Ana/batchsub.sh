
cd /afs/cern.ch/work/y/ymao/CMSSW_6_2_0_pre7/src 
#cmsenv
eval `scramv1 runtime -sh`

cd  /afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/NoUnfolding

HLTTRIG=("Jet100", "Jet80")
for ((i=0; i < ${#HLTTRIG[@]}; i++))
     do
      export TRIG=${HLTTRIG[$i]}
 #  export TRIG=$1
   echo " trigName = $TRIG "
   echo "Processing..."

  echo "root -l -b -q anaJetTrackRpA.C+"


root -b > runjob.log <<EOF
.L anaJetTrackRpA.C+
anaJetTrackRpA()
.q
EOF

      echo "Done for TrigName = $TRIG"
     done
     
echo "Done all jobs!"

#echo "Copying output files to " $destination
