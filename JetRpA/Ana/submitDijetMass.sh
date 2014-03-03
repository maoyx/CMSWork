
cd /afs/cern.ch/work/y/ymao/CMSSW_6_2_0_pre7/src 
#cmsenv
eval `scramv1 runtime -sh`

cd  /afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/NoUnfolding


   export TRIG="Jet100"
   echo " trigName = $TRIG "
   echo "Processing..."

  echo "root -l -b -q anaDijetMassRpA.C+"


root -b > runjob${TRIG}.log <<EOF
.L anaDijetMassRpA.C+
anaDijetMassRpA()
.q
EOF

 echo "Done for TrigName = $TRIG"


echo "Done all jobs!"

#echo "Copying output files to " $destination
