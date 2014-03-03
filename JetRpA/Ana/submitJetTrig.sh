
cd /afs/cern.ch/work/y/ymao/CMSSW_6_2_0_pre7/src 
#cmsenv
eval `scramv1 runtime -sh`

cd  /afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/NoUnfolding

export TRIG=$1

echo "Processing..."
echo "root -l -b -q anaJetTrackRpA.C+"

#Do the analysis
#root -b > runjob.log <<EOF
#.L anaInclusiveJS.C+
#anaInclusiveJS()
#.q
#EOF

root -b > runjob.log <<EOF
.L anaJetTrackRpA.C+
anaJetTrackRpA.C()
.q
EOF

echo "Done!"

#echo "Copying output files to " $destination
