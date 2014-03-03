
cd /afs/cern.ch/work/y/ymao/CMSSW_6_2_0_pre7/src 
#cmsenv
eval `scramv1 runtime -sh`

cd  /afs/cern.ch/work/y/ymao/analysis/AsymmetryPA/NoUnfolding

export PTHAT=$1
export PTMAX=$2

echo "Processing..."
echo "root -l -b -q anaTrackAsy.C+"

#Do the analysis
#root -b > runjob.log <<EOF
#.L anaInclusiveJS.C+
#anaInclusiveJS()
#.q
#EOF

root -b > runjob.log <<EOF
.L anaTrackAsy.C+
anaTrackAsy()
.q
EOF

echo "Done!"

#echo "Copying output files to " $destination
