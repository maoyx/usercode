
export SCRAM_ARCH=slc5_amd64_gcc434
source /osg/app/cmssoft/cms/cmsset_default.sh

cd /net/hidsk0001/d00/scratch/maoyx/JetShape/CMSSW_4_4_2_patch5/src
#cmsenv
eval `scramv1 runtime -sh`

#cd /net/hidsk0001/d00/scratch/maoyx/JetShape/CMSSW_4_4_2_patch5/src/HiForest_2013
cd /net/hidsk0001/d00/scratch/maoyx/JetShape/CMSSW_4_4_2_patch5/src/HiForest_V2_02_16

export PTHAT=$1
export PTMAX=$2

echo "Processing..."
echo "root -l -b -q anaLeadingJS.C+"

#Do the analysis
root -b > runleading.log <<EOF
.L anaLeadingJS.C+
anaLeadingJS()
.q
EOF


echo "Done!"

#echo "Copying output files to " $destination
