#!/bin/sh
root -l -b<<EOF
gSystem->Load("../HiForest_V3/hiForest_h.so")
.x IndResponse.C++(15,"pp", 0);
.q
EOF

root -l -b<<EOF
gSystem->Load("../HiForest_V3/hiForest_h.so")
.x IndResponse.C++(30,"pp", 0);
.q
EOF

root -l -b<<EOF
gSystem->Load("../HiForest_V3/hiForest_h.so")
.x IndResponse.C++(50,"pp", 0);
.q
EOF

root -l -b<<EOF
gSystem->Load("../HiForest_V3/hiForest_h.so")
.x IndResponse.C++(80,"pp", 0);
.q
EOF

root -l -b<<EOF
gSystem->Load("../HiForest_V3/hiForest_h.so")
.x IndResponse.C++(120,"pp", 0);
.q
EOF

root -l -b<<EOF
gSystem->Load("../HiForest_V3/hiForest_h.so")
.x IndResponse.C++(170, "pp", 0);
.q
EOF

root -l -b<<EOF
gSystem->Load("../HiForest_V3/hiForest_h.so")
.x IndResponse.C++(220,"pp", 0);
.q
EOF

root -l -b<<EOF
gSystem->Load("../HiForest_V3/hiForest_h.so")
.x IndResponse.C++(280,"pp",0);
.q
EOF
root -l -b<<EOF
gSystem->Load("../HiForest_V3/hiForest_h.so")
.x IndResponse.C++(370,"pp",0);
.q
EOF
rm -f /net/hidsk0001/d00/scratch/maoyx/pPb/ClosureJEC/CMSSW_5_3_8_HI_patch2/src/JECPawan/Output/pPbResponse_PYTHIA_merged.root
hadd /net/hidsk0001/d00/scratch/maoyx/pPb/ClosureJEC/CMSSW_5_3_8_HI_patch2/src/JECPawan/Output/pPbResponse_PYTHIA_merged.root /net/hidsk0001/d00/scratch/maoyx/pPb/ClosureJEC/CMSSW_5_3_8_HI_patch2/src/JECPawan/Output/Response_pp_pT*.root
