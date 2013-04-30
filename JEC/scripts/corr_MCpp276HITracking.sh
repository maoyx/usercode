root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt15/HiForest_v81_merged01/pt15_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppHITrkJet15.root",30,2.034E-01,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt30/HiForest_v81_merged01/pt30_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppHITrkJet30.root",50,1.075E-02,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt50/HiForest_v81_merged01/pt50_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppHITrkJet50.root",80,1.025E-03,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt80/HiForest_v81_merged01/pt80_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppHITrkJet80.root",120, 9.865E-05,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt120/HiForest_v81_merged01/pt120_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppHITrkJet120.root",170,1.129E-05,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt170/HiForest_v81_merged01/pt170_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppHITrkJet170.root",220,1.465E-06,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt220/HiForest_v81_merged01/pt220_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppHITrkJet220.root",280,2.837E-07,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt280/HiForest_v81_merged01/pt280_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppHITrkJet280.root",9999,5.323E-08,true)
.q
EOF

rm -f /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/jet_Corrpp276HITrack_mergedpthatbins_2013MC.root
hadd /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/jet_Corrpp276HITrack_mergedpthatbins_2013MC.root /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppHITrkJet*.root


