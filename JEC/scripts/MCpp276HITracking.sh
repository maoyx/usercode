root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt15/HiForest_v81_merged01/pt15_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet15.root",30,2.034E-01,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt30/HiForest_v81_merged01/pt30_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet30.root",50,1.075E-02,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt50/HiForest_v81_merged01/pt50_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet50.root",80,1.025E-03,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt80/HiForest_v81_merged01/pt80_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet80.root",120, 9.865E-05,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt120/HiForest_v81_merged01/pt120_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet120.root",170,1.129E-05,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt170/HiForest_v81_merged01/pt170_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet170.root",220,1.465E-06,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt220/HiForest_v81_merged01/pt220_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet220.root",280,2.837E-07,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_AllJet.C+("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod22/JEC_Pythia_pt280/HiForest_v81_merged01/pt280_JEC_ppHiIterativeTrack_P01_prod22_v81_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet280.root",9999,5.323E-08,false)
.q
EOF

rm -f /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/jet_pp276HITrack_mergedpthatbins_2013MC.root
hadd /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/jet_pp276HITrack_mergedpthatbins_2013MC.root /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet*.root


