root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt15/HiForest_v77_merged01/pt15_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppJet15.root",30,2.034E-01,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt30/HiForest_v77_merged01/pt30_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppJet30.root",50,1.075E-02,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt50/HiForest_v77_merged01/pt50_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppJet50.root",80,1.025E-03,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt80/HiForest_v77_merged01/pt80_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppJet80.root",120, 9.865E-05,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt120/HiForest_v77_merged01/pt120_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppJet120.root",170,1.129E-05,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt170/HiForest_v77_merged01/pt170_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppJet170.root",220,1.465E-06,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt220/HiForest_v77_merged01/pt220_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppJet220.root",280,2.837E-07,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt280/HiForest_v77_merged01/pt280_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppJet280.root",9999,5.323E-08,true)
.q
EOF

rm -f /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/jet_Corrpp276_mergedpthatbins_2013MC.root
hadd /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/jet_Corrpp276_mergedpthatbins_2013MC.root /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppJet*.root


root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt15/HiForest_v77_merged01/pt15_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppPuJet15.root",30,2.034E-01,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt30/HiForest_v77_merged01/pt30_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppPuJet30.root",50,1.075E-02,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt50/HiForest_v77_merged01/pt50_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppPuJet50.root",80,1.025E-03,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt80/HiForest_v77_merged01/pt80_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppPuJet80.root",120, 9.865E-05,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt120/HiForest_v77_merged01/pt120_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppPuJet120.root",170,1.129E-05,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt170/HiForest_v77_merged01/pt170_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppPuJet170.root",220,1.465E-06,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt220/HiForest_v77_merged01/pt220_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppPuJet220.root",280,2.837E-07,true)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt280/HiForest_v77_merged01/pt280_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppPuJet280.root",9999,5.323E-08,true)
.q
EOF

rm -f /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/PUJet_Corrpp276_mergedpthatbins_2013MC.root
hadd /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/PUJet_Corrpp276_mergedpthatbins_2013MC.root /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276CorrppPuJet*.root

