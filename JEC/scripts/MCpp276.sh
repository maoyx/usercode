root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt15/HiForest_v77_merged01/pt15_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet15.root",30,2.034E-01,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt30/HiForest_v77_merged01/pt30_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet30.root",50,1.075E-02,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt50/HiForest_v77_merged01/pt50_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet50.root",80,1.025E-03,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt80/HiForest_v77_merged01/pt80_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet80.root",120, 9.865E-05,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt120/HiForest_v77_merged01/pt120_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet120.root",170,1.129E-05,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt170/HiForest_v77_merged01/pt170_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet170.root",220,1.465E-06,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt220/HiForest_v77_merged01/pt220_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet220.root",280,2.837E-07,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt280/HiForest_v77_merged01/pt280_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet280.root",9999,5.323E-08,false)
.q
EOF

rm -f /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/jet_pp276_mergedpthatbins_2013MC.root
hadd /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/jet_pp276_mergedpthatbins_2013MC.root /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppJet*.root


root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt15/HiForest_v77_merged01/pt15_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppPuJet15.root",30,2.034E-01,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt30/HiForest_v77_merged01/pt30_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppPuJet30.root",50,1.075E-02,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt50/HiForest_v77_merged01/pt50_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppPuJet50.root",80,1.025E-03,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt80/HiForest_v77_merged01/pt80_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppPuJet80.root",120, 9.865E-05,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt120/HiForest_v77_merged01/pt120_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppPuJet120.root",170,1.129E-05,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt170/HiForest_v77_merged01/pt170_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppPuJet170.root",220,1.465E-06,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt220/HiForest_v77_merged01/pt220_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppPuJet220.root",280,2.837E-07,false)
.q
EOF
root -l -b <<EOF
.L loadfwlite.C
.x load_JEC.C 
.x weight_PtHatBins_IncPUJet.C+("/mnt/hadoop/cms/store/user/dgulhan/pp2013/P01/prod16/Signal_Pythia_pt280/HiForest_v77_merged01/pt280_P01_prod16_v77_merged_forest_0.root","/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppPuJet280.root",9999,5.323E-08,false)
.q
EOF

rm -f /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/PUJet_pp276_mergedpthatbins_2013MC.root
hadd /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/PUJet_pp276_mergedpthatbins_2013MC.root /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/combinePtHatBins/276ppPuJet*.root

