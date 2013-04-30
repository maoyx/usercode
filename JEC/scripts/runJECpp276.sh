#! /bin/sh -f
echo "Running for non-PU jet........... "
# Generate the response
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input jet_pp276_mergedpthatbins_2013MC.root -useweight true -algs ak3PFJetAnalyzer -drmax 0.3 -output jra_hiF_ak3PF_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input jet_pp276_mergedpthatbins_2013MC.root -useweight true -algs ak4PFJetAnalyzer -drmax 0.4 -output jra_hiF_ak4PF_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input jet_pp276_mergedpthatbins_2013MC.root -useweight true -algs ak5PFJetAnalyzer -drmax 0.5 -output jra_hiF_ak5PF_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input jet_pp276_mergedpthatbins_2013MC.root -useweight true -algs ak6PFJetAnalyzer -drmax 0.6 -output jra_hiF_ak6PF_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input jet_pp276_mergedpthatbins_2013MC.root -useweight true -algs ak2CaloJetAnalyzer -drmax 0.2 -output jra_hiF_ak2Calo_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input jet_pp276_mergedpthatbins_2013MC.root -useweight true -algs ak3CaloJetAnalyzer -drmax 0.3 -output jra_hiF_ak3Calo_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input jet_pp276_mergedpthatbins_2013MC.root -useweight true -algs ak4CaloJetAnalyzer -drmax 0.4 -output jra_hiF_ak4Calo_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input jet_pp276_mergedpthatbins_2013MC.root -useweight true -algs ak5CaloJetAnalyzer -drmax 0.5 -output jra_hiF_ak5Calo_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input jet_pp276_mergedpthatbins_2013MC.root -useweight true -algs ak6CaloJetAnalyzer -drmax 0.6 -output jra_hiF_ak6Calo_dijet.root
# Get the l3 corrections
jet_l3_correction_x -input jra_hiF_ak3PF_dijet.root -era JEC_dijet  -batch true -formats pdf -algs ak3PFJetAnalyzer -output jra_hiF_ak3PF_l3_dijet.root
jet_l3_correction_x -input jra_hiF_ak4PF_dijet.root -era JEC_dijet  -batch true -formats pdf -algs ak4PFJetAnalyzer -output jra_hiF_ak4PF_l3_dijet.root
jet_l3_correction_x -input jra_hiF_ak5PF_dijet.root -era JEC_dijet  -batch true -formats pdf -algs ak5PFJetAnalyzer -output jra_hiF_ak5PF_l3_dijet.root
jet_l3_correction_x -input jra_hiF_ak6PF_dijet.root -era JEC_dijet  -batch true -formats pdf -algs ak6PFJetAnalyzer -output jra_hiF_ak6PF_l3_dijet.root
jet_l3_correction_x -input jra_hiF_ak2Calo_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak2CaloJetAnalyzer -output jra_hiF_ak2Calo_l3_dijet.root
jet_l3_correction_x -input jra_hiF_ak3Calo_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak3CaloJetAnalyzer -output jra_hiF_ak3Calo_l3_dijet.root
jet_l3_correction_x -input jra_hiF_ak4Calo_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak4CaloJetAnalyzer -output jra_hiF_ak4Calo_l3_dijet.root
jet_l3_correction_x -input jra_hiF_ak5Calo_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak5CaloJetAnalyzer -output jra_hiF_ak5Calo_l3_dijet.root
jet_l3_correction_x -input jra_hiF_ak6Calo_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak6CaloJetAnalyzer -output jra_hiF_ak6Calo_l3_dijet.root

# use the l3 output and response output to get l2 corrections
jet_l2_correction_x -input jra_hiF_ak3PF_dijet.root -l3input jra_hiF_ak3PF_l3_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak3PFJetAnalyzer -output jra_hiF_ak3PF_l2_dijet.root
jet_l2_correction_x -input jra_hiF_ak4PF_dijet.root -l3input jra_hiF_ak4PF_l3_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak4PFJetAnalyzer -output jra_hiF_ak4PF_l2_dijet.root
jet_l2_correction_x -input jra_hiF_ak5PF_dijet.root -l3input jra_hiF_ak5PF_l3_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak5PFJetAnalyzer -output jra_hiF_ak5PF_l2_dijet.root
jet_l2_correction_x -input jra_hiF_ak6PF_dijet.root -l3input jra_hiF_ak6PF_l3_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak6PFJetAnalyzer -output jra_hiF_ak6PF_l2_dijet.root
jet_l2_correction_x -input jra_hiF_ak2Calo_dijet.root -l3input jra_hiF_ak2Calo_l3_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak2CaloJetAnalyzer -output jra_hiF_ak2Calo_l2_dijet.root
jet_l2_correction_x -input jra_hiF_ak3Calo_dijet.root -l3input jra_hiF_ak3Calo_l3_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak3CaloJetAnalyzer -output jra_hiF_ak3Calo_l2_dijet.root
jet_l2_correction_x -input jra_hiF_ak4Calo_dijet.root -l3input jra_hiF_ak4Calo_l3_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak4CaloJetAnalyzer -output jra_hiF_ak4Calo_l2_dijet.root
jet_l2_correction_x -input jra_hiF_ak5Calo_dijet.root -l3input jra_hiF_ak5Calo_l3_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak5CaloJetAnalyzer -output jra_hiF_ak5Calo_l2_dijet.root
jet_l2_correction_x -input jra_hiF_ak6Calo_dijet.root -l3input jra_hiF_ak6Calo_l3_dijet.root -era JEC_dijet -batch true -formats pdf -algs ak6CaloJetAnalyzer -output jra_hiF_ak6Calo_l2_dijet.root

#draw l2 corrections
jet_draw_l2_correction_x -alg ak3PFJetAnalyzer -filename jra_hiF_ak3PF_l2_dijet.root -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -outputFormat .pdf
jet_draw_l2_correction_x -alg ak4PFJetAnalyzer -filename jra_hiF_ak4PF_l2_dijet.root -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -outputFormat .pdf
jet_draw_l2_correction_x -alg ak5PFJetAnalyzer -filename jra_hiF_ak5PF_l2_dijet.root -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -outputFormat .pdf
jet_draw_l2_correction_x -alg ak3CaloJetAnalyzer -filename jra_hiF_ak3Calo_l2_dijet.root -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -outputFormat .pdf
jet_draw_l2_correction_x -alg ak4CaloJetAnalyzer -filename jra_hiF_ak4Calo_l2_dijet.root -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -outputFormat .pdf
jet_draw_l2_correction_x -alg ak5CaloJetAnalyzer -filename jra_hiF_ak5Calo_l2_dijet.root -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -outputFormat .pdf

# draw l2l3 corrections
jet_draw_corrections_x -algs ak3pf -path "" -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -era JEC_dijet -useL2Cor true -useL3Cor -true -outputFormat .pdf
jet_draw_corrections_x -algs ak4pf -path "" -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -era JEC_dijet -useL2Cor true -useL3Cor -true -outputFormat .pdf
jet_draw_corrections_x -algs ak5pf -path "" -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -era JEC_dijet -useL2Cor true -useL3Cor -true -outputFormat .pdf
jet_draw_corrections_x -algs ak3calo -path "" -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -era JEC_dijet -useL2Cor true -useL3Cor -true -outputFormat .pdf
jet_draw_corrections_x -algs ak4calo -path "" -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -era JEC_dijet -useL2Cor true -useL3Cor -true -outputFormat .pdf
jet_draw_corrections_x -algs ak5calo -path "" -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -era JEC_dijet -useL2Cor true -useL3Cor -true -outputFormat .pdf

# draw pt closure
jet_draw_closure_pt_x -algs ak3pf -path "" -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -era JEC_dijet -outputFormat .pdf


echo "Running for PuJet........... "
# Generate the response
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input PUJet_pp276_mergedpthatbins_2013MC.root -useweight true -algs akPu3PFJetAnalyzer -drmax 0.3 -output jra_hiF_akPu3PF_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input PUJet_pp276_mergedpthatbins_2013MC.root -useweight true -algs akPu4PFJetAnalyzer -drmax 0.4 -output jra_hiF_akPu4PF_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input PUJet_pp276_mergedpthatbins_2013MC.root -useweight true -algs akPu5PFJetAnalyzer -drmax 0.5 -output jra_hiF_akPu5PF_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input PUJet_pp276_mergedpthatbins_2013MC.root -useweight true -algs akPu6PFJetAnalyzer -drmax 0.6 -output jra_hiF_akPu6PF_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input PUJet_pp276_mergedpthatbins_2013MC.root -useweight true -algs akPu2CaloJetAnalyzer -drmax 0.2 -output jra_hiF_akPu2Calo_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input PUJet_pp276_mergedpthatbins_2013MC.root -useweight true -algs akPu3CaloJetAnalyzer -drmax 0.3 -output jra_hiF_akPu3Calo_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input PUJet_pp276_mergedpthatbins_2013MC.root -useweight true -algs akPu4CaloJetAnalyzer -drmax 0.4 -output jra_hiF_akPu4Calo_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input PUJet_pp276_mergedpthatbins_2013MC.root -useweight true -algs akPu5CaloJetAnalyzer -drmax 0.5 -output jra_hiF_akPu5Calo_dijet.root
jet_response_analyzer_x  largeEta_binning_ptRebin.config -input PUJet_pp276_mergedpthatbins_2013MC.root -useweight true -algs akPu6CaloJetAnalyzer -drmax 0.6 -output jra_hiF_akPu6Calo_dijet.root
# Get the l3 corrections
jet_l3_correction_x -input jra_hiF_akPu3PF_dijet.root -era JEC_akPu3PF  -batch true -formats pdf -algs akPu3PFJetAnalyzer -output jra_hiF_akPu3PF_l3_dijet.root
jet_l3_correction_x -input jra_hiF_akPu4PF_dijet.root -era JEC_akPu4PF  -batch true -formats pdf -algs akPu4PFJetAnalyzer -output jra_hiF_akPu4PF_l3_dijet.root
jet_l3_correction_x -input jra_hiF_akPu5PF_dijet.root -era JEC_akPu5PF  -batch true -formats pdf -algs akPu5PFJetAnalyzer -output jra_hiF_akPu5PF_l3_dijet.root
jet_l3_correction_x -input jra_hiF_akPu6PF_dijet.root -era JEC_akPu6PF  -batch true -formats pdf -algs akPu6PFJetAnalyzer -output jra_hiF_akPu6PF_l3_dijet.root
jet_l3_correction_x -input jra_hiF_akPu2Calo_dijet.root -era JEC_akPu2Calo -batch true -formats pdf -algs akPu2CaloJetAnalyzer -output jra_hiF_akPu2Calo_l3_dijet.root
jet_l3_correction_x -input jra_hiF_akPu3Calo_dijet.root -era JEC_akPu3Calo -batch true -formats pdf -algs akPu3CaloJetAnalyzer -output jra_hiF_akPu3Calo_l3_dijet.root
jet_l3_correction_x -input jra_hiF_akPu4Calo_dijet.root -era JEC_akPu4Calo -batch true -formats pdf -algs akPu4CaloJetAnalyzer -output jra_hiF_akPu4Calo_l3_dijet.root
jet_l3_correction_x -input jra_hiF_akPu5Calo_dijet.root -era JEC_akPu5Calo -batch true -formats pdf -algs akPu5CaloJetAnalyzer -output jra_hiF_akPu5Calo_l3_dijet.root
jet_l3_correction_x -input jra_hiF_akPu6Calo_dijet.root -era JEC_akPu6Calo -batch true -formats pdf -algs akPu6CaloJetAnalyzer -output jra_hiF_akPu6Calo_l3_dijet.root

# use the l3 output and response output to get l2 corrections
jet_l2_correction_x -input jra_hiF_akPu3PF_dijet.root -l3input jra_hiF_akPu3PF_l3_dijet.root -era JEC_akPu3PF -batch true -formats pdf -algs akPu3PFJetAnalyzer -output jra_hiF_akPu3PF_l2_dijet.root
jet_l2_correction_x -input jra_hiF_akPu4PF_dijet.root -l3input jra_hiF_akPu4PF_l3_dijet.root -era JEC_akPu4PF -batch true -formats pdf -algs akPu4PFJetAnalyzer -output jra_hiF_akPu4PF_l2_dijet.root
jet_l2_correction_x -input jra_hiF_akPu5PF_dijet.root -l3input jra_hiF_akPu5PF_l3_dijet.root -era JEC_akPu5PF -batch true -formats pdf -algs akPu5PFJetAnalyzer -output jra_hiF_akPu5PF_l2_dijet.root
jet_l2_correction_x -input jra_hiF_akPu6PF_dijet.root -l3input jra_hiF_akPu6PF_l3_dijet.root -era JEC_akPu6PF -batch true -formats pdf -algs akPu6PFJetAnalyzer -output jra_hiF_akPu6PF_l2_dijet.root
jet_l2_correction_x -input jra_hiF_akPu2Calo_dijet.root -l3input jra_hiF_akPu2Calo_l3_dijet.root -era JEC_akPu2Calo -batch true -formats pdf -algs akPu2CaloJetAnalyzer -output jra_hiF_akPu2Calo_l2_dijet.root
jet_l2_correction_x -input jra_hiF_akPu3Calo_dijet.root -l3input jra_hiF_akPu3Calo_l3_dijet.root -era JEC_akPu3Calo -batch true -formats pdf -algs akPu3CaloJetAnalyzer -output jra_hiF_akPu3Calo_l2_dijet.root
jet_l2_correction_x -input jra_hiF_akPu4Calo_dijet.root -l3input jra_hiF_akPu4Calo_l3_dijet.root -era JEC_akPu4Calo -batch true -formats pdf -algs akPu4CaloJetAnalyzer -output jra_hiF_akPu4Calo_l2_dijet.root
jet_l2_correction_x -input jra_hiF_akPu5Calo_dijet.root -l3input jra_hiF_akPu5Calo_l3_dijet.root -era JEC_akPu5Calo -batch true -formats pdf -algs akPu5CaloJetAnalyzer -output jra_hiF_akPu5Calo_l2_dijet.root
jet_l2_correction_x -input jra_hiF_akPu6Calo_dijet.root -l3input jra_hiF_akPu6Calo_l3_dijet.root -era JEC_akPu6Calo -batch true -formats pdf -algs akPu6CaloJetAnalyzer -output jra_hiF_akPu6Calo_l2_dijet.root

#draw l2 corrections
jet_draw_l2_correction_x -alg akPu3PFJetAnalyzer -filename jra_hiF_akPu3PF_l2_dijet.root -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -outputFormat .pdf
jet_draw_l2_correction_x -alg akPu4PFJetAnalyzer -filename jra_hiF_akPu4PF_l2_dijet.root -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -outputFormat .pdf
jet_draw_l2_correction_x -alg akPu5PFJetAnalyzer -filename jra_hiF_akPu5PF_l2_dijet.root -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -outputFormat .pdf
jet_draw_l2_correction_x -alg akPu3CaloJetAnalyzer -filename jra_hiF_akPu3Calo_l2_dijet.root -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -outputFormat .pdf
jet_draw_l2_correction_x -alg akPu4CaloJetAnalyzer -filename jra_hiF_akPu4Calo_l2_dijet.root -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -outputFormat .pdf
jet_draw_l2_correction_x -alg akPu5CaloJetAnalyzer -filename jra_hiF_akPu5Calo_l2_dijet.root -outputDir /net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/testCorrections/ -outputFormat .pdf



