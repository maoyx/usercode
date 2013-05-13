//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jul  8 19:18:37 2012 by ROOT version 5.27/06b
// from TTree HltTree/
// found on file: /d102/yjlee/hiForest2MC/Pythia80_HydjetDrum_mix01_HiForest2_v20.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

class Skims {
public :
   Skims(){};
   ~Skims(){};

   // Declaration of leaf types
   Int_t           L1simulation_step;
   Int_t           reco_extra;
   Int_t           reco_extra_jet;
   Int_t           gen_step;
   Int_t           pat_step;
   Int_t           extrapatstep;
   Int_t           ana_step;
   Int_t           phltJetHI;
   Int_t           pcollisionEventSelection;
   Int_t           pPAcollisionEventSelectionPA;
   Int_t           pHBHENoiseFilter;
   Int_t           phiEcalRecHitSpikeFilter;
   Int_t           phfCoincFilter;
   Int_t           ppurityFractionFilter;

   // List of branches
   TBranch        *b_L1simulation_step;   //!
   TBranch        *b_reco_extra;   //!
   TBranch        *b_reco_extra_jet;   //!
   TBranch        *b_gen_step;   //!
   TBranch        *b_pat_step;   //!
   TBranch        *b_extrapatstep;   //!
   TBranch        *b_ana_step;   //!
   TBranch        *b_phltJetHI;   //!
   TBranch        *b_pcollisionEventSelection;   //!
   TBranch        *b_pPAcollisionEventSelectionPA;   //!
   TBranch        *b_pHBHENoiseFilter;   //!
   TBranch        *b_phiEcalRecHitSpikeFilter;   //!
   TBranch        *b_phfCoincFilter;   //!
   TBranch        *b_ppurityFractionFilter;   //!

};


void setupSkimTree(TTree *t,Skims &tSkims,bool doCheck = 1)
{
   // Set branch addresses and branch pointers
   if (t->GetBranch("L1simulation_step")) t->SetBranchAddress("L1simulation_step", &tSkims.L1simulation_step, &tSkims.b_L1simulation_step);
   if (t->GetBranch("reco_extra")) t->SetBranchAddress("reco_extra", &tSkims.reco_extra, &tSkims.b_reco_extra);
   if (t->GetBranch("reco_extra_jet")) t->SetBranchAddress("reco_extra_jet", &tSkims.reco_extra_jet, &tSkims.b_reco_extra_jet);
   if (t->GetBranch("gen_step")) t->SetBranchAddress("gen_step", &tSkims.gen_step, &tSkims.b_gen_step);
   if (t->GetBranch("pat_step")) t->SetBranchAddress("pat_step", &tSkims.pat_step, &tSkims.b_pat_step);
   if (t->GetBranch("extrapatstep")) t->SetBranchAddress("extrapatstep", &tSkims.extrapatstep, &tSkims.b_extrapatstep);
   if (t->GetBranch("ana_step")) t->SetBranchAddress("ana_step", &tSkims.ana_step, &tSkims.b_ana_step);
   if (t->GetBranch("phltJetHI")) t->SetBranchAddress("phltJetHI", &tSkims.phltJetHI, &tSkims.b_phltJetHI);
   if (t->GetBranch("pcollisionEventSelection")) t->SetBranchAddress("pcollisionEventSelection", &tSkims.pcollisionEventSelection, &tSkims.b_pcollisionEventSelection);
   if (t->GetBranch("pPAcollisionEventSelectionPA"))  t->SetBranchAddress("pPAcollisionEventSelectionPA", &tSkims.pPAcollisionEventSelectionPA, &tSkims.b_pPAcollisionEventSelectionPA);
   if (t->GetBranch("pHBHENoiseFilter")) t->SetBranchAddress("pHBHENoiseFilter", &tSkims.pHBHENoiseFilter, &tSkims.b_pHBHENoiseFilter);
   if (t->GetBranch("phiEcalRecHitSpikeFilter")) t->SetBranchAddress("phiEcalRecHitSpikeFilter", &tSkims.phiEcalRecHitSpikeFilter, &tSkims.b_phiEcalRecHitSpikeFilter);
   if (t->GetBranch("phfCoincFilter")) t->SetBranchAddress("phfCoincFilter", &tSkims.phfCoincFilter, &tSkims.b_phfCoincFilter);
   if (t->GetBranch("ppurityFractionFilter")) t->SetBranchAddress("ppurityFractionFilter", &tSkims.ppurityFractionFilter, &tSkims.b_ppurityFractionFilter);
   if (doCheck) {
   }
}

