#ifndef hiForest_h
#define hiForest_h

#include <iostream>
#include <vector>
#include <algorithm>

#include "commonTool.h"
#include "SetupPhotonTree.h"
#include "SetupJetTree.h"
#include "SetupHltTree.h"
#include "SetupSkimTree.h"
#include "SetupTrackTree.h"
#include "SetupHitTree.h"
#include "SetupEvtTree.h"
#include "SetupMetTree.h"
#include "SetupGenpTree.h"
#include "SetupPFTree.h"
#include "SetupGenParticleTree.h"
#include "TrackingParam.h"
//#include "TrackingCorrectionsv6.h"
#include "TrackingCorrections4Cent.h"
#include <TTree.h>
#include <TFile.h>
#include <TString.h>
#include <TF1.h>
#include <TCut.h>

#include "DummyJetCorrector.h"

// ==========================================================
// Main class which can be used to read the hiForest trees
//
// Author: Yen-Jie Lee
//
// ==========================================================

namespace names{
  enum Algo{
  };
  string AlgoRename[100] = {
  };
  string AlgoAnalyzer[100] = {
  };
}

class HiForest : public TNamed
{

  public: 
   HiForest(const char *file, const char *name="forest", bool ispp = 0, bool ismc = 0, bool isrecorrected = 0);
  ~HiForest();

  //==================================================================================================================================
  // Utility functions
  //==================================================================================================================================
  void CheckTree(TTree *t,const char *title);			// Check the status of a tree
  void CheckArraySizes();					// Check if array size is large enough
  void GetEntry(int i);
  int  GetEntries();  						// Get the number of entries 
  void GetEnergyScaleTable(char *fileNameTable);                // Get photon energy scale table
  void InitTree();						// Initialize track correction
  void PrintStatus();						// Print the status of the hiForest
  void SetOutputFile(const char *name);               		// Set output file name for skim
  void AddCloneTree(TTree* t, const char *dirName, const char *treeName);   // Add a clone tree to the clone forest
  void FillOutput();						// Fill output forest  
  
  Long64_t Draw(const char* varexp, const char* selection, Option_t* option = "", Long64_t nentries = 1000000000, Long64_t firstentry = 0){
     return tree->Draw(varexp,selection,option,nentries,firstentry);
  }

  //==================================================================================================================================
  // Event filtering utility functions
  //==================================================================================================================================
  bool selectEvent();
  TCut eventSelection();

  //==================================================================================================================================
  // Photon utility functions
  //==================================================================================================================================
  bool isSpike(int i);                          		// return true if it is considered as a spike candidate
  bool isLooseEGamma(int i);
  bool isLoosePhoton(int i);
  bool isGoodPhoton(int i);                   			// return true if it is considered as a hiGoodPhoton candidate
  bool isIsolatedPhoton(int i);                     		// return true if it is considered as a hiGoodPhoton candidate
  bool isMCSignal(int i);
  bool isDirectPhoton(int i);
  bool isFragPhoton(int i);
  float getCorrEt(int j);                                       // Get Corrected photon Et

  //==================================================================================================================================
  // Jet utility functions
  //==================================================================================================================================
  void sortJets(TTree* jetTree, Jets& jets, double etaMax = 2, double ptMin = 40, bool allEvents = 1, int smearType = -1);
  int leadingJet();
  int subleadingJet();
  int thirdJet();
  double deltaPhiDijet(Jets& jets);
  bool hasDiJet(Jets& jets, double pt1 = 100, double pt2 = 40, double dphiMin = 2.*3.1415926/3.);
  void fakeRejection(TTree *jetTree, Jets &jets, bool allEvents);

  //==================================================================================================================================
  // Track utility functions
  //==================================================================================================================================
  int getMatchedCaloTowerAllowReuse(int j);
  int getMatchedHBHEAllowReuse(int j);
  void matchTrackCalo(bool allEvents = 1);
  double getTrackCorrectionPara(int j);
  double getTrackCorrection(int j);
  
  //==================================================================================================================================
  // Get track-jet correlated variables. Not needed if correlatePF is run.
  //==================================================================================================================================
  void correlateTracks(TTree* jetTree, Jets& jets, bool allEvents = 1, bool smeared = 0);

  // TFile
  TFile *inf; 					// Input file 
  TFile *outf;                                  // Output file if we want to export the forest

  // Trees
  TTree *photonTree;				// Photon Tree, see branches in SetupPhotonTree.h
  TTree *icPu5jetTree;				// Jet Tree with icPu5 algorithm, see branches in SetupJetTree.h
  TTree *akPu2jetTree;				// Jet Tree with akPu2PF algorithm, see branches in SetupJetTree.h
  TTree *akPu3jetTree;				// Jet Tree with akPu3PF algorithm, see branches in SetupJetTree.h
  TTree *ak3jetTree;				// Jet Tree with ak3PF algorithm, see branches in SetupJetTree.h
  TTree *akPu4jetTree;				// Jet Tree with akPu4PF algorithm, see branches in SetupJetTree.h
  TTree *akPu2CaloJetTree;			// Jet Tree with akPu2Calo algorithm, see branches in SetupJetTree.h
  TTree *akPu3CaloJetTree;			// Jet Tree with akPu3Calo algorithm, see branches in SetupJetTree.h
  TTree *akPu4CaloJetTree;		        // Jet Tree with akPu4Calo algorithm, see branches in SetupJetTree.h
  TTree *hltTree;				// OpenHLT Tree, see branches in SetupHltTree.h
  TTree *trackTree;				// Track Tree, see branches in SetupTrackTree.h
  TTree *pixtrackTree;				// Track Tree, see branches in SetupTrackTree.h
  TTree *skimTree;				// Skim Tree, contains event selection info, see branches in SetupSkimTree.h
  TTree *towerTree;                             // Tower Tree
  TTree *hbheTree;                              // HCAL HBHE Tree
  TTree *ebTree;                                // ECAL eb Tree
  TTree *evtTree;                               // Event Tree
  TTree *metTree;                               // MET Tree
  TTree *pfTree;                                // PF candidate Tree, see branches in SetupPFTree.h
  TTree *genpTree;                              // Gen photon of the signal event Tree
  TTree *genParticleTree;                       // Stable Gen particles
  
  TTree *tree;					// Pointer to the available tree, all trees in the forest are friended to each other

  vector<TTree*> jetTrees;
  vector<TTree*> cloneForest;                   // Vector of clones for skim

  TF1* fGauss;

  // Branches
  Hlts hlt;
  Skims skim;

  vector<Jets> alljets;

  Jets icPu5;
  Jets akPu2PF;
  Jets akPu3PF;
  Jets ak3PF;
  Jets akPu4PF;
  Jets akPu2Calo;
  Jets akPu3Calo;
  Jets akPu4Calo;

  Photons photon;
  Tracks track;
  Tracks pixtrack;
  Hits tower;
  Hits hbhe;
  Hits eb;
  Evts evt;
  Mets met;
  PFs pf;
  Genps genp;
  GenParticles genparticle;
    
  // Booleans
  bool hasPhotonTree;
  bool hasEvtTree;
  bool hasMetTree;
  bool hasPFTree;
  bool hasIcPu5JetTree;
  bool hasAkPu2JetTree;
  bool hasAkPu3JetTree;
  bool hasAk3JetTree;
  bool hasAkPu4JetTree;
  bool hasAkPu2CaloJetTree;
  bool hasAkPu3CaloJetTree;
  bool hasAkPu4CaloJetTree;
  bool hasHltTree;
  bool hasTrackTree;
  bool hasPixTrackTree;
  bool hasSkimTree;
  bool hasTowerTree;
  bool hasHbheTree;
  bool hasEbTree;
  bool hasGenpTree;
  bool hasGenParticleTree;

  bool setupOutput;
  bool verbose;
  bool pp;
  bool mc;
  bool doJetCorrection;
  bool doTrackCorrections;
  bool doTrackingSeparateLeadingSubleading;

  // Extra variables
  Float_t* towerEt;
  Float_t* towerdR;
  Float_t* hbheEt;
  Float_t* hbhedR;

  Float_t* jtChg;
  Float_t* jtNeut;
  Float_t* jtEM;

  Float_t* jtChgGen;
  Float_t* jtNeutGen;
  Float_t* jtEMGen;

  Float_t* jtPtMax;
  Float_t* jtPtMean;
  Float_t* jtPtMeanW;

  Int_t* jtLeadType;

  Int_t jtLead;
  Int_t jtSubLead;
  bool jtHasDijet;
  bool jtHasLeadingJet;

  Float_t* tjDeltaEtaLead;
  Float_t* tjDeltaPhiLead;
  Float_t* zLead;

  Float_t* tjDeltaEtaSubLead;
  Float_t* tjDeltaPhiSubLead;
  Float_t* zSubLead;

  Float_t* zOldLead;
  Float_t* zOldSubLead;

  Float_t* zSingleLead;
  Float_t* zLabLead;
  Float_t* zSingleSubLead;
  Float_t* zLabSubLead;
  Float_t* tjDeltaThetaLead;
  Float_t* tjDeltaThetaLabLead;
  Float_t* tjDeltaThetaSingleLead;

  Float_t* tjDeltaThetaSubLead;
  Float_t* tjDeltaThetaLabSubLead;
  Float_t* tjDeltaThetaSingleSubLead;

  Float_t* corrLead;
  Float_t* corrSubLead;

  Float_t minJetPtForTrkCor;
  Float_t leadingJetPtForTrkCor;
  Float_t subleadingJetPtForTrkCor;
  Float_t leadingJetEtaForTrkCor;
  Float_t subleadingJetEtaForTrkCor;
  Float_t leadingJetPhiForTrkCor;
  Float_t subleadingJetPhiForTrkCor;

  int nEntries;
  int currentEvent;
  double cone;
  bool initialized;

  vector<JetCorrectorParameters> vpar_HI310x;
  FactorizedJetCorrector *_JEC_HI310X;

//  TrackingParam trackCorrFromParam("TrkCorr_3DHistos_pThat80_Inclusive.root","CentralityWeights.root","PtResidualWeights.root");
  TrackingParam *trackCorrFromParam;

  vector<TrackingCorrections*> trackCorrections;
  TF1* fEnergyScale[2][10];  // [a][b],  a =0 for unconverted,  a=1 for converted.   b: 1,2,3 is centrality bin. b=0 is empty                                                                                       
  
 private:
  
  
  
};

HiForest::HiForest(const char *infName, const char* name, bool ispp, bool ismc, bool recjec):
   tree(0),
   fGauss(0),
   verbose(0),
   pp(ispp),
   mc(ismc),
   nEntries(0),
   currentEvent(0)
{
   tree = new TTree("tree","");
  SetName(name);
  // Input file
  inf = TFile::Open(infName);

  cone = 0.3;
  doTrackCorrections = 0;

  // Track correction initialized?
  initialized = 0;

  // Load trees. Hard coded for the moment
  hltTree          = (TTree*) inf->Get("hltanalysis/HltTree");
  skimTree         = (TTree*) inf->Get("skimanalysis/HltTree");
  photonTree       = (TTree*) inf->Get("multiPhotonAnalyzer/photon");
  trackTree        = (TTree*) inf->Get("anaTrack/trackTree");
//  trackTree        = (TTree*) inf->Get("ppTrack/trackTree");
  towerTree        = (TTree*) inf->Get("rechitanalyzer/tower");
  icPu5jetTree     = (TTree*) inf->Get("icPu5JetAnalyzer/t");
  akPu2jetTree     = (TTree*) inf->Get("akPu2PFJetAnalyzer/t");
  akPu3jetTree     = (TTree*) inf->Get("akPu3PFJetAnalyzer/t");
  ak3jetTree     = (TTree*) inf->Get("ak3PFJetAnalyzer/t");
  akPu4jetTree     = (TTree*) inf->Get("akPu4PFJetAnalyzer/t");
  akPu2CaloJetTree = (TTree*) inf->Get("akPu2CaloJetAnalyzer/t");
  akPu3CaloJetTree = (TTree*) inf->Get("akPu3CaloJetAnalyzer/t");
  akPu4CaloJetTree = (TTree*) inf->Get("akPu4CaloJetAnalyzer/t");
  hbheTree         = (TTree*) inf->Get("rechitanalyzer/hbhe");
  ebTree           = (TTree*) inf->Get("rechitanalyzer/eb");
  evtTree          = (TTree*) inf->Get("hiEvtAnalyzer/HiTree");
  metTree          = (TTree*) inf->Get("anaMET/metTree");
  pfTree           = (TTree*) inf->Get("pfcandAnalyzer/pfTree");
  genpTree         = (TTree*) inf->Get("genpana/photon");
  genParticleTree  = (TTree*) inf->Get("HiGenParticleAna/hi");
  
  // doesn't load genParticle by default
  //genParticleTree  = (TTree*) inf->Get("HiGenParticleAna/hi");

  // Check the validity of the trees.
  hasPhotonTree        = (photonTree       	!= 0);
  hasPFTree            = (pfTree     		!= 0);
  hasEvtTree           = (evtTree      		!= 0);
  hasMetTree           = (metTree      		!= 0);
  hasIcPu5JetTree      = (icPu5jetTree 		!= 0);
  hasAkPu2JetTree      = (akPu2jetTree 		!= 0);
  hasAkPu3JetTree      = (akPu3jetTree 		!= 0);
  hasAk3JetTree      = (ak3jetTree 		!= 0);
  hasAkPu4JetTree      = (akPu4jetTree 		!= 0);
  hasAkPu2CaloJetTree  = (akPu2CaloJetTree 	!= 0);
  hasAkPu3CaloJetTree  = (akPu3CaloJetTree 	!= 0);
  hasAkPu4CaloJetTree  = (akPu4CaloJetTree 	!= 0);
  hasTrackTree     = (trackTree    		!= 0);
  hasHltTree       = (hltTree    		!= 0);
  hasSkimTree      = (skimTree   		!= 0);
  hasTowerTree     = (towerTree    		!= 0);
  hasHbheTree      = (hbheTree     		!= 0);
  hasEbTree        = (ebTree       		!= 0);
  hasGenpTree	   = (genpTree     		!=0);
  hasGenParticleTree = (genParticleTree   	!=0);
  setupOutput = false;
  
  // Setup branches. See also Setup*.h
  if (hasPhotonTree) {
    photonTree->SetName("photon");
    photonTree->SetAlias("swiss","1-(eRight+eLeft+eTop+eBottom)/eMax");
    if (tree == 0) tree = photonTree; else tree->AddFriend(photonTree);
    setupPhotonTree(photonTree,photon);
  }

  if (hasPFTree) {
    pfTree->SetName("pf");
    if (tree == 0) tree = pfTree; else tree->AddFriend(pfTree);
    setupPFTree(pfTree,pf);
  }


  if (hasEvtTree) {
    evtTree->SetName("evtTree");
    if (tree == 0) tree = evtTree; else tree->AddFriend(evtTree);
    setupEvtTree(evtTree,evt);
  }

  if (hasMetTree) {
    metTree->SetName("met");
    if (tree == 0) tree = metTree; else tree->AddFriend(metTree);
    setupMetTree(metTree,met);
  }

  if (hasHltTree) {
    hltTree->SetName("hltTree");
    if (tree == 0) tree = hltTree; else tree->AddFriend(hltTree);
    setupHltTree(hltTree,hlt);
  }

  if (hasIcPu5JetTree) {
    icPu5jetTree->SetName("icPu5");
    if (tree == 0) tree = icPu5jetTree; else tree->AddFriend(icPu5jetTree);
    setupJetTree(icPu5jetTree,icPu5);
  }

  if (hasAkPu2JetTree) {
    akPu2jetTree->SetName("akPu2PF");
    if (tree == 0) tree = akPu2jetTree; else tree->AddFriend(akPu2jetTree);
    setupJetTree(akPu2jetTree,akPu2PF);
  }

  if (hasAkPu3JetTree) {
    akPu3jetTree->SetName("akPu3PF");
    if (tree == 0) tree = akPu3jetTree; else tree->AddFriend(akPu3jetTree);
    setupJetTree(akPu3jetTree,akPu3PF);
  }
  if (hasAk3JetTree) {
    ak3jetTree->SetName("ak3PF");
    if (tree == 0) tree = ak3jetTree; else tree->AddFriend(ak3jetTree);
    setupJetTree(ak3jetTree,ak3PF);
  }
  if (hasAkPu4JetTree) {
    akPu4jetTree->SetName("akPu4PF");
    if (tree == 0) tree = akPu4jetTree; else tree->AddFriend(akPu4jetTree);
    setupJetTree(akPu4jetTree,akPu4PF);
  }
  if (hasAkPu2CaloJetTree) {
    akPu2CaloJetTree->SetName("akPu2Calo");
    if (tree == 0) tree = akPu2CaloJetTree; else tree->AddFriend(akPu2CaloJetTree);
    setupJetTree(akPu2CaloJetTree,akPu2Calo);
  }

  if (hasAkPu3CaloJetTree) {
    akPu3CaloJetTree->SetName("akPu3Calo");
    if (tree == 0) tree = akPu3CaloJetTree; else tree->AddFriend(akPu3CaloJetTree);
    setupJetTree(akPu3CaloJetTree,akPu3Calo);
  }

  if (hasAkPu4CaloJetTree) {
    akPu4CaloJetTree->SetName("akPu4Calo");
    if (tree == 0) tree = akPu4CaloJetTree; else tree->AddFriend(akPu4CaloJetTree);
    setupJetTree(akPu4CaloJetTree,akPu4Calo);
  }


  if (hasTrackTree) {
    trackTree->SetName("track");
    if (tree == 0) tree = trackTree; else tree->AddFriend(trackTree);
    trackTree->SetAlias("mergedGeneral","(trkAlgo<4||(highPurity))");
    trackTree->SetAlias("mergedSelected","(trkAlgo<4||(highPurity&&trkAlgo==4)))");
    setupTrackTree(trackTree,track);
  }
   
  if (hasPixTrackTree) {
    pixtrackTree->SetName("pixtrack");
    if (tree == 0) tree = pixtrackTree; else tree->AddFriend(pixtrackTree);
    setupTrackTree(pixtrackTree,pixtrack);
  }
   
  if (hasSkimTree) {
    skimTree->SetName("skim");
    if (tree == 0) tree = skimTree; else tree->AddFriend(skimTree);
    setupSkimTree(skimTree,skim);
  }

  if (hasTowerTree) {
    towerTree->SetName("tower");
    if (tree == 0) tree = towerTree; else tree->AddFriend(towerTree);
    setupHitTree(towerTree,tower);
  }

  if (hasHbheTree) {
    hbheTree->SetName("hbhe");
    if (tree == 0) tree = hbheTree; else tree->AddFriend(hbheTree);
    setupHitTree(hbheTree,hbhe);
  }

  if (hasEbTree) {
    ebTree->SetName("eb");
    if (tree == 0) tree = ebTree; else tree->AddFriend(ebTree);
    setupHitTree(ebTree,eb);
  }
  
  if (hasGenpTree) {
    genpTree->SetName("genp");
    if (tree == 0) tree = genpTree; else tree->AddFriend(genpTree);
    setupGenpTree(genpTree,genp);
  }

  if (hasGenParticleTree) {
    genParticleTree->SetName("genParticle");
    if (tree == 0) tree = genParticleTree; else tree->AddFriend(genParticleTree);
    setupGenParticleTree(genParticleTree,genparticle);
  }

  tree->SetMarkerStyle(20);

  // Print the status of thre forest
  PrintStatus();

}

HiForest::~HiForest()
{
  if (setupOutput) {
    for (unsigned int i=0; i<cloneForest.size(); i++)
    {
      cloneForest[i]->AutoSave();
    }
    outf->Close();
  }
}

void HiForest::GetEntry(int i)
{

   currentEvent = i;
  // get the entry of the available trees
  if (hasPhotonTree)   photonTree   ->GetEntry(i);
  if (hasPFTree)       pfTree   ->GetEntry(i);
  if (hasEvtTree)      evtTree   ->GetEntry(i);
  if (hasMetTree)      metTree   ->GetEntry(i);
  if (hasHltTree)      hltTree      ->GetEntry(i);
  if (hasSkimTree)     skimTree     ->GetEntry(i);
  if (hasIcPu5JetTree) icPu5jetTree ->GetEntry(i);
  if (hasAkPu2JetTree) akPu2jetTree ->GetEntry(i);
  if (hasAkPu3JetTree) akPu3jetTree ->GetEntry(i);
  if (hasAk3JetTree) ak3jetTree ->GetEntry(i);
  if (hasAkPu4JetTree) akPu4jetTree ->GetEntry(i);
  if (hasAkPu2CaloJetTree) akPu2CaloJetTree ->GetEntry(i);
  if (hasAkPu3CaloJetTree) akPu3CaloJetTree ->GetEntry(i);
  if (hasAkPu4CaloJetTree) akPu4CaloJetTree ->GetEntry(i);
  if (hasTrackTree)    trackTree    ->GetEntry(i);
  if (hasPixTrackTree) pixtrackTree ->GetEntry(i);
  if (hasTowerTree)    towerTree    ->GetEntry(i);
  if (hasHbheTree)     hbheTree     ->GetEntry(i);
  if (hasEbTree)       ebTree     ->GetEntry(i);
  if (hasGenpTree)     genpTree   ->GetEntry(i);
  if (hasGenParticleTree) genParticleTree   ->GetEntry(i);

  minJetPtForTrkCor = 40;
  leadingJetPtForTrkCor = -100;
  subleadingJetPtForTrkCor = -100;
  leadingJetEtaForTrkCor = -100;
  subleadingJetEtaForTrkCor = -100;
  leadingJetPhiForTrkCor = -100;
  subleadingJetPhiForTrkCor = -100;
}

int HiForest::GetEntries()
{
  // get the entries of the available trees
  return nEntries;
}


void HiForest::InitTree()
{
   // Setup Track Corrections 	 
   if(doTrackCorrections){

      trackCorrFromParam = new TrackingParam();

      if (!pp) {
         trackCorrections.push_back(new TrackingCorrections("Forest2STAv14","Forest2_MergedGeneral_jetfine"));
      } else {
         trackCorrections.push_back(new TrackingCorrections("Forest2STApp","Forest2_MergedGeneral_jetfine"));
      }
//       trackCorrections.push_back(new TrackingCorrections("Forest2STAv12","Forest2_MergedGeneral_j1"));
//       trackCorrections.push_back(new TrackingCorrections("Forest2STAv12","Forest2_MergedGeneral_j2"));

      for(int i = 0; i < trackCorrections.size(); ++i){
          if (!pp) {
              trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/TrkCorrtest_hy18dj50_akPu3PF_100_-1_-1000_genJetMode0.root",50);
              trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/TrkCorrtest_hy18dj80_akPu3PF_100_-1_-1000_genJetMode0.root",80);
              trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/TrkCorrtest_hy18dj100_akPu3PF_100_-1_-1000_genJetMode0.root",100);
              trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/TrkCorrtest_hy18dj120_akPu3PF_100_-1_-1000_genJetMode0.root",120);
              trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/TrkCorrtest_hy18dj170_akPu3PF_100_-1_-1000_genJetMode0.root",170);
              trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/TrkCorrtest_hy18dj200_akPu3PF_100_-1_-1000_genJetMode0.root",200);
              trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/TrkCorrtest_hy18dj250_akPu3PF_100_-1_-1000_genJetMode0.root",250);
             trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec/TrkCorrtest_hy18dj300_akPu3PF_100_-1_-1000_genJetMode0.root",300);
              // v12: frozen for preapproval
              //            trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv12XSec/IterTrkCorrv12XSec_hy18dj100_akPu3PF_100_40_2749_genJetMode0.root",100);
              //            trackCorrections[i]->AddNormFile("trkcorr/IterTrkCorrv12XSec/IterTrkCorrv12XSec_hy18dj100_akPu3PF_-1_-1_-1000_genJetMode0.root");
              //            trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv12XSec/IterTrkCorrv12XSec_hy18dj200_akPu3PF_100_40_2749_genJetMode0.root",200);
              //            trackCorrections[i]->AddNormFile("trkcorr/IterTrkCorrv12XSec/IterTrkCorrv12XSec_hy18dj200_akPu3PF_-1_-1_-1000_genJetMode0.root");
          } else {
     /*      trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_ppSigdj50_ak3PF_100_-1_-1000_genJetMode0.root",50);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_ppSigdj80_ak3PF_100_-1_-1000_genJetMode0.root",80);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_ppSigdj120_ak3PF_100_-1_-1000_genJetMode0.root",120);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_ppSigdj170_ak3PF_100_-1_-1000_genJetMode0.root",170);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_ppSigdj200_ak3PF_100_-1_-1000_genJetMode0.root",200);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_ppSigdj250_ak3PF_100_-1_-1000_genJetMode0.root",250);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_ppSigdj300_ak3PF_100_-1_-1000_genJetMode0.root",300);
         
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/IterTrkCorrv14XSec_pp_sigdj80to120_akPu3PF_100_-1_-1000_genJetMode0.root",80);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/IterTrkCorrv14XSec_pp_sigdj120to170_akPu3PF_100_-1_-1000_genJetMode0.root",120);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/IterTrkCorrv14XSec_pp_sigdj170to200_akPu3PF_100_-1_-1000_genJetMode0.root",170);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/IterTrkCorrv14XSec_pp_sigdj200to250_akPu3PF_100_-1_-1000_genJetMode0.root",200);
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/IterTrkCorrv14XSec_pp_sigdj250to9999_akPu3PF_100_-1_-1000_genJetMode0.root",250); 
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp_Sigdj30_ak3PF_100_-1_-1000_genJetMode0.root", 30) ;
       
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp_Sigdj50_ak3PF_100_-1_-1000_genJetMode0.root", 50) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp_Sigdj80_ak3PF_100_-1_-1000_genJetMode0.root", 80) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp_Sigdj120_ak3PF_100_-1_-1000_genJetMode0.root", 120) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp_Sigdj170_ak3PF_100_-1_-1000_genJetMode0.root", 170) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp_Sigdj220_ak3PF_100_-1_-1000_genJetMode0.root", 220) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp_Sigdj280_ak3PF_100_-1_-1000_genJetMode0.root", 280) ;
      
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp2012HITracking_Sigdj50_ak3PF_100_-1_-1000_genJetMode0.root", 50) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp2012HITracking_Sigdj80_ak3PF_100_-1_-1000_genJetMode0.root", 80) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp2012HITracking_Sigdj120_ak3PF_100_-1_-1000_genJetMode0.root", 120) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp2012HITracking_Sigdj170_ak3PF_100_-1_-1000_genJetMode0.root", 170) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp2012HITracking_Sigdj200_ak3PF_100_-1_-1000_genJetMode0.root", 200) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp2012HITracking_Sigdj250_ak3PF_100_-1_-1000_genJetMode0.root", 250) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp2012HITracking_Sigdj300_ak3PF_100_-1_-1000_genJetMode0.root", 300) ;
      */ 
          trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp2013HITracking_Sigdj50_ak3PF_100_-1_-1000_genJetMode0.root", 50) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp2013HITracking_Sigdj80_ak3PF_100_-1_-1000_genJetMode0.root", 80) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp2013HITracking_Sigdj120_ak3PF_100_-1_-1000_genJetMode0.root", 120) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp2013HITracking_Sigdj170_ak3PF_100_-1_-1000_genJetMode0.root", 170) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp2013HITracking_Sigdj220_ak3PF_100_-1_-1000_genJetMode0.root", 220) ;
           trackCorrections[i]->AddSample("trkcorr/IterTrkCorrv14XSec_pp/TrkCorrtest_pp2013HITracking_Sigdj280_ak3PF_100_-1_-1000_genJetMode0.root", 280) ;
        }
         trackCorrections[i]->weightSamples_ = true;
         trackCorrections[i]->smoothLevel_ = 0;
         trackCorrections[i]->trkPhiMode_ = false;
         trackCorrections[i]->ppMode_ = pp;
         trackCorrections[i]->Init();
      }
      minJetPtForTrkCor = 40;
      initialized = 1;
      if (doTrackingSeparateLeadingSubleading&&trackCorrections.size()<3) {
         cout << "Fatal Error: Not enough correction tables to do separate leading/subleading jet tracking correction" << endl;
         exit(1);
      }
   }
}

void HiForest::CheckTree(TTree *t,const char *title)
{
  int entries = t->GetEntries();
  if (nEntries==0) nEntries = entries;
  cout <<title<<": "<<entries<<" entries loaded.";
  if (entries != nEntries) {
    cout <<" Inconsistent number of entries!!"<<endl;
  } else {
    cout <<endl;
  }
}


void HiForest::CheckArraySizes(){

  vector<int> trackOverflow(0);
  vector<int> objectOverflow(0);

  for(int ie = 0; ie < nEntries; ++ie){
    GetEntry(ie);
    if(track.nTrk > maxEntryTrack) trackOverflow.push_back(track.nTrk);
    if(akPu3PF.nref > maxEntry) objectOverflow.push_back(akPu3PF.nref);
    if(ak3PF.nref > maxEntry) objectOverflow.push_back(ak3PF.nref);
    if(icPu5.nref > maxEntry) objectOverflow.push_back(icPu5.nref);
  }

  for(int i = 0; i < trackOverflow.size(); ++i){
    cout<<trackOverflow[i]<<endl;
  }

  if(trackOverflow.size() == 0) cout<<"Track sizes OK"<<endl;
  else cout<<"tracks crash"<<endl; // TODO : really crash

  for(int i = 0; i < objectOverflow.size(); ++i){
    cout<<objectOverflow[i]<<endl;
  }

  if(objectOverflow.size() == 0) cout<<"Object sizes OK"<<endl;
  else cout<<"objects crash"<<endl; // TODO : really crash

}



void HiForest::PrintStatus()
{
  if (hasHltTree)      CheckTree(hltTree,      "HltTree");
  if (hasSkimTree)     CheckTree(skimTree,     "SkimTree");
  if (hasEvtTree)      CheckTree(evtTree,      "EvtTree");
  if (hasIcPu5JetTree) CheckTree(icPu5jetTree, "IcPu5jetTree");
  if (hasAkPu2JetTree) CheckTree(akPu2jetTree, "AkPu2jetTree");
  if (hasAkPu3JetTree) CheckTree(akPu3jetTree, "AkPu3jetTree");
  if (hasAk3JetTree) CheckTree(ak3jetTree, "Ak3jetTree");
  if (hasAkPu4JetTree) CheckTree(akPu4jetTree, "AkPu4jetTree");
  if (hasAkPu2CaloJetTree) CheckTree(akPu2CaloJetTree, "AkPu2CaloJetTree");
  if (hasAkPu3CaloJetTree) CheckTree(akPu3CaloJetTree, "AkPu3CaloJetTree");
  if (hasAkPu4CaloJetTree) CheckTree(akPu4CaloJetTree, "AkPu4CaloJetTree");
  if (hasTrackTree)    CheckTree(trackTree,    "TrackTree");
  if (hasPixTrackTree) CheckTree(pixtrackTree, "PixTrackTree");
  if (hasPhotonTree)   CheckTree(photonTree,   "PhotonTree");
  if (hasPFTree)       CheckTree(pfTree,   "PFTree");
  if (hasMetTree)      CheckTree(metTree,   "MetTree");
  if (hasTowerTree)    CheckTree(towerTree,    "TowerTree");
  if (hasHbheTree)     CheckTree(hbheTree,     "HbheTree");
  if (hasEbTree)       CheckTree(ebTree,     "EbTree");
  if (hasGenpTree)      CheckTree(genpTree,   "GenpTree");
  if (hasGenParticleTree)      CheckTree(genParticleTree,   "GenParticleTree");

}

void HiForest::SetOutputFile(const char *name)
{
  outf = new TFile(name,"recreate");
  if (hasHltTree)      AddCloneTree(hltTree,      "hltanalysis",        "HltTree");
  if (hasSkimTree)     AddCloneTree(skimTree,     "skimanalysis",       "HltTree");
  if (hasIcPu5JetTree) AddCloneTree(icPu5jetTree, "icPu5JetAnalyzer",   "t");
  if (hasAkPu2JetTree) AddCloneTree(akPu2jetTree, "akPu2PFJetAnalyzer", "t");
  if (hasAkPu3JetTree) AddCloneTree(akPu3jetTree, "akPu3PFJetAnalyzer", "t");
  if (hasAk3JetTree) AddCloneTree(ak3jetTree, "ak3PFJetAnalyzer", "t");
  if (hasAkPu4JetTree) AddCloneTree(akPu4jetTree, "akPu4PFJetAnalyzer", "t");
  if (hasAkPu2CaloJetTree) AddCloneTree(akPu2CaloJetTree, "akPu2CaloJetAnalyzer", "t");
  if (hasAkPu3CaloJetTree) AddCloneTree(akPu3CaloJetTree, "akPu3CaloJetAnalyzer", "t");
  if (hasAkPu4CaloJetTree) AddCloneTree(akPu4CaloJetTree, "akPu4CaloJetAnalyzer", "t");
  if (hasTrackTree)    AddCloneTree(trackTree,    "anaTrack",           "trackTree");
//  if (hasTrackTree)    AddCloneTree(trackTree,    "ppTrack",           "trackTree");
  if (hasPixTrackTree) AddCloneTree(pixtrackTree, "anaPixTrack",        "trackTree");
  if (hasPhotonTree)   AddCloneTree(photonTree,   "multiPhotonAnalyzer",            "photon");
  if (hasPFTree)   AddCloneTree(pfTree,   "pfcandAnalyzer",            "pfTree");
  if (hasEvtTree)      AddCloneTree(evtTree,      "hiEvtAnalyzer",            "HiTree");
  if (hasMetTree)      AddCloneTree(metTree,      "anaMET",            "metTree");
  if (hasTowerTree)    AddCloneTree(towerTree,    "rechitanalyzer",              "tower");
  if (hasHbheTree)     AddCloneTree(hbheTree,     "rechitanalyzer",               "hbhe");
  if (hasEbTree)       AddCloneTree(ebTree,       "rechitanalyzer",               "eb");
  setupOutput = true;
}

void HiForest::AddCloneTree(TTree* t, const char *dirName, const char *treeName)
{
  // Make directory
  outf->cd();
  outf->mkdir(dirName);
  outf->cd(dirName);

  // Add a clone tree to the clone forest
  TTree *tClone = t->CloneTree(0);
  tClone->SetMaxTreeSize(4000000000);
  tClone->SetName(treeName);
  
  cloneForest.push_back(tClone);
}

void HiForest::FillOutput()
{
  if (setupOutput) {
     for (unsigned int i=0; i<cloneForest.size(); i++)
     {
       cloneForest[i]->Fill();
     } 
  } else {
       cout <<"ERROR: Specify an output file by hiForest.SetOutputFile(filename)!"<<endl;
  }
}

// ====================== Event Utilities ========================

bool HiForest::selectEvent(){
  /*
   bool select = skim.phbheReflagNewTimeEnv 
      && 
      skim.phcalTimingFilter 
      && 
      skim.pHBHENoiseFilter 
      && 
      skim.phiEcalRecHitSpikeFilter;
  */
  bool select = skim.pHBHENoiseFilter;
   if(!pp){
      select = select && skim.pcollisionEventSelection;
   }else{
     //      select = select && skim.phfCoincFilter && skim.ppurityFractionFilter;
   }
   return select;
}

TCut HiForest::eventSelection(){
  //   TCut select("skim.phbheReflagNewTimeEnv && skim.phcalTimingFilter && skim.pHBHENoiseFilter && skim.phiEcalRecHitSpikeFilter");
  TCut select("skim.pHBHENoiseFilter");
  if(!pp){
    select = select && "skim.pcollisionEventSelection";
  }else{
    //      select = select && "skim.phfCoincFilter && skim.ppurityFractionFilter";
  }
   return select;
}


void HiForest::GetEnergyScaleTable(char *fileNameTable) {
   TFile* f = new TFile(fileNameTable);
   fEnergyScale[0][1] = (TF1*)f->Get("fit_hscale_r9gt94_1");
   fEnergyScale[0][2] = (TF1*)f->Get("fit_hscale_r9gt94_2");
   fEnergyScale[0][3] = (TF1*)f->Get("fit_hscale_r9gt94_3");
   fEnergyScale[1][1] = (TF1*)f->Get("fit_hscale_r9lt94_1");
   fEnergyScale[1][2] = (TF1*)f->Get("fit_hscale_r9lt94_2");
   fEnergyScale[1][3] = (TF1*)f->Get("fit_hscale_r9lt94_3");
}


// ====================== Track Utilities ========================
#include "TrackUtilities.C"

// ======================= Jet Utilities =========================
#include "JetUtilities.C"

// ====================== Photon Utilities ========================
#include "PhotonUtilities.C"


#endif
