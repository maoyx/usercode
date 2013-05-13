//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jul  8 19:18:37 2012 by ROOT version 5.27/06b
// from TTree HiTree/
// found on file: /d102/yjlee/hiForest2MC/Pythia80_HydjetDrum_mix01_HiForest2_v20.root
//////////////////////////////////////////////////////////
#include "commonSetup.h"
#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>

class Evts {
public :
   Evts(){};
   ~Evts(){};

   // Declaration of leaf types
   Int_t           run;
   Int_t           evt;
   Int_t           lumi;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Int_t           hiBin;
   Float_t         hiHF;
   Float_t         hiHFplus;
   Float_t         hiHFminus;
   Float_t         hiZDC;
   Float_t         hiZDCplus;
   Float_t         hiZDCminus;
   Float_t         hiHFhit;
   Float_t         hiHFhitPlus;
   Float_t         hiHFhitMinus;
   Float_t         hiET;
   Float_t         hiEE;
   Float_t         hiEB;
   Float_t         hiEEplus;
   Float_t         hiEEminus;
   Int_t           hiNpix;
   Int_t           hiNpixelTracks;
   Int_t           hiNtracks;
   Int_t           hiNtracksPtCut;
   Int_t           hiNtracksEtaCut;
   Int_t           hiNtracksEtaPtCut;
   Int_t           hiNevtPlane;
   Float_t         hiEvtPlanes[126];   //[hiNevtPlane]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_hiBin;   //!
   TBranch        *b_hiHF;   //!
   TBranch        *b_hiHFplus;   //!
   TBranch        *b_hiHFminus;   //!
   TBranch        *b_hiZDC;   //!
   TBranch        *b_hiZDCplus;   //!
   TBranch        *b_hiZDCminus;   //!
   TBranch        *b_hiHFhit;   //!
   TBranch        *b_hiHFhitPlus;   //!
   TBranch        *b_hiHFhitMinus;   //!
   TBranch        *b_hiET;   //!
   TBranch        *b_hiEE;   //!
   TBranch        *b_hiEB;   //!
   TBranch        *b_hiEEplus;   //!
   TBranch        *b_hiEEminus;   //!
   TBranch        *b_hiNpix;   //!
   TBranch        *b_hiNpixelTracks;   //!
   TBranch        *b_hiNtracks;   //!
   TBranch        *b_hiNtracksPtCut;   //!
   TBranch        *b_hiNtracksEtaCut;   //!
   TBranch        *b_hiNtracksEtaPtCut;   //!
   TBranch        *b_hiNevtPlane;   //!
   TBranch        *b_hiEvtPlanes;   //!

};


void setupEvtTree(TTree *t,Evts &tEvts,bool doCheck = 1)
{
   // Set branch addresses and branch pointers
   if (t->GetBranch("run")) t->SetBranchAddress("run", &tEvts.run, &tEvts.b_run);
   if (t->GetBranch("evt")) t->SetBranchAddress("evt", &tEvts.evt, &tEvts.b_evt);
   if (t->GetBranch("lumi")) t->SetBranchAddress("lumi", &tEvts.lumi, &tEvts.b_lumi);
   if (t->GetBranch("vx")) t->SetBranchAddress("vx", &tEvts.vx, &tEvts.b_vx);
   if (t->GetBranch("vy")) t->SetBranchAddress("vy", &tEvts.vy, &tEvts.b_vy);
   if (t->GetBranch("vz")) t->SetBranchAddress("vz", &tEvts.vz, &tEvts.b_vz);
   if (t->GetBranch("hiBin")) t->SetBranchAddress("hiBin", &tEvts.hiBin, &tEvts.b_hiBin);
   if (t->GetBranch("hiHF")) t->SetBranchAddress("hiHF", &tEvts.hiHF, &tEvts.b_hiHF);
   if (t->GetBranch("hiHFplus")) t->SetBranchAddress("hiHFplus", &tEvts.hiHFplus, &tEvts.b_hiHFplus);
   if (t->GetBranch("hiHFminus")) t->SetBranchAddress("hiHFminus", &tEvts.hiHFminus, &tEvts.b_hiHFminus);
   if (t->GetBranch("hiZDC")) t->SetBranchAddress("hiZDC", &tEvts.hiZDC, &tEvts.b_hiZDC);
   if (t->GetBranch("hiZDCplus")) t->SetBranchAddress("hiZDCplus", &tEvts.hiZDCplus, &tEvts.b_hiZDCplus);
   if (t->GetBranch("hiZDCminus")) t->SetBranchAddress("hiZDCminus", &tEvts.hiZDCminus, &tEvts.b_hiZDCminus);
   if (t->GetBranch("hiHFhit")) t->SetBranchAddress("hiHFhit", &tEvts.hiHFhit, &tEvts.b_hiHFhit);
   if (t->GetBranch("hiHFhitPlus")) t->SetBranchAddress("hiHFhitPlus", &tEvts.hiHFhitPlus, &tEvts.b_hiHFhitPlus);
   if (t->GetBranch("hiHFhitMinus")) t->SetBranchAddress("hiHFhitMinus", &tEvts.hiHFhitMinus, &tEvts.b_hiHFhitMinus);
   if (t->GetBranch("hiET")) t->SetBranchAddress("hiET", &tEvts.hiET, &tEvts.b_hiET);
   if (t->GetBranch("hiEE")) t->SetBranchAddress("hiEE", &tEvts.hiEE, &tEvts.b_hiEE);
   if (t->GetBranch("hiEB")) t->SetBranchAddress("hiEB", &tEvts.hiEB, &tEvts.b_hiEB);
   if (t->GetBranch("hiEEplus")) t->SetBranchAddress("hiEEplus", &tEvts.hiEEplus, &tEvts.b_hiEEplus);
   if (t->GetBranch("hiEEminus")) t->SetBranchAddress("hiEEminus", &tEvts.hiEEminus, &tEvts.b_hiEEminus);
   if (t->GetBranch("hiNpix")) t->SetBranchAddress("hiNpix", &tEvts.hiNpix, &tEvts.b_hiNpix);
   if (t->GetBranch("hiNpixelTracks")) t->SetBranchAddress("hiNpixelTracks", &tEvts.hiNpixelTracks, &tEvts.b_hiNpixelTracks);
   if (t->GetBranch("hiNtracks")) t->SetBranchAddress("hiNtracks", &tEvts.hiNtracks, &tEvts.b_hiNtracks);
   if (t->GetBranch("hiNtracksPtCut")) t->SetBranchAddress("hiNtracksPtCut", &tEvts.hiNtracksPtCut, &tEvts.b_hiNtracksPtCut);
   if (t->GetBranch("hiNtracksEtaCut")) t->SetBranchAddress("hiNtracksEtaCut", &tEvts.hiNtracksEtaCut, &tEvts.b_hiNtracksEtaCut);
   if (t->GetBranch("hiNtracksEtaPtCut")) t->SetBranchAddress("hiNtracksEtaPtCut", &tEvts.hiNtracksEtaPtCut, &tEvts.b_hiNtracksEtaPtCut);
   if (t->GetBranch("hiNevtPlane")) t->SetBranchAddress("hiNevtPlane", &tEvts.hiNevtPlane, &tEvts.b_hiNevtPlane);
   if (t->GetBranch("hiEvtPlanes")) t->SetBranchAddress("hiEvtPlanes", tEvts.hiEvtPlanes, &tEvts.b_hiEvtPlanes);
   if (doCheck) {
      if (t->GetMaximum("hiNevtPlane")>126) { cout <<"FATAL ERROR: Arrary size of hiNevtPlane too small!!!  "<<t->GetMaximum("hiNevtPlane")<<endl; exit(0);
 }   }
}

