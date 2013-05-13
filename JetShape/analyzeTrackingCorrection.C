#include "../HiForest_HEAD/hiForest.h"
#include <TFile.h>
#include <TH1D.h>
#include <TNtuple.h>
#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include "TChain.h"
#include "TMath.h"
#include "analyzeTrackingCorrection.h"
using namespace std;
using namespace TMath;

void analyzeTrackingCorrection(
   TString jetAlgo = "akPu3PF",
   TString trkCol = "anaTrack",
   TString inname="/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/Dijet${pthat}_HydjetDrum_v27_mergedV1.root",
   TString outname="trkcorr/test/TrkCorrtest_hy18dj100.root",
   bool isPP=false,
   double samplePtHat=0,
   double ptHatMax=9999,
   double sampleWeight = 1,
   double cutPtTrk=0.5,
   double vzMax = 15,
   int maxEntries = -1,
   double leadingJetPtMin=-1,
   double subleadingJetPtMin=-1,
   double sigDPhi=-1,
   bool genJetMode=false
) {
   TH1::SetDefaultSumw2();
   TString tag=Form("%s_%.0f_%.0f_%.0f_genJetMode%d",jetAlgo.Data(),leadingJetPtMin,subleadingJetPtMin,sigDPhi*1000,genJetMode);
   outname.ReplaceAll(".root",Form("_%s.root",tag.Data()));
   cout << "Input: " << inname << endl;
   cout << "Sample pthat = " << samplePtHat << " ptHatMax = " << ptHatMax << endl;
   cout << "Track pt min = " << cutPtTrk << endl;
   cout << "skim: leading Jet > " << leadingJetPtMin << " subleading > " << subleadingJetPtMin << " dphi > " << sigDPhi  << endl;
   cout << "Genjet mode: " << genJetMode << endl;
   cout << "Output: " << outname << endl;

   ///////////////////////////////////////////////////
   // Setup Analysis
   ///////////////////////////////////////////////////
   double cutjetPt = 40;
   double cutjetEta = 2;
   double cutEtaTrk = 2.4;
   
   // Define the input file and HiForest
//    HiForest * c = new HiForest(inname,"mergedTrack");
   HiForest * c = new HiForest(inname,trkCol.Data());
   c->doTrackCorrections = false;

   // intialize jet variables
   Jets * anajet = 0;
   if (jetAlgo=="akPu3PF") anajet = &(c->akPu3PF);
   else if (jetAlgo=="akPu3Calo") anajet = &(c->akPu3Calo);
   else {
      cout << "Fatal Error: jetalgo " << jetAlgo << " not defined" << endl;
      exit(1);
   }
   
   // Output file
   cout << "Output: " << outname << endl;
   TFile *output = new TFile(outname,"recreate");

   EvtSel evt;
   DiJet gj;
//    TrkReso trkRes1(1.6);

   ///////////////////////////////////////////////////
   // Book Histograms
   ///////////////////////////////////////////////////
   // Tracking Corrections
   TrkCorrHisAna effMergedGeneral("Forest2_MergedGeneral",output,40,isPP);
   effMergedGeneral.DeclareHistograms();

   TrkCorrHisAna effMergedGeneral_j1("Forest2_MergedGeneral_j1",output,40,isPP);
   effMergedGeneral_j1.DeclareHistograms();

   TrkCorrHisAna effMergedGeneral_j2("Forest2_MergedGeneral_j2",output,40,isPP);
   effMergedGeneral_j2.DeclareHistograms();

   TrkCorrHisAna effMergedGeneral_jetfine("Forest2_MergedGeneral_jetfine",output,40,isPP);
   effMergedGeneral_jetfine.jetBins.clear();
   float fineJetBins[10] = {0,40,80,100,120,160,200,250,500,1000};
   effMergedGeneral_jetfine.jetBins.insert(effMergedGeneral_jetfine.jetBins.begin(),fineJetBins,fineJetBins+9);
   effMergedGeneral_jetfine.DeclareHistograms();

   // basics
   output->cd();
   TH1D * hCent = new TH1D("hCent","",40,0,40);
   TH1D * hVz = new TH1D("hVz","",60,-30,30);
   TH1D * hPtHatBeforeSel = new TH1D("hPtHatBeforeSel","",200,0,1000);
   TH1D * hPtHat = new TH1D("hPtHat","",200,0,1000);
   TH2D * hJetPt2D = new TH2D("hJetPt2D","",100,0,500,100,0,500);
   TH1D * hJDPhi = new TH1D("hJDPhi","",40,0,Pi());
   TH1D * hAj = new TH1D("hAj","",32,0,0.8);

   // In separate centrality bins
   vector<TH1D*> vhCent;
   vector<TH1D*> vhPtHat;
   vector<TH1D*> vhPtHatBeforeSel;
   vector<TH2D*> vhJetPt2D;
   for (int ib=0; ib<effMergedGeneral.centBins.size()-1; ++ib) {
      vhCent.push_back(new TH1D(Form("hCent_c%d",ib),
      Form("cbin%dto%d;Centrality Bin;",effMergedGeneral.centBins[ib],effMergedGeneral.centBins[ib+1]),40,0,40));
      vhPtHat.push_back(new TH1D(Form("hPtHat_c%d",ib),"",200,0,1000));
      vhPtHatBeforeSel.push_back(new TH1D(Form("hPtHatBeforeSel_c%d",ib),"",200,0,1000));
      vhJetPt2D.push_back(new TH2D(Form("hJetPt2D_c%d",ib),"",100,0,500,100,0,500));
   }

   ///////////////////////////////////////////////////
   // Main loop
   ///////////////////////////////////////////////////
   if (maxEntries<0) maxEntries = c->GetEntries();
   for (int i=0;i<maxEntries;i++) {
      c->GetEntry(i);
      evt.cBin = c->evt.hiBin;
      evt.evtPlane = c->evt.hiEvtPlanes[21];
      evt.nJ = anajet->nref;
      evt.nT = c->track.nTrk;
      evt.offlSel = (c->skim.pcollisionEventSelection > 0);
      if (!c->hasSkimTree) evt.offlSel = (c->evt.hiNtracks>0 && c->evt.hiHFplus>=4 && c->evt.hiHFminus>=4);
      evt.vz = c->track.vz[1];
      evt.sampleWeight = sampleWeight/maxEntries; // for different mc sample, 1 for data
      evt.pthat = anajet->pthat;
      evt.samplePtHat = samplePtHat;

      if (i%1000==0) cout <<i<<" / "<< maxEntries << " pthat: " << evt.pthat << " cBin: " << evt.cBin << " nT: " << evt.nT <<endl;

      // initialize
      int leadingIndex=-1,genLeadingIndex=-1;
      int awayIndex=-1,genAwayIndex=-1;
      gj.clear();
      
      ///////////////////////////////////////////////////////
      // Find Leading jets
      ///////////////////////////////////////////////////////
      for (int j=0;j<anajet->nref;j++) {
         if (anajet->jtpt[j]<cutjetPt) continue;
         if (fabs(anajet->jteta[j])>cutjetEta) continue;
         if (anajet->jtpt[j] > gj.pt1) {
            gj.pt1 = anajet->jtpt[j];
            leadingIndex = j;
         }
      }

      // Found a leading jet which passed basic quality cut!
      if (leadingIndex!=-1) {
         // set leading jet
         gj.pt1raw=anajet->rawpt[leadingIndex];
         gj.eta1=anajet->jteta[leadingIndex];
         gj.phi1=anajet->jtphi[leadingIndex];
         gj.ref1pt = anajet->refpt[leadingIndex];
         gj.ref1eta = anajet->refeta[leadingIndex];
         gj.ref1phi = anajet->refphi[leadingIndex];
         gj.ref1partonpt = anajet->refparton_pt[leadingIndex];
         gj.ref1partonflavor = anajet->refparton_flavor[leadingIndex];
         
         // Loop over jet tree to find a away side leading jet
         gj.nJet=0;
         for (int j=0;j<anajet->nref;j++) {
            if (anajet->jtpt[j]<cutjetPt) continue;
            if (fabs(anajet->jteta[j])>cutjetEta) continue;
            gj.inclJetPt[gj.nJet] = anajet->jtpt[j];
            gj.inclJetEta[gj.nJet] = anajet->jteta[j];
            gj.inclJetPhi[gj.nJet] = anajet->jtphi[j];
            gj.inclJetRefPt[gj.nJet] = anajet->refpt[j];
            gj.inclJetRefPartonPt[gj.nJet] = anajet->refparton_pt[j];
            if (j!=leadingIndex&&anajet->jtpt[j] > gj.pt2) {
               gj.pt2 = anajet->jtpt[j];
               awayIndex = j;
            }
            ++gj.nJet;
         }

         if (awayIndex !=-1) { // Found an away jet!
            gj.pt2raw = anajet->rawpt[awayIndex];
            gj.eta2   = anajet->jteta[awayIndex];
            gj.phi2   = anajet->jtphi[awayIndex];
            gj.deta   = gj.eta2 - gj.eta1;
            gj.dphi   = deltaPhi(gj.phi2,gj.phi1);
            gj.Aj     = (gj.pt1-gj.pt2)/(gj.pt1+gj.pt2);
            gj.ref2pt = anajet->refpt[awayIndex];
            gj.ref2eta = anajet->refeta[awayIndex];
            gj.ref2phi = anajet->refphi[awayIndex];
            gj.ref2partonpt = anajet->refparton_pt[awayIndex];
            gj.ref2partonflavor = anajet->refparton_flavor[awayIndex];
         }         
      } // end of if leadingIndex
      
      ////////////////////////////////////////
      // Find Leading Genjets
      ////////////////////////////////////////
      // Loop over genjets to look for leading genjet candidate in the event
      for (int j=0;j<anajet->ngen;j++) {
         if (anajet->genpt[j]<cutjetPt) continue;
         if (fabs(anajet->geneta[j])>cutjetEta) continue;
         if (anajet->genpt[j] > gj.genjetpt1) {
            gj.genjetpt1=anajet->genpt[j];
            genLeadingIndex=j;
         }
      }

      if (genLeadingIndex!=-1) {
         gj.genjeteta1=anajet->geneta[genLeadingIndex];
         gj.genjetphi1=anajet->genphi[genLeadingIndex];
      }
      // subleading
      gj.nGenJet=0;
      for (int j=0;j<anajet->ngen;j++) {
         if (anajet->genpt[j]<cutjetPt) continue;
         if (fabs(anajet->geneta[j])>cutjetEta) continue;
         gj.inclGenJetPt[gj.nGenJet] = anajet->genpt[j];
         gj.inclGenJetEta[gj.nGenJet] = anajet->geneta[j];
         gj.inclGenJetPhi[gj.nGenJet] = anajet->genphi[j];
         if (j!=genLeadingIndex && anajet->genpt[j]>gj.genjetpt2) {
            gj.genjetpt2=anajet->genpt[j];
            gj.genjeteta2=anajet->geneta[j];
            gj.genjetphi2=anajet->genphi[j];
         }
         ++gj.nGenJet;
      }
      if (genJetMode) {
         gj.pt1 = gj.genjetpt1;
         gj.eta1 = gj.genjeteta1;
         gj.phi1 = gj.genjetphi1;
         gj.pt2 = gj.genjetpt2;
         gj.eta2 = gj.genjeteta2;
         gj.phi2 = gj.genjetphi2;
         gj.dphi   = deltaPhi(gj.phi2,gj.phi1);
      }

      ///////////////////////////////////////////////////////
      // Skim
      ///////////////////////////////////////////////////////
      if (!evt.offlSel) continue;
      if (samplePtHat>0 && evt.pthat>=ptHatMax) continue;
      if (vzMax>0 && fabs(evt.vz)>vzMax) continue;
      // protection against high pt jet from background event
      if (leadingIndex>=0&&anajet->subid[leadingIndex]>0) continue;
      hPtHatBeforeSel->Fill(evt.pthat);
      for (int ib=0; ib<effMergedGeneral.centBins.size(); ++ib) {
         if(evt.cBin>=effMergedGeneral.centBins[ib] && evt.cBin<effMergedGeneral.centBins[ib+1]){
            vhPtHatBeforeSel[ib]->Fill(evt.pthat);
         }
      }
      // ensure jet distribution unbiased by pthat turn on
      if (leadingJetPtMin>0&&gj.pt1<leadingJetPtMin) continue;
      if (subleadingJetPtMin>0&&gj.pt2<subleadingJetPtMin) continue;
      if (sigDPhi>0&&fabs(gj.dphi)<sigDPhi) continue;
      
      // Fill Baisc Event info
      hCent->Fill(evt.cBin);
      hVz->Fill(evt.vz);
      hPtHat->Fill(evt.pthat);
      hJetPt2D->Fill(gj.pt1,gj.pt2);
      hJDPhi->Fill(fabs(gj.dphi));
      hAj->Fill(gj.Aj);
      for (int ib=0; ib<effMergedGeneral.centBins.size(); ++ib) {
         if(evt.cBin>=effMergedGeneral.centBins[ib] && evt.cBin<effMergedGeneral.centBins[ib+1]){
            vhCent[ib]->Fill(evt.cBin);
            vhPtHat[ib]->Fill(evt.pthat);
            vhJetPt2D[ib]->Fill(gj.pt1,gj.pt2);
         }
      }

      ///////////////////////////////////////////////////////
      // Tracks
      ///////////////////////////////////////////////////////
      // Full Tracks, Pixel Tracks
      Tracks * anaTrks[1] = {&(c->track)};
      for (int it=0; it<anaTrks[0]->nTrk; ++it) {
         // Kinematic Selection
         if (anaTrks[0]->trkPt[it] < cutPtTrk) continue;
         if (fabs(anaTrks[0]->trkEta[it]) > cutEtaTrk) continue;
         if (anaTrks[0]->trkAlgo[it]>=4&&!anaTrks[0]->highPurity[it]) continue; // quality selection

         RecTrack_t r;
         r.ptr = anaTrks[0]->trkPt[it];
         r.etar = anaTrks[0]->trkEta[it];
         r.phir = anaTrks[0]->trkPhi[it];
         r.algo = anaTrks[0]->trkAlgo[it];
         r.nsim = !anaTrks[0]->trkFake[it];
         r.status = 1; // for now correct all tracks
         float dr1 = deltaR(r.etar,r.phir,gj.eta1,gj.phi1);
         float dr2 = deltaR(r.etar,r.phir,gj.eta2,gj.phi2);
         r.jet  = 0;
         r.jeta = -99;
         r.jdr  = -99;
         if (samplePtHat>0) {
           if (dr1<0.5&&gj.pt1>=40) {
              r.jet = gj.pt1;
              r.jeta = gj.eta1;
              r.jdr = dr1;
              effMergedGeneral_j1.FillRecHistograms(evt,gj,r);
           } else if (dr2<0.5&&gj.pt2>=40) {
              r.jet = gj.pt2;
              r.jeta = gj.eta2;
              r.jdr = dr2;
              effMergedGeneral_j2.FillRecHistograms(evt,gj,r);
           }
         }

         // Fill
         effMergedGeneral.FillRecHistograms(evt,gj,r);
         effMergedGeneral_jetfine.FillRecHistograms(evt,gj,r);
//          effMergedGeneral_trkPhi.FillRecHistograms(evt,gj,r);
//          if (r.jet>=120) effMergedGeneral_trkPhi_jet_120to999.FillRecHistograms(evt,gj,r);
//          else if (r.jet>=cutjetPt) effMergedGeneral_trkPhi_jet_50to120.FillRecHistograms(evt,gj,r);
//          else effMergedGeneral_trkPhi_noJet.FillRecHistograms(evt,gj,r);
//          effMergedGeneral_SimSmear.FillRecHistograms(evt,gj,r);
      }

      ///////////////////////////////////////////////////////
      // SimTracks
      ///////////////////////////////////////////////////////
      for (int ip=0; ip<anaTrks[0]->nParticle; ++ip) {
         if (anaTrks[0]->pPt[ip] < cutPtTrk) continue;
         if (fabs(anaTrks[0]->pEta[ip]) > cutEtaTrk) continue;

         SimTrack_t s;
         s.pts = anaTrks[0]->pPt[ip];
         s.etas = anaTrks[0]->pEta[ip];
         s.phis = anaTrks[0]->pPhi[ip];
         s.status = 1; // for now assume all simtracks as signal tracks
         s.nrec = (anaTrks[0]->pNRec[ip])*(int)(anaTrks[0]->mtrkAlgo[ip]<4||anaTrks[0]->mtrkQual[ip]>0);
//          if (s.pts<1.5) s.nrec = anaTrks[0]->pNRec[ip];
//          else s.nrec = (anaTrks[0]->pNRec[ip])*(anaTrks[0]->mtrkQual[ip]);
         s.acc = anaTrks[0]->pAcc[ip];
         s.jet = gj.pt1;
         s.jeta = gj.eta1;
         float dr1 = deltaR(s.etas,s.phis,gj.eta1,gj.phi1);
         float dr2 = deltaR(s.etas,s.phis,gj.eta2,gj.phi2);
         s.jet = 0;
         s.jeta = -99;
         s.jdr = -99;
         if (samplePtHat>0) {
           if (dr1<0.5&&gj.pt1>=40) {
              s.jet = gj.pt1;
              s.jeta = gj.eta1;
              s.jdr = dr1;
              effMergedGeneral_j1.FillSimHistograms(evt,gj,s);
           } else if (dr2<0.5&&gj.pt2>=40) {
              s.jet = gj.pt2;
              s.jeta = gj.eta2;
              s.jdr = dr2;
              effMergedGeneral_j2.FillSimHistograms(evt,gj,s);
           }
         }
         // Fill
         effMergedGeneral.FillSimHistograms(evt,gj,s);
         effMergedGeneral_jetfine.FillSimHistograms(evt,gj,s);
//          effMergedGeneral_trkPhi.FillSimHistograms(evt,gj,s);
//          if (s.jet>=120) effMergedGeneral_trkPhi_jet_120to999.FillSimHistograms(evt,gj,s);
//          else if (s.jet>=cutjetPt) effMergedGeneral_trkPhi_jet_50to120.FillSimHistograms(evt,gj,s);
//          else effMergedGeneral_trkPhi_noJet.FillSimHistograms(evt,gj,s);
         // Sim Smearing
//          if (s.pts>=1&&s.pts<2) {
//             float trksm = trkRes1.GetSmear(s.pts);
//             s.pts*=trksm;
//          }
//          effMergedGeneral_SimSmear.FillSimHistograms(evt,gj,s);
      }
      // All done
   }

//    effMergedGeneral.WriteHistograms();
   output->Write();
   output->Close();
   delete c;
}

