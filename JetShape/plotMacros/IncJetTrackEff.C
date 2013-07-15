/*
 *  Display.C
 *  
 *
 *  Created by ymao on 12/28/11.
 *  Copyright 2011 Central China Normal University  (CN). All rights reserved.
 *
 */
#include <iostream>
#include "TCanvas.h"
#include "TError.h"
#include "TPad.h"
#include "TString.h"
#include "TRandom.h"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include <TF1.h>
#include <TMath.h>
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TString.h"
#include "TCut.h"
#include <TMath.h>
#include <TSystem.h>
#include "THStack.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include <iostream>
using namespace std;

//define the kinematics cuts for the analysis
//define the kinematics cuts for the analysis
const double conesize = 0.3 ;
const double deltacone = 0.05 ;

//const double radius = 0.025;

//const int nrbin = 10 ;
//double rbin[nrbin+1]={0.0,0.05,0.10, 0.15, 0.20,0.25,0.30,0.35,0.40,0.45,0.50};

const int nrbin = 6 ;
double rbin[nrbin+1]={0.0,0.05,0.10, 0.15, 0.20,0.25,0.30};

const double etacut = 2.0 ;
const double dphicut = TMath::TwoPi()/3. ;
const double leadingjetcut = 100. ;
const double subleadjetcut = 40. ;
const double trackcut = 1.0;

const Double_t etalimit = 0.3 ;

bool SavePlot = kFALSE ;
bool SaveFile = kFALSE ;

const bool IsMC=kTRUE ;

bool DoSmear=kFALSE ;
bool DoPtWeight=kFALSE;

const int nptbin = 1 ;
//const double pt[]={100., 120., 140., 160.,200., 300., 500.};
const double pt[]={50., 500.};
//const double trkpt[]={1.,2., 4., 8., 16., 32., 300.};
//const double trkpt[]={1.,1.5,2., 2.5, 4.,8., 16., 32., 500.};
const double trkpt[]={1.,300.};

const Double_t TrkBin[]={0, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.5, 3.0, 4., 5., 6., 8.,10.,12.,15.,20.,25.,30.,45.,60.,80.,100., 120.,150., 180.,300.,500.};
const int nTrkBin = sizeof(TrkBin)/sizeof(Double_t)-1 ;

const Double_t pbpbLumi = 150 ;
const Double_t ppLumi = 231 ;

enum Option_t {kDist, kRatio} ; 
enum Display_t {kPt, kRadius} ;
enum Data_t {kPP, kHI} ;
enum Algo_t {kPF, kCalo} ;
enum Origin_t {kJet, kJetBkg, kInclusive} ;
//---------------------------------------------------------------------

void IncJetTrackEff(Option_t opt=kDist, Origin_t org =kJet, Display_t dis=kRadius, Data_t data = kPP, Algo_t algo = kPF, Int_t trackbin = 0){
    
    TString plotsdir = "/Users/ymao/group/CMS/plots/finalJS/CorrectDR";
    TDatime d ;
    int date = d.GetDate();    
  //  const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/FullStatPbPb" ;      
if(IsMC)
 //   const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/CorrectdR/NoMatchGen/NoTrkAlgoCut" ;
    const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/pPb" ;
else
    const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/CorrectdR" ;
    TString coll ;
    if(data==kPP) coll="PP";
    else coll = "HI" ;
    TString type ;
    if(IsMC) type="MC";
    else type = "DATA" ;
    TString met ;
    if(algo==kPF) met="PF";
    else met = "Calo" ;
    TString dep ;
    if(dis==kPt) dep = "Pt";
    else dep = "Radius" ;
    TString trig ;
    if(org==kJet) trig = "Jet";
    else trig = "Bkg" ;
    TString lev ;
    if(DoPtWeight)lev="wt";
    else lev="";
    TString effTab = "HistIterTrkCorrtest";  //"HistIterTrkCorrv14fXSecMBEff"; // "HistIterTrkCorrv14fXSec"; // "HistIterTrkCorrtest" ; //
    
    if(data==kPP) {
        ////if it is pp, no centrality bins, only one
        const int nbin = 1 ;
        const int centr[] ={0,100};  
    }
    else {
        // for HI centrality bin
        const int nbin = 4 ;
        const int centr[] ={0,10,30,50,100}; 
    }
    
    TString fileName ;
    TFile * f ;
    TString CorrfileName ;
    TFile * Corrf ;
    
    TH2F * recopt[nbin][nptbin];
    TH2F * genpt[nbin][nptbin];
    TH2F * ppgenpt[nbin][nptbin];
    TH2F * matchpt[nbin][nptbin];
    TH2F * corrpt[nbin][nptbin];

    TH2F * sumjetptbinrecopt[nbin];
    TH2F * sumjetptbingenpt[nbin];
    TH2F * ppsumjetptbingenpt[nbin];
    TH2F * sumjetptbinmatchpt[nbin];
    TH2F * sumjetptbincorrpt[nbin];
    
    TH1F * track[nbin];
    TH1F * ppgentrack[nbin];
    TH1F * gentrack[nbin];
    TH1F * matchtrack[nbin]];
    TH1F * corrtrack[nbin];
    
    TH1F * genratio[nbin];
    TH1F * recoEff[nbin]; //matchtrack divide gen track
//    TH1F * matchEff[nbin][nptbin]; //match track divide reco track
//
//    TH2F * bkgrecopt[nbin][nptbin];
//    TH2F * bkggenpt[nbin][nptbin];
//    TH2F * bkgmatchpt[nbin][nptbin];
//    TH2F * bkgcorrpt[nbin][nptbin];
//    
//    TH1F * bkgcorrtrack[nbin][nptbin];
//    TH1F * bkgtrack[nbin][nptbin];
//    TH1F * genbkgtrack[nbin][nptbin];
//    TH1F * matchbkgtrack[nbin][nptbin];
//    
//    TH1F * bkgrecoEff[nbin][nptbin]; //matchtrack divide gen track
//    TH1F * bkgmatchEff[nbin][nptbin]; //match track divide reco track
//
    TH1F * corrClosure[nbin]; //matchtrack divide gen track
//    TH1F * bkgClosure[nbin][nptbin]; //match track divide reco track

    TH1::SetDefaultSumw2();
    
    TString canv_name = "c1";
    //    canv_name += parti;
    //    canv_name += cent;
    if(data==kHI){
        const Double_t kw = 1000;
        const Double_t kh = 560;
        c1 = new TCanvas(canv_name," ",10,10,kw,kh);
        //    makeMultiPanelCanvas(c1,4,1,0.0,0.0,0.2,0.2,0.02);
        if(IsMC) makeMultiPanelCanvas(c1,nbin,2,0.0,0.0,0.2,0.2,0.02);
        else makeMultiPanelCanvas(c1,nbin,1,0.0,0.0,0.2,0.2,0.02);
    }
    else {
        const Double_t kw = 800;
        const Double_t kh = 800;
        c1 = new TCanvas(canv_name," ",10,10,kw,kh);
        //    makeMultiPanelCanvas(c1,4,1,0.0,0.0,0.2,0.2,0.02);
        if(IsMC) makeMultiPanelCanvas(c1,nbin,2,0.0,0.0,0.1,0.2,0.05);  
        else makeMultiPanelCanvas(c1,nbin,1,0.0,0.0,0.1,0.2,0.05);  
    }
    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);
    
    
    TLegend *t1=new TLegend(0.35,0.7,0.6,0.92);
    t1->SetFillColor(0);
    t1->SetBorderSize(0);
    t1->SetFillStyle(0);
    t1->SetTextFont(63);
    t1->SetTextSize(17);
    TLegend *t2=new TLegend(0.25,0.55,0.6,0.7);
    t2->SetFillColor(0);
    t2->SetBorderSize(0);
    t2->SetFillStyle(0);
    t2->SetTextFont(63);
    t2->SetTextSize(17);
 
    // open the data file
    switch(data){
        case kPP:
            if(IsMC){
//                fileName=Form("mergedCSdiff_MCPP2013_Ak3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_ppHiIterativeTrack_merged.root", trackcut);
//                
//                CorrfileName=Form("mergedCSdiff_MCPP2013_Ak3PFIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_ppHiIterativeTrack_merged.root",effTab.Data(), trackcut);
                
//                fileName=Form("mergedCSdiff_MCPP2011_Ak3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_merged.root", trackcut);
//                
//                CorrfileName=Form("mergedCSdiff_MCPP2011_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_merged.root", effTab.Data(), trackcut);
                fileName=Form("mergedCSdiff_MCPP_Ak3PFIncJetPt50_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin4_HP04_hiforest77_hiSignal5TeV.root", trackcut);
                
                CorrfileName=Form("mergedCSdiff_MCPP_Ak3PFIncJetPt50_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin4_HP04_hiforest77_hiSignal5TeV.root",effTab.Data(), trackcut);

            }
            else {
//                fileName=Form("DATAPP_Ak3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_pp_merged_full.root", trackcut);   
//                
//                CorrfileName=Form("DATAPP_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_pp_merged_full.root", effTab.Data(), trackcut);    
                fileName=Form("DATAPP_Ak3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_JetDB_forestv78.root", trackcut);
                
                CorrfileName=Form("DATAPP_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_JetDB_forestv78.root", effTab.Data(), trackcut);
            }
            break ;
        case kHI:
            if(IsMC){
                fileName=Form("mergedCSdiff_MCHI_AkPu3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_Dijet_HydjetDrum_v27v28.root", trackcut);
                CorrfileName=Form("mergedCSdiff_MCHI_AkPu3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_Dijet_HydjetDrum_v27v28.root", effTab.Data(), trackcut);
            }
            else {
                fileName=Form("DATAHI_AkPu3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_promptskim-hihighpt-hltjet80-pt90-v20.root", trackcut);
                CorrfileName=Form("DATAHI_AkPu3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin3_promptskim-hihighpt-hltjet80-pt90-v20.root", effTab.Data(), trackcut, nbin);
            }
            break ;
    }            
    f = TFile::Open(Form("%s/%s", kHomeDir, fileName.Data()), "readonly");
 //   f = TFile::Open(Form("/Users/ymao/group/CMS/anaOutputs/JetShape/CorrectdR/%s", fileName.Data()), "readonly");
//    if(IsMC)Corrf = TFile::Open(Form("%s/NoMatchGen/%s", kHomeDir, CorrfileName.Data()), "readonly");
//    else
    if(data==kPP)
 //       Corrf = TFile::Open(Form("%s/MyselfTrkEff/%s", kHomeDir, CorrfileName.Data()), "readonly");
        Corrf = TFile::Open(Form("%s/%s", kHomeDir, CorrfileName.Data()), "readonly");
    else
       Corrf = TFile::Open(Form("%s/%s", kHomeDir, CorrfileName.Data()), "readonly");
    for(int ibin = 0 ; ibin <nbin; ibin++){

        if(org==kInclusive){
            track[ibin]=(TH1F*)f->Get(Form("inclusivetrackpt_%d-%d%%",centr[ibin],centr[ibin+1]));
            track[ibin]->Sumw2();
            corrtrack[ibin]=(TH1F*)Corrf->Get(Form("inclusivetrackpt_%d-%d%%",centr[ibin],centr[ibin+1]));
            corrtrack[ibin]->Sumw2();
            if(IsMC){
                gentrack[ibin]=(TH1F*)f->Get(Form("inclusivegenpartpt_%d-%d%%",centr[ibin],centr[ibin+1]));
                gentrack[ibin]->Sumw2();
                ppgentrack[ibin]=(TH1F*)Corrf->Get(Form("inclusivegenpartpt_%d-%d%%",centr[ibin],centr[ibin+1]));
                ppgentrack[ibin]->Sumw2();
                matchtrack[ibin]=(TH1F*)f->Get(Form("inclusivegenmatchpt_%d-%d%%",centr[ibin],centr[ibin+1]));
                matchtrack[ibin]->Sumw2();
                for(Int_t i = 1 ; i < =nTrkBin ; i++){
                    ppgentrack[ibin]->SetBinContent(i,ppgentrack[ibin]->GetBinContent(i)/ppgentrack[ibin]->GetBinWidth(i));
                    ppgentrack[ibin]->SetBinError(i,ppgentrack[ibin]->GetBinError(i)/ppgentrack[ibin]->GetBinWidth(i));
                    gentrack[ibin]->SetBinContent(i,gentrack[ibin]->GetBinContent(i)/gentrack[ibin]->GetBinWidth(i));
                    gentrack[ibin]->SetBinError(i,gentrack[ibin]->GetBinError(i)/gentrack[ibin]->GetBinWidth(i));
                    matchtrack[ibin]->SetBinContent(i,matchtrack[ibin]->GetBinContent(i)/matchtrack[ibin]->GetBinWidth(i));
                    matchtrack[ibin]->SetBinError(i,matchtrack[ibin]->GetBinError(i)/matchtrack[ibin]->GetBinWidth(i));
                }
                
            }
            for(Int_t i = 1 ; i < =nTrkBin ; i++){
                track[ibin]->SetBinContent(i,track[ibin]->GetBinContent(i)/track[ibin]->GetBinWidth(i));
                track[ibin]->SetBinError(i,track[ibin]->GetBinError(i)/track[ibin]->GetBinWidth(i));
                corrtrack[ibin]->SetBinContent(i,corrtrack[ibin]->GetBinContent(i)/corrtrack[ibin]->GetBinWidth(i));
                corrtrack[ibin]->SetBinError(i,corrtrack[ibin]->GetBinError(i)/corrtrack[ibin]->GetBinWidth(i));
            }
            
        }  //for inclusive track histograms
        else { //for jet track histograms
            for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
                switch(org){
                    case kJet:
                        recopt[ibin][ipt]=(TH2F*)f->Get(Form("reco%sTrackJetPt%.f_%.f_%d-%d%%",lev.Data(),pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                        recopt[ibin][ipt]->Sumw2();
                        if(IsMC){
                            ppgenpt[ibin][ipt] = (TH2F*)f->Get(Form("%sgenPartJetPt%.f_%.f_%d-%d%%",lev.Data(),pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                            ppgenpt[ibin][ipt]->Sumw2();
                            genpt[ibin][ipt] = (TH2F*)Corrf->Get(Form("%sgenPartJetPt%.f_%.f_%d-%d%%",lev.Data(),pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                            genpt[ibin][ipt]->Sumw2();
                            matchpt[ibin][ipt]=(TH2F*)Corrf->Get(Form("%smatchPartJetPt%.f_%.f_%d-%d%%",lev.Data(),pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                            matchpt[ibin][ipt]->Sumw2();
                        }
                        corrpt[ibin][ipt]=(TH2F*)Corrf->Get(Form("reco%sTrackJetPt%.f_%.f_%d-%d%%",lev.Data(),pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                        corrpt[ibin][ipt]->Sumw2();
                        
                        break ;
                    case kJetBkg:
                        recopt[ibin][ipt]=(TH2F*)f->Get(Form("recobkg%sTrackJetPt%.f_%.f_%d-%d%%",lev.Data(),pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                        recopt[ibin][ipt]->Sumw2();
                        corrpt[ibin][ipt]=(TH2F*)Corrf->Get(Form("recobkg%sTrackJetPt%.f_%.f_%d-%d%%",lev.Data(),pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                        corrpt[ibin][ipt]->Sumw2();
                        if(IsMC){
                            ppgenpt[ibin][ipt] = (TH2F*)f->Get(Form("%sgenbkgPartJetPt%.f_%.f_%d-%d%%",lev.Data(),pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                            ppgenpt[ibin][ipt]->Sumw2();
                            genpt[ibin][ipt] = (TH2F*)Corrf->Get(Form("%sgenbkgPartJetPt%.f_%.f_%d-%d%%",lev.Data(),pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                            genpt[ibin][ipt]->Sumw2();
                            matchpt[ibin][ipt]=(TH2F*)Corrf->Get(Form("%smatchbkgPartJetPt%.f_%.f_%d-%d%%",lev.Data(),pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                            matchpt[ibin][ipt]->Sumw2();
                        }
                        
                        break ;
                }
                if(ipt==0){
                    ppsumjetptbingenpt[ibin]=(TH2F*)ppgenpt[ibin][ipt]->Clone(Form("%s%sppgenPartCent%d-%d%%",trig.Data(),lev.Data(),centr[ibin],centr[ibin+1]));
                    sumjetptbinrecopt[ibin]=(TH2F*)recopt[ibin][ipt]->Clone(Form("%s%srecoTrackCent%d-%d%%",trig.Data(), lev.Data(),centr[ibin],centr[ibin+1]));
                    sumjetptbincorrpt[ibin]=(TH2F*)corrpt[ibin][ipt]->Clone(Form("%s%scorrTrackCent%d-%d%%",trig.Data(),lev.Data(),centr[ibin],centr[ibin+1]));
                    if(IsMC){
                        sumjetptbingenpt[ibin]=(TH2F*)genpt[ibin][ipt]->Clone(Form("%s%sgenPartCent%d-%d%%",trig.Data(),lev.Data(),centr[ibin],centr[ibin+1]));
                        sumjetptbinmatchpt[ibin]=(TH2F*)matchpt[ibin][ipt]->Clone(Form("%s%smatchTrackCent%d-%d%%",trig.Data(),lev.Data(),centr[ibin],centr[ibin+1]));
                    }
                }
                else {
                    sumjetptbinrecopt[ibin]->Add(recopt[ibin][ipt],1);
                    sumjetptbincorrpt[ibin]->Add(corrpt[ibin][ipt],1);
                    if(IsMC){
                        ppsumjetptbingenpt[ibin]->Add(ppgenpt[ibin][ipt],1);
                        sumjetptbingenpt[ibin]->Add(genpt[ibin][ipt],1);
                        sumjetptbinmatchpt[ibin]->Add(matchpt[ibin][ipt],1);
                    }
                    
                }
            } //finish jet pt bin loop
            switch(dis){
                case kPt:
                    track[ibin]=(TH1F*)sumjetptbinrecopt[ibin]->ProjectionX(Form("%sreco%s%sPt_%d-%d%%",lev.Data(),dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]), 0, -1, "");
                    for(Int_t i = 1 ; i < =nTrkBin ; i++){
                        track[ibin]->SetBinContent(i,track[ibin]->GetBinContent(i)/track[ibin]->GetBinWidth(i));
                        track[ibin]->SetBinError(i,track[ibin]->GetBinError(i)/track[ibin]->GetBinWidth(i));
                    }
                    
                    if(IsMC){
                        ppgentrack[ibin]=(TH1F*)ppsumjetptbingenpt[ibin]->ProjectionX(Form("%sppgen%s%sPt_%d-%d%%",lev.Data(),dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]), 0, -1, "");
                        ppgentrack[ibin]->Sumw2();
                        gentrack[ibin]=(TH1F*)sumjetptbingenpt[ibin]->ProjectionX(Form("%sgen%s%sPt_%d-%d%%",lev.Data(),dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]), 0, -1, "");
                        gentrack[ibin]->Sumw2();
                        matchtrack[ibin]=(TH1F*)sumjetptbinmatchpt[ibin]->ProjectionX(Form("%smatch%s%sPt_%d-%d%%",lev.Data(),dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]), 0, -1, "");
                        matchtrack[ibin]->Sumw2();
                        for(Int_t i = 1 ; i < =nTrkBin ; i++){
                            ppgentrack[ibin]->SetBinContent(i,ppgentrack[ibin]->GetBinContent(i)/ppgentrack[ibin]->GetBinWidth(i));
                            ppgentrack[ibin]->SetBinError(i,ppgentrack[ibin]->GetBinError(i)/ppgentrack[ibin]->GetBinWidth(i));
                            gentrack[ibin]->SetBinContent(i,gentrack[ibin]->GetBinContent(i)/gentrack[ibin]->GetBinWidth(i));
                            gentrack[ibin]->SetBinError(i,gentrack[ibin]->GetBinError(i)/gentrack[ibin]->GetBinWidth(i));
                            matchtrack[ibin]->SetBinContent(i,matchtrack[ibin]->GetBinContent(i)/matchtrack[ibin]->GetBinWidth(i));
                            matchtrack[ibin]->SetBinError(i,matchtrack[ibin]->GetBinError(i)/matchtrack[ibin]->GetBinWidth(i));
                        }
                    }
                    corrtrack[ibin]=(TH1F*)sumjetptbincorrpt[ibin]->ProjectionX(Form("%scorr%s%sPt_%d-%d%%",lev.Data(),dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]), 0, -1, "");
                    corrtrack[ibin]->Sumw2();
                    for(Int_t i = 1 ; i < =nTrkBin ; i++){
                        //                    track[ibin]->SetBinContent(i,track[ibin]->GetBinContent(i)/track[ibin]->GetBinWidth(i));
                        //                    track[ibin]->SetBinError(i,track[ibin]->GetBinError(i)/track[ibin]->GetBinWidth(i));
                        corrtrack[ibin]->SetBinContent(i,corrtrack[ibin]->GetBinContent(i)/corrtrack[ibin]->GetBinWidth(i));
                        corrtrack[ibin]->SetBinError(i,corrtrack[ibin]->GetBinError(i)/corrtrack[ibin]->GetBinWidth(i));
                    }
                    
                    break ;
                case kRadius:
                    //                        track[ibin]=(TH1F*)sumjetptbinrecopt[ibin]->ProjectionY(Form("reco%s%sPt_%d-%d%%",dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]), sumjetptbinrecopt[ibin]->GetXaxis()->FindBin(tracklimit), -1, "");
                    //                        gentrack[ibin]=(TH1F*)sumjetptbingenpt[ibin]->ProjectionY(Form("gen%s%sPt_%d-%d%%",dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]), sumjetptbingenpt[ibin]->GetXaxis()->FindBin(tracklimit), -1, "");
                    //                        gentrack[ibin]->Sumw2();
                    //                        matchtrack[ibin]=(TH1F*)sumjetptbinmatchpt[ibin]->ProjectionY(Form("match%s%sPt_%d-%d%%",dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]), sumjetptbinmatchpt[ibin]->GetXaxis()->FindBin(tracklimit), -1, "");
                    //                        matchtrack[ibin]->Sumw2();
                    //
                    //                        corrtrack[ibin]=(TH1F*)sumjetptbincorrpt[ibin]->ProjectionY(Form("corr%s%sPt_%d-%d%%",dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]), sumjetptbincorrpt[ibin]->GetXaxis()->FindBin(tracklimit), -1, "");
                    //                        corrtrack[ibin]->Sumw2();
                    track[ibin]=(TH1F*)sumjetptbinrecopt[ibin]->ProjectionY(Form("%sreco%s%s_%d-%d%%",lev.Data(),dep.Data(),trig.Data(), centr[ibin],centr[ibin+1]), sumjetptbinrecopt[ibin]->GetXaxis()->FindBin(trkpt[trackbin]), sumjetptbinrecopt[ibin]->GetXaxis()->FindBin(trkpt[trackbin+1]), "");
                    if(IsMC){
                        ppgentrack[ibin]=(TH1F*)ppsumjetptbingenpt[ibin]->ProjectionY(Form("%sppgen%s%s_%d-%d%%",lev.Data(),dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]), ppsumjetptbingenpt[ibin]->GetXaxis()->FindBin(trkpt[trackbin]), ppsumjetptbingenpt[ibin]->GetXaxis()->FindBin(trkpt[trackbin+1]), "");
                        ppgentrack[ibin]->Sumw2();
                        gentrack[ibin]=(TH1F*)sumjetptbingenpt[ibin]->ProjectionY(Form("%sgen%s%s_%d-%d%%",lev.Data(),dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]), sumjetptbingenpt[ibin]->GetXaxis()->FindBin(trkpt[trackbin]), sumjetptbingenpt[ibin]->GetXaxis()->FindBin(trkpt[trackbin+1]), "");
                        gentrack[ibin]->Sumw2();
                        matchtrack[ibin]=(TH1F*)sumjetptbinmatchpt[ibin]->ProjectionY(Form("%smatch%s%s_%d-%d%%",lev.Data(),dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]), sumjetptbinmatchpt[ibin]->GetXaxis()->FindBin(trkpt[trackbin]), sumjetptbinmatchpt[ibin]->GetXaxis()->FindBin(trkpt[trackbin+1]), "");
                        matchtrack[ibin]->Sumw2();
                    }
                    corrtrack[ibin]=(TH1F*)sumjetptbincorrpt[ibin]->ProjectionY(Form("%scorr%s%s_%d-%d%%",lev.Data(),dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]), sumjetptbincorrpt[ibin]->GetXaxis()->FindBin(trkpt[trackbin]), sumjetptbincorrpt[ibin]->GetXaxis()->FindBin(trkpt[trackbin+1]), "");
                    corrtrack[ibin]->Sumw2();  
                    
                    break ;
            }
        }
        if(IsMC){
//            gentrack[ibin]->Scale(1./jetpt[ibin]->Integral());
//            ppgentrack[ibin]->Scale(1./jetpt2[ibin]->Integral());
            genratio[ibin]=(TH1F*)gentrack[ibin]->Clone(Form("%sTrkCut%.f_%.fGenRatio%s%s_Cen%d-%d%%",lev.Data(),trkpt[trackbin], trkpt[trackbin+1],dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]));
            genratio[ibin]->Divide(ppgentrack[ibin]);
            recoEff[ibin]=(TH1F*)matchtrack[ibin]->Clone(Form("%sTrkCut%.f_%.fEff%s%s_Cen%d-%d%%",lev.Data(),trkpt[trackbin], trkpt[trackbin+1],dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]));
            
            recoEff[ibin]->Divide(gentrack[ibin]);
            
//            recoEff[ibin]=(TH1F*)track[ibin]->Clone(Form("%sTrkCut%.f_%.fEff%s%s_Cen%d-%d%%",lev.Data(),trkpt[trackbin], trkpt[trackbin+1],dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]));
//            
//            recoEff[ibin]->Divide(matchtrack[ibin]);
            
            corrClosure[ibin]=(TH1F*)corrtrack[ibin]->Clone(Form("%sTrkCut%.f_%.fClosure%s%s_Cen%d-%d%%",lev.Data(),trkpt[trackbin], trkpt[trackbin+1], dep.Data(), trig.Data(), centr[ibin],centr[ibin+1]));
            
            corrClosure[ibin]->Divide(gentrack[ibin]);
            
            for(Int_t ir = 0 ; ir < =corrClosure[ibin]->GetNbinsX(); ir++){
                cout <<"ir =" << ir << "ratio =" <<TMath::Abs(corrClosure[ibin]->GetBinContent(ir)-1)*100<<endl ; 
            }
        }
        
    }    
//    if(SaveFile){
//        ofstream myfile;
//        myfile.open("sysTrkNonClosure.txt");
//        double nonclosure[nbin][nrbin]={{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
//
//        for(int ibin = 0 ; ibin <nbin; ibin++){
//            
//            for(Int_t i = 1 ; i < =corrClosure[ibin]->GetNbinsX(); i++){
//                nonclosure[ibin][i-1] = 100*TMath::Abs(corrClosure[ibin]->GetBinContent(i)-1);  
//                cout <<"ibin =" <<ibin<<"sys =" <<nonclosure[ibin][i-1]<<endl;
//            }
//            myfile << nonclosure[ibin][0] <<", "<< nonclosure[ibin][1]<<", "<< nonclosure[ibin][2]<<", "<< nonclosure[ibin][3]<<", "<< nonclosure[ibin][4]<<", "<< nonclosure[ibin][5]<<"\n";
//        }
//        myfile.close();
//    }
    
    if(dis==kPt){
        TH1D * hFrame = new TH1D("hFrame","frame",500,0.,500.);
        hFrame->SetAxisRange(1.,100.,"X");
        hFrame->SetXTitle("p_{T}^{trk} (GeV/c)");
    }
    else {
        TH1D * hFrame = new TH1D("hFrame","frame",6,0.,conesize);
        hFrame->SetAxisRange(0.,conesize,"X");
        hFrame->SetXTitle("radius r");
    }
    hFrame->SetTitle("");
    fixedFontHist(hFrame,1.5, 2.0);
    
    switch (opt){
        case kDist:
            if(dis==kPt){
                hFrame->SetYTitle("dN/dp_{T}");
                if(org==kInclusive)
                    hFrame-> SetAxisRange(8.e-9,1.e,"Y");
                else 
                    hFrame-> SetAxisRange(8.e-9,8.e-4,"Y");
//                else   
//                //    hFrame-> SetAxisRange(1.,1.e6,"Y");
//                hFrame-> SetAxisRange(1.e-4,1.e2,"Y");

            }
            else {
                if(IsMC)hFrame-> SetAxisRange(1.e-7,3.e-2,"Y");
                else hFrame-> SetAxisRange(1.e3,1.e8,"Y");
                if(DoPtWeight)hFrame->SetYTitle("weighted dN/dr");  
                 else hFrame->SetYTitle("dN/dr");  
            }

            for(int ipad = 1 ; ipad <=nbin; ipad++){
                c1->cd(ipad);
                c1->cd(ipad)->SetLogy();
                if(dis==kPt)c1->cd(ipad)->SetLogx();
                hFrame->DrawCopy();
//                if(dis==kPt){
//                    if(IsMC){
//                        gentrack[nbin-ipad]->Rebin(10);
//                        gentrack[nbin-ipad]->Scale(1./10.);
//                        matchtrack[nbin-ipad]->Rebin(10);
//                        matchtrack[nbin-ipad]->Scale(1./10.);   
//                    }
//                    track[nbin-ipad]->Rebin(10);
//                    track[nbin-ipad]->Scale(1./10.);
//                    corrtrack[nbin-ipad]->Rebin(10);
//                    corrtrack[nbin-ipad]->Scale(1./10.);
//
//                }
                
                if(IsMC){
                    matchtrack[nbin-ipad]->SetMarkerStyle(1);
                    matchtrack[nbin-ipad]->SetMarkerColor(2);
                    matchtrack[nbin-ipad]->SetMarkerSize(1.5);
                    matchtrack[nbin-ipad]->SetLineColor(2);
                    matchtrack[nbin-ipad]->SetFillStyle(3004);
                    matchtrack[nbin-ipad]->SetLineColor(2);
                    matchtrack[nbin-ipad]->SetFillColor(2);
                    matchtrack[nbin-ipad]->DrawCopy("same HIST");   
                    gentrack[nbin-ipad]->SetMarkerStyle(25);
                    gentrack[nbin-ipad]->SetMarkerColor(1);
                    gentrack[nbin-ipad]->SetLineColor(1);
                    gentrack[nbin-ipad]->SetMarkerSize(1.5);
                    gentrack[nbin-ipad]->SetFillColor(1);
                    gentrack[nbin-ipad]->DrawCopy("same PE");
                }
                    track[nbin-ipad]->SetMarkerStyle(20);
                    track[nbin-ipad]->SetMarkerSize(1.5);
                    track[nbin-ipad]->SetMarkerColor(1);
                    track[nbin-ipad]->SetLineColor(1);
                    track[nbin-ipad]->DrawCopy("same PE");
 
                        corrtrack[nbin-ipad]->SetMarkerStyle(24);
                        corrtrack[nbin-ipad]->SetMarkerSize(1.5);
                        corrtrack[nbin-ipad]->SetLineColor(1);
                corrtrack[nbin-ipad]->SetMarkerColor(1);
                        corrtrack[nbin-ipad]->DrawCopy("same PE"); 
                   // if(nbin>1)drawText(Form("%d-%d%%",centr[ipad-1],centr[ipad]),0.6,0.65, 15);
             //   if(ipad==1)drawCMS(0.5,0.85,pbpbLumi);
                if(ipad==1){
                    if(IsMC)drawCMSmc(0.5,0.9,coll.Data());
                    else {
                        if(coll==kPP) drawCMSpp(0.65,0.9,ppLumi);
                        else drawCMS(0.65,0.9,pbpbLumi);
                    }
                }
                if(ipad==2){
                    if(IsMC){
                    t1->AddEntry(gentrack[nbin-ipad],"gen", "PL");
                    t1->AddEntry(matchtrack[nbin-ipad],"gen_matched", "LF2");
                    }
                    t1->AddEntry(track[nbin-ipad],"reco", "PL");
                    t1->AddEntry(corrtrack[nbin-ipad],"reco_corrected", "PL");
                    t1->Draw("same");     
                }
                if(nbin>1){
                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f (GeV/c)", leadingjetcut),0.25,0.9,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.25,0.8,17);
                    if(ipad==3)drawText("Ak PF, R=0.3",0.35,0.9,17);
             //       if(ipad==3)drawText(Form("%s %.1f<p_{T}^{trk}<%.1f GeV/c",lev.Data(), trkpt[trackbin], trkpt[trackbin+1]),0.35,0.8,17);
                    if(ipad==3)drawText(Form("%s p_{T}^{trk}>%.f GeV/c",lev.Data(), trkpt[trackbin]),0.35,0.8,17);

                }
                else {
                    if(IsMC){
                    t1->AddEntry(gentrack[nbin-ipad],"gen", "PL");
                    t1->AddEntry(matchtrack[nbin-ipad],"gen_matched", "LF2");
                    }
                    t1->AddEntry(track[nbin-ipad],"reco", "PL");
                    t1->AddEntry(corrtrack[nbin-ipad],"reco_corrected", "PL");
                    t1->Draw("same");     

                    drawText("Ak PF, R=0.3",0.65,0.78,17);
                    drawText(Form("p_{T}^{jet}>%.f (GeV/c)", leadingjetcut),0.65,0.72,17);
                    drawText(Form("%s p_{T}^{trk}>%.f GeV/c",lev.Data(), trkpt[trackbin]),0.65,0.66,17);
             //       drawText(Form("%s %.1f<p_{T}^{trk}<%.1f GeV/c",lev.Data(), trkpt[trackbin], trkpt[trackbin+1]),0.65,0.68,17);
                }

            }
            if(IsMC){    
                for(int ibin = 0 ; ibin <nbin; ibin++){
                 //   if(data==kPP)
                        c1->cd(ibin+nbin+1);
                 //   else c1->cd(ibin+5);
                    if(dis==kPt)c1->cd(ibin+nbin+1)->SetLogx();
                    
                    hFrame->SetYTitle("ratio");
                    //        hFrame->SetYTitle("corrected/gen");
                           hFrame-> SetAxisRange(0.4,1.22,"Y");
                  //  hFrame-> SetAxisRange(0.42,1.62,"Y");
                    hFrame->DrawCopy();
//                    if(dis==kPt){
//                        recoEff[nbin-ibin-1]->Rebin(10);
//                        recoEff[nbin-ibin-1]->Scale(1./10.);
//                        corrClosure[nbin-ibin-1]->Rebin(10);
//                        corrClosure[nbin-ibin-1]->Scale(1./10.);
//                    }
                    recoEff[nbin-ibin-1]->SetMarkerStyle(20);
                    recoEff[nbin-ibin-1]->SetMarkerColor(1);
                    recoEff[nbin-ibin-1]->SetMarkerSize(1.5);
                    recoEff[nbin-ibin-1]->SetLineColor(1);
                    recoEff[nbin-ibin-1]->DrawCopy("same PE");
                    
                    corrClosure[nbin-ibin-1]->SetMarkerStyle(24);
                    corrClosure[nbin-ibin-1]->SetMarkerColor(1);
                    corrClosure[nbin-ibin-1]->SetMarkerSize(1.5);
                    corrClosure[nbin-ibin-1]->SetLineColor(1);
                    corrClosure[nbin-ibin-1]->DrawCopy("same PE");
                  //  corrClosure[nbin-ibin-1]->Print();
                    if(dis==kPt)regSun(0.,1.,200.,1.,1, 1);
                    else regSun(0.,1.,0.3,1.,1, 1);
                    if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.25,0.25, 17);
                    
                }            
                t2->AddEntry(corrClosure[nbin-1],"reco_corrected/gen", "PL"); 
 //               t2->AddEntry(recoEff[nbin-1],"reco/gen_matched", "PL");
                t2->AddEntry(recoEff[nbin-1],"gen_matched/gen", "PL");
                t2->Draw("same");
            }
            if(SavePlot)
            //    c1->Print(Form("%s/%s%s%sPt%.f_%.fEtaCut%.f%sTrackPtCut%.f%s_%d.eps",plotsdir.Data(),coll.Data(),met.Data(),trig.Data(),pt[0],pt[1], etacut,select.Data(),tracklimit, dep.Data(),date));               
               c1->Print(Form("%s/%s%s%sPtEtaCut%.f%s%s2013ppHiIterativeTrackPtCut%.f_%.f%s.pdf",plotsdir.Data(),coll.Data(),met.Data(),trig.Data(),etacut,lev.Data(),effTab.Data(), trkpt[trackbin]*10,trkpt[trackbin+1]*10, dep.Data()));   
            if(SaveFile){
                TFile * outf = new TFile(Form("%s%s%s%sTrackRadius.root", type.Data(), coll.Data(), effTab.Data(), lev.Data()), "UPDATE");
                for(int ibin = 0 ; ibin <nbin; ibin++){
                    track[ibin]->Write();
                    corrtrack[ibin]->Write();
                    corrClosure[ibin]->Write();
                }
                outf->Close();
            }

            c1->Update();   
            break ;
        case kRatio:
            if(IsMC){
                if(dis==kPt){
                    hFrame->SetXTitle("p_{T}^{gen} (GeV/c)");
                    hFrame->SetYTitle("dN/dp_{T}");
                    // hFrame->SetYTitle("1/N_{jets}dN/dp_{T}");
                    hFrame-> SetAxisRange(1.e-9,1.e2,"Y");
                }
                else {
                    hFrame-> SetAxisRange(1.e3,1.e8,"Y");
                    if(DoPtWeight)hFrame->SetYTitle("weighted dN/dr");
                    else hFrame->SetYTitle("1/N_{jets}dN/dr");
                }
                
                for(int ipad = 1 ; ipad <=nbin; ipad++){
                    c1->cd(ipad);
                    c1->cd(ipad)->SetLogy();
                    if(dis==kPt)c1->cd(ipad)->SetLogx();
                    hFrame->DrawCopy();
                    ppgentrack[nbin-ipad]->SetMarkerStyle(24);
                    ppgentrack[nbin-ipad]->SetMarkerColor(1);
                    ppgentrack[nbin-ipad]->SetLineColor(1);
                    ppgentrack[nbin-ipad]->SetMarkerSize(1.5);
                    ppgentrack[nbin-ipad]->SetFillColor(1);
                    ppgentrack[nbin-ipad]->DrawCopy("same PE");
                    gentrack[nbin-ipad]->SetMarkerStyle(29);
                    gentrack[nbin-ipad]->SetMarkerColor(2);
                    gentrack[nbin-ipad]->SetLineColor(2);
                    gentrack[nbin-ipad]->SetMarkerSize(1.5);
                    gentrack[nbin-ipad]->SetFillColor(2);
                    gentrack[nbin-ipad]->DrawCopy("same PE");
                    // if(nbin>1)drawText(Form("%d-%d%%",centr[ipad-1],centr[ipad]),0.6,0.65, 15);
                    //   if(ipad==1)drawCMS(0.5,0.85,pbpbLumi);
                    if(ipad==1){
                        drawCMSmc(0.65,0.9,coll.Data());
                    }
                    if(ipad==2){
                        t1->AddEntry(gentrack[nbin-ipad]," Corr gen", "PL");
                        t1->AddEntry(ppgentrack[nbin-ipad],"Raw gen", "PL");
                        t1->Draw("same");
                    }
                    if(nbin>1){
                        if(ipad==4)drawText(Form("p_{T}^{jet}>%.f (GeV/c)", leadingjetcut),0.25,0.9,17);
                        if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.25,0.8,17);
                        if(ipad==3)drawText("Ak PF, R=0.3",0.35,0.9,17);
                              if(ipad==3)drawText(Form("%s %.1f<p_{T}^{trk}<%.1f GeV/c",lev.Data(), trkpt[trackbin], trkpt[trackbin+1]),0.35,0.8,17);
                     //   if(ipad==3)drawText(Form("%s p_{T}^{trk}>%.f GeV/c",lev.Data(), trkpt[trackbin]),0.35,0.8,17);
                        
                    }
                    else {
                        t1->AddEntry(gentrack[nbin-ipad],"Corr gen", "PL");
                        t1->AddEntry(ppgentrack[nbin-ipad],"Raw gen", "PL");
                        t1->Draw("same");
                        
                        drawText("Ak PF, R=0.3",0.65,0.78,17);
                        drawText(Form("p_{T}^{jet}>%.f (GeV/c)", leadingjetcut),0.65,0.72,17);
                        drawText(Form("%s p_{T}^{trk}>%.f GeV/c",lev.Data(), trkpt[trackbin]),0.65,0.68,17);
                        //    drawText(Form("%s %.1f<p_{T}^{trk}<%.1f GeV/c",lev.Data(), trkpt[trackbin], trkpt[trackbin+1]),0.65,0.68,17);
                    }
                    
                }
                //         if(IsMC){
                for(int ibin = 0 ; ibin <nbin; ibin++){
                    if(data==kPP){
                        if(dis==kPt) c1->cd(ibin+2)->SetLogx();
                        else c1->cd(ibin+2);
                    }
                    else {
                        if(dis==kPt)  c1->cd(ibin+5)->SetLogx();
                        else c1->cd(ibin+5);
                    }
                    
                    //    hFrame->SetYTitle("ratio");
                    hFrame->SetYTitle(Form("Corr/Raw Gen"));
                    //                    hFrame->SetYTitle(Form("Gen ppTracking/HITracking"));
                    hFrame-> SetAxisRange(0.8,1.12,"Y");
                    hFrame->DrawCopy();
                    genratio[nbin-ibin-1]->SetMarkerStyle(20);
                    genratio[nbin-ibin-1]->SetMarkerColor(1);
                    genratio[nbin-ibin-1]->SetMarkerSize(1.5);
                    genratio[nbin-ibin-1]->SetLineColor(1);
                    genratio[nbin-ibin-1]->DrawCopy("same PE");
                    if(dis==kPt)regSun(0.,1.,200.,1.,1, 1);
                    else regSun(0.,1.,0.3,1.,1, 1);
                    if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.25,0.25, 17);
                    
                }
                if(SavePlot)
                    //    c1->Print(Form("%s/%s%s%sPt%.f_%.fEtaCut%.f%sTrackPtCut%.f%s_%d.eps",plotsdir.Data(),coll.Data(),met.Data(),trig.Data(),pt[0],pt[1], etacut,select.Data(),tracklimit, dep.Data(),date));
                    c1->Print(Form("%s/%s%s%s%sPtEtaCut%.f%s%sGenDistDiffFilePtCut%.f_%.f%s.pdf",plotsdir.Data(),type.Data(), coll.Data(),met.Data(),trig.Data(),etacut,lev.Data(),effTab.Data(), trkpt[trackbin]*10,trkpt[trackbin+1]*10, dep.Data()));
                c1->Update();
            }
            break ;

    }
}
//---------------------------------------------------
void drawText(const char *text, float xp, float yp, int textSize=15){
    TLatex *tex = new TLatex(xp,yp,text);
    tex->SetTextFont(63);
    //tex->SetTextSize(20);
    tex->SetTextSize(textSize);
    //tex->SetTextSize(0.05);                                                                   
    tex->SetTextColor(kBlack);
    tex->SetLineWidth(1);
    tex->SetNDC();
    tex->Draw();
}

void drawText2(const char *text, float xp, float yp, int textSize=18){
    TLatex *tex = new TLatex(xp,yp,text);
    tex->SetTextFont(63);
    tex->SetTextSize(textSize);
    tex->SetTextColor(kBlack);
    tex->SetLineWidth(1);
    tex->Draw();
}

void makeMultiPanelCanvas(TCanvas*& canv,
                          const Int_t columns,
                          const Int_t rows,
                          const Float_t leftOffset,
                          const Float_t bottomOffset,
                          const Float_t leftMargin,
                          const Float_t bottomMargin,
                          const Float_t edge) {
    if (canv==0) {
        Error("makeMultiPanelCanvas","Got null canvas.");
        return;
    }
    canv->Clear();
    
    TPad* pad[columns][rows];
    
    Float_t Xlow[columns];
    Float_t Xup[columns];
    Float_t Ylow[rows];
    Float_t Yup[rows];
    Float_t PadWidth = 
    (1.0-leftOffset)/((1.0/(1.0-leftMargin)) +
                      (1.0/(1.0-edge))+(Float_t)columns-2.0);
    Float_t PadHeight =
    (1.0-bottomOffset)/((1.0/(1.0-bottomMargin)) +
                        (1.0/(1.0-edge))+(Float_t)rows-2.0);
    Xlow[0] = leftOffset;
    Xup[0] = leftOffset + PadWidth/(1.0-leftMargin);
    Xup[columns-1] = 1;
    Xlow[columns-1] = 1.0-PadWidth/(1.0-edge);
    
    Yup[0] = 1;
    Ylow[0] = 1.0-PadHeight/(1.0-edge);
    Ylow[rows-1] = bottomOffset;
    Yup[rows-1] = bottomOffset + PadHeight/(1.0-bottomMargin);
    
    for(Int_t i=1;i<columns-1;i++) {
        Xlow[i] = Xup[0] + (i-1)*PadWidth;
        Xup[i] = Xup[0] + (i)*PadWidth;
    }
    Int_t ct = 0;
    for(Int_t i=rows-2;i>0;i--) {
        Ylow[i] = Yup[rows-1] + ct*PadHeight;
        Yup[i] = Yup[rows-1] + (ct+1)*PadHeight;
        ct++;
    }
    
    TString padName;
    for(Int_t i=0;i<columns;i++) {
        for(Int_t j=0;j<rows;j++) {
            canv->cd();
            padName = Form("p_%d_%d",i,j);
            pad[i][j] = new TPad(padName.Data(),padName.Data(),
                                 Xlow[i],Ylow[j],Xup[i],Yup[j]);
            if(i==0) pad[i][j]->SetLeftMargin(leftMargin);
            else pad[i][j]->SetLeftMargin(0);
            
            if(i==(columns-1)) pad[i][j]->SetRightMargin(edge);
            else pad[i][j]->SetRightMargin(0);
            
            if(j==0) pad[i][j]->SetTopMargin(edge);
            else pad[i][j]->SetTopMargin(0);
            
            if(j==(rows-1)) pad[i][j]->SetBottomMargin(bottomMargin);
            else pad[i][j]->SetBottomMargin(0);
            
            pad[i][j]->Draw();
            pad[i][j]->cd();
            pad[i][j]->SetNumber(columns*j+i+1);
        }
    }
}
void drawCMSmc(float px, float py, TString coll) {
    TLatex *cms = new TLatex(px,py,"CMS Preliminary");
    cms->SetTextFont(63);
    cms->SetTextSize(17);
    cms->SetNDC();
    cms->Draw();
    if(coll=="HI")TLatex *lumi = new TLatex(px,py-0.08,"PYTHIA+HYDJET");
    else TLatex *lumi = new TLatex(px,py-0.05,"PYTHIA 2013");
//    else TLatex *lumi = new TLatex(px,py-0.05,"PYTHIA pp-tracking");
    lumi->SetTextFont(63);
    lumi->SetTextSize(15);
    lumi->SetNDC();
    lumi->Draw();
}

void drawCMS(float px, float py, float nLumi) {
    TLatex *cms = new TLatex(px,py,"CMS Preliminary");
    cms->SetTextFont(63);
    cms->SetTextSize(17);
    cms->SetNDC();
    cms->Draw();
    TLatex *lumi = new TLatex(px,py-0.10,Form("#intL dt = %.1f #mub^{-1}",nLumi));
    lumi->SetTextFont(63);
    lumi->SetTextSize(15);
    lumi->SetNDC();
    lumi->Draw();
}

void drawCMSpp(float px, float py, float nLumi) {
    TLatex *cms = new TLatex(px,py,"CMS Preliminary");
    cms->SetTextFont(63);
    cms->SetTextSize(15);
    cms->SetNDC();
    cms->Draw();
    TLatex *lumi = new TLatex(px,py-0.05,Form("#intL dt = %.1f nb^{-1}",nLumi));
    lumi->SetTextFont(63);
    lumi->SetTextSize(15);
    lumi->SetNDC();
    lumi->Draw();
}
void fixedFontHist(TH1 * h, Float_t xoffset=1.3, Float_t yoffset=1.2)
{
    h->SetLabelFont(43,"X");
    h->SetLabelFont(43,"Y");
    //h->SetLabelOffset(0.01);
    h->SetLabelSize(22);
    h->SetTitleFont(43);
    h->SetTitleSize(16);
    h->SetLabelSize(18,"Y");
    h->SetLabelSize(18,"X");
    h->SetTitleFont(43,"Y");
    h->SetTitleSize(20,"Y");
    h->SetTitleSize(20,"X");
    h->SetTitleOffset(xoffset,"X");
    h->SetTitleOffset(yoffset,"Y");
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->CenterTitle();
}

void drawSys(TGraph *h, double *sys, double width=5, int theColor= kYellow, int fillStyle = -1, int lineStyle = -1)
{
    for (int i=0;i<h->GetN();i++)
    {
        double val;
        double theX;
        h->GetPoint(i,theX,val);
        double err = val * sys[i]/100;
        TBox *b = new TBox(theX-width,val-err,theX+width,val+err);
        
        b->SetLineColor(theColor);
        b->SetFillColor(theColor);
        if ( fillStyle > -1 ) b->SetFillStyle(fillStyle);
        if ( lineStyle > -1 ) b->SetLineStyle(lineStyle);
        
        b->Draw();
    }
}
void onSun(double x1=0,double y1=0,double x2=1,double y2=1,int color=1, double width=1)
{
    TLine* t1 = new TLine(x1,y1,x2,y2);
    t1->SetLineWidth(width);
    t1->SetLineStyle(1);
    t1->SetLineColor(color);
    t1->Draw();
}
void regSun(double x1=0,double y1=0,double x2=1,double y2=1,int color=1, double width=1)
{
    TLine* t1 = new TLine(x1,y1,x2,y2);
    t1->SetLineWidth(width);
    t1->SetLineStyle(3);
    t1->SetLineColor(color);
    t1->Draw();
}
