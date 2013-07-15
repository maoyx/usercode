/*
 *  Rocket plot from data and MC in pp and HI
 *
 *
 *  Created by ymao on 07/12/13.
 *  Copyright 2013 Vanderbilt University (US). All rights reserved.
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
//#include "./rootlogon.h"    //// A good drawing style is defined here
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<iostream>
#include<fstream>
#include<string>
#include <iostream>
using namespace std;

//define the kinematics cuts for the analysis
const double conesize = 0.3 ;
const double deltacone = 0.05 ;
const double etalimit = 0.3 ;
const double etacut = 2.0 ;
const double dphicut = 7*TMath::Pi()/8. ;
const double leadingjetcut = 100. ;
const double subleadjetcut = 40. ;
const double trackcut = 1.0;
bool SavePlot = kTRUE ;

bool DoSmear=kFALSE ;

const double pt[]={100., 500.};
//const double pt[]={100., 120., 140., 160., 200., 300., 500.};

const int nptbin = 4 ;
//const double trkpt[]={1., 2., 4., 8., 20., 50., 500.};

//const double pt[]={100., 120., 140.,500.};
const double subpt[]={40., 50., 60., 80., 100., 120., 200.};

const Double_t pbpbLumi = 150 ;
//const Double_t ppLumi = 230 ;
const Double_t ppLumi = 5.3 ;

enum Level_t {kRAW, kBkg, kSub, kSigBkg} ;

//---------------------------------------------------------------------

void ppHIRocketPlot(Level_t lev =kSub){
    if(conesize==0.3){
        const int nrbin = 6 ;
        double rbin[nrbin+1]={0.0,0.05,0.10, 0.15, 0.20,0.25,0.30};
  
    }
    else {
        const int nrbin = 10 ;
        double rbin[nrbin+1]={0.0,0.05,0.10, 0.15, 0.20,0.25,0.30, 0.35, 0.4, 0.45, 0.5};

    }
    double rad[nrbin];
    double ratio[nrbin];
    double errR[nrbin]; 
    for(int i = 0 ; i <nrbin; i++){
        //        if(opt==kInt) rad[i]=rbin[i+1] ;
        //        else         
        rad[i]=deltacone/2.+i*deltacone ;
        ratio[i]=rad[i]/conesize;
      //  errR[i]=deltacone/2. ;
        errR[i]=0 ;
    }
    TString Norm="NormJet" ;    
    TString plotsdir = "/Users/ymao/group/CMS/plots/finalJS/CorrectDR";
    TDatime d ;
    int date = d.GetDate();  
        // for HI centrality bin or smeared pp
        const int nbin = 4 ;
        const int centr[] ={0,10,30,50,100}; 

    if(nptbin==4) {
        const double trkpt[]={1., 2., 4., 10., 300.};
    }
    else if(nptbin==1) { 
          const double trkpt[]={1., 300.};
    }
    else {
          const double trkpt[]={1., 2., 4., 8., 16., 32., 500.};
    }

//    const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/FullStatPbPb/RebinJetPt/TrkCorr4CentBin" ;
    const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/CorrectdR/NoMatchGen/NoTrkAlgoCut/MyselfTrkEff/AfterCWR" ;
    TString fileName ;
    TFile * f ;

    TString mcfileName ;
    TFile * mcf ;

    TString ppfileName ;
    TFile * ppf ;
    
    TString mcppfileName ;
    TFile * mcppf ;

    TH1F * jetpt[nbin];
    TH1F * mcjetpt[nbin];

    TH1F * ppjetpt;
    TH1F * mcppjetpt;
    
    TH2F * jetdrSumTrkPt[nbin][nptbin];

    TH2F * bkgjetdrSumTrkPt[nbin][nptbin];
    TH2F * mcjetdrSumTrkPt[nbin][nptbin];
    
    TH2F * mcbkgjetdrSumTrkPt[nbin][nptbin];
 
    TH2F * ppjetdrSumTrkPt[nptbin];
    TH2F * ppbkgjetdrSumTrkPt[nptbin];
    TH2F * mcppjetdrSumTrkPt[nptbin];
    TH2F * mcppbkgjetdrSumTrkPt[nptbin];
    
//    TProfile * SumPtdR[nbin][nptbin];
//    TProfile * bkgSumPtdR[nbin][nptbin];
    TH1F * SumPtdR[nbin][nptbin];
    TH1F * bkgSumPtdR[nbin][nptbin];
    
    TH1F * SumPtdRsub[nbin][nptbin];

    TH1F * mcSumPtdR[nbin][nptbin];
    TH1F * mcbkgSumPtdR[nbin][nptbin];
    
    TH1F * mcSumPtdRsub[nbin][nptbin];
// TProfile * SumPtdRsub[nbin][nptbin];
    TH1F * ppSumPtdR[nptbin];
    TH1F * ppbkgSumPtdR[nptbin];
    TH1F * ppSumPtdRsub[nptbin];
    TH1F * mcppSumPtdR[nptbin];
    TH1F * mcppbkgSumPtdR[nptbin];    
    TH1F * mcppSumPtdRsub[nptbin];

    TH1F * ratioSum[nbin][nptbin];
    TH1F * ratiobkgSum[nbin][nptbin];
    TH1F * ratioSubSum[nbin][nptbin];

    TH1F * ratioSigBkgSum[nbin][nptbin];
    TH1F * mcratioSigBkgSum[nbin][nptbin];

    TH1F * SumPtBindR[nbin][nptbin];
    TH1F * bkgSumPtBindR[nbin][nptbin];
    TH1F * mcSumPtBindR[nbin][nptbin];
    TH1F * mcbkgSumPtBindR[nbin][nptbin];
    
    TH1F * ratioSigBkg[nbin][nptbin];
    TH1F * mcratioSigBkg[nbin][nptbin];

    TProfile * StackSumPt[nbin][nptbin];

    TH1::SetDefaultSumw2();

//    TCanvas  * c1 = new TCanvas("c1", "c1", kw, kh);
//    c1->SetFillColor(0);
//    c1->SetBorderSize(0);
//    c1->SetFrameBorderMode(0); 
//    gStyle->SetOptStat(0);    

    TString canv_name = "c1";
//    if(nbin>1){
//        const Double_t kw = 1000;
//        const Double_t kh = 350;
//        c1 = new TCanvas(canv_name," ",10,10,kw,kh);
//        makeMultiPanelCanvas(c1,4,1,0.0,0.0,0.2,0.2,0.02); 
//    }
//    else {
//        c1 = new TCanvas(canv_name," ",560, 560);
//    }
//    if(nbin>1){
////        if(lev==kSub){
////            const Double_t kw = 1000;
////            const Double_t kh = 350;
////            c1 = new TCanvas(canv_name," ",10,10,kw,kh);
////            makeMultiPanelCanvas(c1,nbin,1,0.0,0.0,0.2,0.2,0.02); 
////        }
////        else {
//        const Double_t kw = 1000;
//        const Double_t kh = 560;
//        c1 = new TCanvas(canv_name," ",10,10,kw,kh);
//        makeMultiPanelCanvas(c1,4,2,0.0,0.0,0.2,0.2,0.02);
////        }
//    }
//    else {
    c1 = new TCanvas(canv_name," ",500, 500);
//        c1 = new TCanvas(canv_name," ",1000, 510);
      //  c1->Divide(1, 2, 0, 0);
//        makeMultiPanelCanvas(c1,2,1, 0.0,0.0,0.15,0.1,0.02);
//    }

    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);

    
//    if(nbin>1)TLegend *t1=new TLegend(0.2,0.6,0.85,0.95);
//    else
    TLegend *t1=new TLegend(0.65,0.68,0.85,0.92);
    t1->SetFillColor(0);
    t1->SetBorderSize(0);
    t1->SetFillStyle(0);
    t1->SetTextFont(63);
    t1->SetTextSize(17);
    TLegend *t2=new TLegend(0.05,0.85,0.35,0.92);
    t2->SetFillColor(0);
    t2->SetBorderSize(0);
    t2->SetFillStyle(0);
    t2->SetTextFont(16);
    t2->SetTextSize(16);
    TLegend *t3=new TLegend(0.25,0.85,0.5,0.92);
    t3->SetFillColor(0);
    t3->SetBorderSize(0);
    t3->SetFillStyle(0);
    t3->SetTextFont(16);
    t3->SetTextSize(16);
    // open the data file
    mcppfileName = Form("mergedCSdiff_MCPP_Ak3PFIncJetPt50_2013HistIterTrkCorrtest%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin%d_pp2013_P01_prod22_v81_merged.root", trackcut, nptbin);
    ppfileName =  Form("DATAPP_Ak3PFIncJetPt100_2013HistIterTrkCorrtest%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin%d_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", trackcut, nptbin);

    mcfileName=Form("mergedCSdiff_MCHI_AkPu3PFIncJetPt100_2013HistIterTrkCorrtest%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%d_Dijet_HydjetDrum_v27v28.root", trackcut, nbin,nptbin);
    fileName=Form("DATAHI_AkPu3PFIncJetPt100_2013HistIterTrkCorrtest%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin%d_promptskim-hihighpt-hltjet80-pt90-v20.root", trackcut, nbin,nptbin);


    ppf = TFile::Open(Form("%s/%s", kHomeDir, ppfileName.Data()), "readonly");
    mcppf = TFile::Open(Form("%s/%s", kHomeDir, mcppfileName.Data()), "readonly");
    f = TFile::Open(Form("%s/%s", kHomeDir, fileName.Data()), "readonly");
    mcf = TFile::Open(Form("%s/%s", kHomeDir, mcfileName.Data()), "readonly");

    ppjetpt=(TH1F*)ppf->Get(Form("jetpt_%d-%d%%",centr[0],centr[nbin]));
    ppjetpt->Sumw2();
    mcppjetpt=(TH1F*)mcppf->Get(Form("jetpt_%d-%d%%",centr[0],centr[nbin]));
    mcppjetpt->Sumw2();
    for(int ipt =0 ; ipt <nptbin; ipt++){
        ppSumPtdR[ipt]=(TH1F*)ppf->Get(Form("JetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[0],centr[nbin]));
        ppSumPtdR[ipt]->Sumw2();
        
        ppbkgSumPtdR[ipt]=(TH1F*)ppf->Get(Form("bkgJetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[0],centr[nbin]));
        ppbkgSumPtdR[ipt]->Sumw2();
        
        mcppSumPtdR[ipt]=(TH1F*)mcppf->Get(Form("JetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[0],centr[nbin]));
        mcppSumPtdR[ipt]->Sumw2();
        mcppbkgSumPtdR[ipt]=(TH1F*)mcppf->Get(Form("bkgJetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[0],centr[nbin]));
        mcppbkgSumPtdR[ipt]->Sumw2();
        
        ppSumPtdR[ipt]->Scale(1./ppjetpt->Integral());
        ppbkgSumPtdR[ipt]->Scale(1./ppjetpt->Integral());
        mcppSumPtdR[ipt]->Scale(1./mcppjetpt->Integral());
        mcppbkgSumPtdR[ipt]->Scale(1./mcppjetpt->Integral());
        
        ppSumPtdRsub[ipt]=(TH1F*)ppSumPtdR[ipt]->Clone(Form("ppJetPt%.f_%.fSumTrk%.f_%.f_SubdRCen%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[0],centr[nbin]));
        ppSumPtdRsub[ipt]->Add(ppbkgSumPtdR[ipt], -1);
        mcppSumPtdRsub[ipt]=(TH1F*)mcppSumPtdR[ipt]->Clone(Form("mcppJetPt%.f_%.fSumTrk%.f_%.f_SubdRCen%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[0],centr[nbin]));
        mcppSumPtdRsub[ipt]->Add(mcppbkgSumPtdR[ipt], -1);

    }
    
    for(int ibin = 0 ; ibin <nbin; ibin++){
        jetpt[ibin]=(TH1F*)f->Get(Form("jetpt_%d-%d%%",centr[ibin],centr[ibin+1]));
        jetpt[ibin]->Sumw2();
        mcjetpt[ibin]=(TH1F*)mcf->Get(Form("jetpt_%d-%d%%",centr[ibin],centr[ibin+1])); 
        mcjetpt[ibin]->Sumw2();
        for(int ipt =0 ; ipt <nptbin; ipt++){
            SumPtdR[ibin][ipt]=(TH1F*)f->Get(Form("JetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[ibin],centr[ibin+1])); 
            SumPtdR[ibin][ipt]->Sumw2();
          //              SumPtdR[ibin][ipt]->Scale(1./SumPtdR[ibin][ipt]->Integral());
          //  SumPtdR[ibin][ipt]->Scale(6./SumPtdR[ibin][ipt]->GetEntries());
            
            bkgSumPtdR[ibin][ipt]=(TH1F*)f->Get(Form("bkgJetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[ibin],centr[ibin+1])); 
            bkgSumPtdR[ibin][ipt]->Sumw2();
           //             bkgSumPtdR[ibin][ipt]->Scale(1./bkgSumPtdR[ibin][ipt]->Integral());
        //    bkgSumPtdR[ibin][ipt]->Scale(6./bkgSumPtdR[ibin][ipt]->GetEntries());
            
            mcSumPtdR[ibin][ipt]=(TH1F*)mcf->Get(Form("JetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[ibin],centr[ibin+1])); 
            mcSumPtdR[ibin][ipt]->Sumw2();
             //   mcSumPtdR[ibin][ipt]->Scale(1./mcSumPtdR[ibin][ipt]->Integral());
         //   mcSumPtdR[ibin][ipt]->Scale((6.)/mcSumPtdR[ibin][ipt]->GetEntries());
            mcbkgSumPtdR[ibin][ipt]=(TH1F*)mcf->Get(Form("bkgJetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[ibin],centr[ibin+1])); 
            mcbkgSumPtdR[ibin][ipt]->Sumw2();
           //             mcbkgSumPtdR[ibin][ipt]->Scale(1./mcbkgSumPtdR[ibin][ipt]->Integral());
         //   mcbkgSumPtdR[ibin][ipt]->Scale((6.)/(bkgSumPtdR[ibin][ipt]->GetEntries()));
          

            SumPtdR[ibin][ipt]->Scale(1./jetpt[ibin]->Integral());
            bkgSumPtdR[ibin][ipt]->Scale(1./jetpt[ibin]->Integral());
            mcSumPtdR[ibin][ipt]->Scale(1./mcjetpt[ibin]->Integral());
            mcbkgSumPtdR[ibin][ipt]->Scale(1./mcjetpt[ibin]->Integral());

            SumPtdRsub[ibin][ipt]=(TH1F*)SumPtdR[ibin][ipt]->Clone(Form("JetPt%.f_%.fSumTrk%.f_%.f_SubdRCen%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[ibin],centr[ibin+1]));
            SumPtdRsub[ibin][ipt]->Add(bkgSumPtdR[ibin][ipt], -1);
            mcSumPtdRsub[ibin][ipt]=(TH1F*)mcSumPtdR[ibin][ipt]->Clone(Form("mcJetPt%.f_%.fSumTrk%.f_%.f_SubdRCen%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[ibin],centr[ibin+1]));
            mcSumPtdRsub[ibin][ipt]->Add(mcbkgSumPtdR[ibin][ipt], -1);
//            SumPtdRsub[ibin][ipt]->Scale(1./jetpt[ibin]->Integral());
//            mcSumPtdRsub[ibin][ipt]-> Scale(1./mcjetpt[ibin]->Integral());           
            ratioSubSum[ibin][ipt]=(TH1F*)SumPtdRsub[ibin][ipt]->Clone(Form("JetPt%.f_%.fRatioSumTrk%.f_%.f_SubdRCen%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[ibin],centr[ibin+1]));
            ratioSubSum[ibin][ipt]->Divide(mcSumPtdRsub[ibin][ipt]);
                
                ratioSum[ibin][ipt]=(TH1F*)SumPtdR[ibin][ipt]->Clone(Form("JetPt%.f_%.fRatioSumTrk%.f_%.f_dRCen%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[ibin],centr[ibin+1]));
            ratioSum[ibin][ipt]->Divide(mcSumPtdR[ibin][ipt]);
            
            ratiobkgSum[ibin][ipt]=(TH1F*)bkgSumPtdR[ibin][ipt]->Clone(Form("JetPt%.f_%.fRatiobkgSumTrk%.f_%.f_dRCen%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[ibin],centr[ibin+1]));
            ratiobkgSum[ibin][ipt]->Divide(mcbkgSumPtdR[ibin][ipt]);
                ratioSigBkgSum[ibin][ipt]=(TH1F*)SumPtdRsub[ibin][ipt]->Clone(Form("JetPt%.f_%.fRatioBkgSigSumTrk%.f_%.f_dRCen%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[ibin],centr[ibin+1]));
                ratioSigBkgSum[ibin][ipt]->Divide(bkgSumPtdR[ibin][ipt]);
                
                mcratioSigBkgSum[ibin][ipt]=(TH1F*)mcSumPtdRsub[ibin][ipt]->Clone(Form("mcJetPt%.f_%.fRatioBkgSigSumTrk%.f_%.f_dRCen%d-%d%%",pt[0],pt[0+1],trkpt[ipt],trkpt[ipt+1], centr[ibin],centr[ibin+1]));
                mcratioSigBkgSum[ibin][ipt]->Divide(mcbkgSumPtdR[ibin][ipt]);
        } 
    }


    TH1F * dummy = new TH1F("dummy", "dummy", 100, -0.005, 1.);
    dummy->SetTitle("") ;
    dummy->SetStats(kFALSE) ;
    dummy->SetAxisRange(-0.005, 0.307, "X") ;
    //    dummy->SetAxisRange(0., 1., "X") ;
    dummy->GetXaxis()->SetTitle("r");
    dummy->GetYaxis()->SetTitle("#Sigma p_{T}^{track}/p_{T}^{jet}");
    if(nbin==1) fixedFontHist(dummy,1.1, 1.3);
    else fixedFontHist(dummy,1.2, 1.4);
    switch(lev){
        case kRAW:
            c1->cd(1);
            c1->cd(1)->SetLogy();
            dummy->SetAxisRange(8.e-4, 5., "Y") ;
            dummy->DrawCopy();
            for(int ipt = nptbin-1 ; ipt >=0; ipt--){
                mcppSumPtdR[ipt]->SetMarkerStyle(0);
                mcppSumPtdR[ipt]->SetMarkerSize(1.5);
                if(ipt==2){
                    mcppSumPtdR[ipt]->SetMarkerColor(ipt+2);
                    mcppSumPtdR[ipt]->SetFillColor(ipt+2);
                    mcppSumPtdR[ipt]->SetFillStyle(3001);
                    mcppSumPtdR[ipt]->SetLineColor(ipt+2);
                    ppSumPtdR[ipt]->SetMarkerStyle(20);
                    ppSumPtdR[ipt]->SetMarkerSize(1.5);
                    ppSumPtdR[ipt]->SetMarkerColor(ipt+2);
                    ppSumPtdR[ipt]->SetFillColor(ipt+2);
                    ppSumPtdR[ipt]->SetLineColor(ipt+1);
                }
                else {
                    mcppSumPtdR[ipt]->SetMarkerColor(ipt+1);
                    mcppSumPtdR[ipt]->SetFillColor(ipt+1);
                    mcppSumPtdR[ipt]->SetFillStyle(3001);
                    mcppSumPtdR[ipt]->SetLineColor(ipt+1);
                    ppSumPtdR[ipt]->SetMarkerStyle(20);
                    ppSumPtdR[ipt]->SetMarkerSize(1.5);
                    ppSumPtdR[ipt]->SetMarkerColor(ipt+1);
                    ppSumPtdR[ipt]->SetFillColor(ipt+1);
                    ppSumPtdR[ipt]->SetLineColor(ipt+1);
                }
                mcppSumPtdR[ipt]->DrawCopy("same eHIST");
                ppSumPtdR[ipt]->DrawCopy("same PE");
                drawCMSpp(0.2,0.9,ppLumi);
                
                    t1->AddEntry(ppSumPtdR[ipt],Form("%.f<p_{T}^{trk}<%.f (GeV/c)", trkpt[ipt],trkpt[ipt+1]), "PL");
                    if(ipt==0){
                        t2->AddEntry(mcppSumPtdR[ipt], "MC", "LF2");
                        t2->AddEntry(ppSumPtdR[ipt], "DATA", "PL");
                    }

                
                cout <<"ipt =" << ipt<<"Int ="<<ppSumPtdR[ipt]->Integral()<<endl ;
                
            }
            c1->cd(2);
            c1->cd(2)->SetLogy();
            dummy->SetAxisRange(8.e-4, 5., "Y") ;
            dummy->DrawCopy();
                for(int ipt = nptbin-1 ; ipt >=0; ipt--){
                    mcSumPtdR[0][ipt]->SetMarkerStyle(0);  
                    mcSumPtdR[0][ipt]->SetMarkerSize(1.5);
                    if(ipt==2){
                        mcSumPtdR[0][ipt]->SetMarkerColor(ipt+2);
                        mcSumPtdR[0][ipt]->SetFillColor(ipt+2);
                        mcSumPtdR[0][ipt]->SetFillStyle(3001);
                        mcSumPtdR[0][ipt]->SetLineColor(ipt+2);
                        SumPtdR[0][ipt]->SetMarkerStyle(20);
                        SumPtdR[0][ipt]->SetMarkerSize(1.5);
                        SumPtdR[0][ipt]->SetMarkerColor(ipt+2);
                        SumPtdR[0][ipt]->SetFillColor(ipt+2);
                        SumPtdR[0][ipt]->SetLineColor(ipt+1);
                    }
                    else {
                    mcSumPtdR[0][ipt]->SetMarkerColor(ipt+1);   
                    mcSumPtdR[0][ipt]->SetFillColor(ipt+1); 
                    mcSumPtdR[0][ipt]->SetFillStyle(3001);
                    mcSumPtdR[0][ipt]->SetLineColor(ipt+1);
                        SumPtdR[0][ipt]->SetMarkerStyle(20);
                        SumPtdR[0][ipt]->SetMarkerSize(1.5);
                        SumPtdR[0][ipt]->SetMarkerColor(ipt+1);
                        SumPtdR[0][ipt]->SetFillColor(ipt+1);
                        SumPtdR[0][ipt]->SetLineColor(ipt+1);
                    }
                    mcSumPtdR[0][ipt]->DrawCopy("same eHIST");
                    SumPtdR[0][ipt]->DrawCopy("same PE"); 
                    
                        t1->AddEntry(SumPtdR[0][ipt],Form("%.f<p_{T}^{track}<%.f (GeV/c)", trkpt[ipt],trkpt[ipt+1]), "PL");   
                        if(ipt==0){
                            t2->AddEntry(mcSumPtdR[0][ipt], "MC", "LF2");      
                            t2->AddEntry(SumPtdR[0][ipt], "DATA", "PL");
                        }

                    
                    
                }
                drawText(Form("%d-%d%%",centr[0],centr[1]),0.25,0.92, 17);
                        drawCMS(0.2,0.9,pbpbLumi);

                    t1->Draw("same");
                    drawText(Form("p_{T}^{jet}>%.f GeV/c", pt[0]),0.68,0.6,17);
                    drawText(Form(" %.1f < |#eta|_{jet} < %.f",etalimit, etacut),0.68,0.55,17);
                    drawText("Ak PF, R=0.3",0.68,0.5,17);
                    t2->Draw("same");
            
            if(SavePlot)
                //    c1->Print(Form("%s/%sJetPt%.f_%.fSumTrkFraction%s.eps",plotsdir.Data(),type.Data(), pt[0],pt[0+1],level.Data()));               
                c1->Print(Form("%s/ppHIJetPt%.f_%.fTrkBin%dLogYRawSumPtRadius.gif",plotsdir.Data(),pt[0],pt[0+1],nptbin));
            c1->Update();   
            break ;
        case kSub:
            c1->cd(1);
            c1->cd(1)->SetLogy();
            dummy->SetAxisRange(5.e-5, 10., "Y") ;
//            dummy->SetAxisRange(8.e-4, 20., "Y") ;
            dummy->DrawCopy();
            for(int ipt = nptbin-1 ; ipt >=0; ipt--){
                mcppSumPtdR[ipt]->SetMarkerStyle(0);
                mcppSumPtdR[ipt]->SetMarkerSize(1.5);
//                if(ipt>=2){
//                    mcppSumPtdRsub[ipt]->SetMarkerColor(ipt+2);
//                    mcppSumPtdRsub[ipt]->SetFillColor(ipt+2);
//                    mcppSumPtdRsub[ipt]->SetFillStyle(3001);
//                    mcppSumPtdRsub[ipt]->SetLineColor(ipt+2);
//                    ppSumPtdRsub[ipt]->SetMarkerStyle(20);
//                    ppSumPtdRsub[ipt]->SetMarkerSize(1.5);
//                    ppSumPtdRsub[ipt]->SetMarkerColor(ipt+2);
//                    ppSumPtdRsub[ipt]->SetFillColor(ipt+2);
//                    ppSumPtdRsub[ipt]->SetLineColor(ipt+1);
//                }
//                else {
                    mcppSumPtdRsub[ipt]->SetMarkerColor(ipt+1);
                     mcppSumPtdRsub[ipt]->SetMarkerStyle(1);
                    mcppSumPtdRsub[ipt]->SetFillColor(ipt+1);
                    mcppSumPtdRsub[ipt]->SetFillStyle(3004);
                    mcppSumPtdRsub[ipt]->SetLineColor(ipt+1);
                    ppSumPtdRsub[ipt]->SetMarkerStyle(20+ipt);
                    ppSumPtdRsub[ipt]->SetMarkerSize(1.5);
                    ppSumPtdRsub[ipt]->SetMarkerColor(ipt+1);
                  //  ppSumPtdRsub[ipt]->SetFillColor(ipt+1);
                    ppSumPtdRsub[ipt]->SetLineColor(ipt+1);
                if(ipt==1){
                    mcppSumPtdRsub[ipt]->SetMarkerColor(3);
                    mcppSumPtdRsub[ipt]->SetLineColor(3);
                    mcppSumPtdRsub[ipt]->SetFillColor(3);
                    ppSumPtdRsub[ipt]->SetMarkerColor(3);
                    ppSumPtdRsub[ipt]->SetLineColor(3);
                }
                if(ipt==2){
                    mcppSumPtdRsub[ipt]->SetMarkerColor(2);
                    mcppSumPtdRsub[ipt]->SetLineColor(2);
                    mcppSumPtdRsub[ipt]->SetFillColor(2);
                    ppSumPtdRsub[ipt]->SetMarkerColor(2);
                    ppSumPtdRsub[ipt]->SetLineColor(2);
                    ppSumPtdRsub[ipt]->SetMarkerStyle(24);
                }
                if(ipt==3){
                    ppSumPtdRsub[ipt]->SetMarkerStyle(25);
                }
                
//                }
                mcppSumPtdRsub[ipt]->DrawCopy("same eHIST");
                ppSumPtdRsub[ipt]->DrawCopy("same PE");
                drawCMSpp(0.3,0.85,ppLumi);
                
                t1->AddEntry(ppSumPtdRsub[ipt],Form("%.f<p_{T}^{track}<%.f ", trkpt[ipt],trkpt[ipt+1]), "PL");
//                t1->AddEntry(ppSumPtdRsub[ipt],Form("%.f<p_{T}^{trk}<%.f (GeV/c)", trkpt[ipt],trkpt[ipt+1]), "PL");
                t1->Draw("same");
                if(ipt==0){
                    t3->AddEntry(mcppSumPtdRsub[ipt], "PYTHIA", "LF2");
                    t3->AddEntry(ppSumPtdRsub[ipt], "Data, #sqrt{s} = 2.76 TeV", "PL");
                }
                t3->Draw("same");
 
                
                cout <<"ipt =" << ipt<<"Int ="<<ppSumPtdR[ipt]->Integral()<<endl ;
                
            }
//            c1->cd(2);
//            c1->cd(2)->SetLogy();
//            dummy->SetAxisRange(1.e-5, 10., "Y") ;
//            dummy->DrawCopy();
//            for(int ipt = nptbin-1 ; ipt >=0; ipt--){
//                mcSumPtdRsub[0][ipt]->SetMarkerStyle(0);
//                mcSumPtdRsub[0][ipt]->SetMarkerSize(1.5);
////                if(ipt==2){
////                    mcSumPtdRsub[0][ipt]->SetMarkerColor(ipt+2);
////                    mcSumPtdRsub[0][ipt]->SetFillColor(ipt+2);
////                    mcSumPtdRsub[0][ipt]->SetFillStyle(3001);
////                    mcSumPtdRsub[0][ipt]->SetLineColor(ipt+2);
////                    SumPtdRsub[0][ipt]->SetMarkerStyle(20);
////                    SumPtdRsub[0][ipt]->SetMarkerSize(1.5);
////                    SumPtdRsub[0][ipt]->SetMarkerColor(ipt+2);
////                    SumPtdRsub[0][ipt]->SetFillColor(ipt+2);
////                    SumPtdRsub[0][ipt]->SetLineColor(ipt+1);
////                }
////                else {
//                    mcSumPtdRsub[0][ipt]->SetMarkerColor(ipt+1);
//                    mcSumPtdRsub[0][ipt]->SetFillColor(ipt+1);
//                    mcSumPtdRsub[0][ipt]->SetFillStyle(3001);
//                    mcSumPtdRsub[0][ipt]->SetLineColor(ipt+1);
//                    SumPtdRsub[0][ipt]->SetMarkerStyle(20);
//                    SumPtdRsub[0][ipt]->SetMarkerSize(1.5);
//                    SumPtdRsub[0][ipt]->SetMarkerColor(ipt+1);
//                    SumPtdRsub[0][ipt]->SetFillColor(ipt+1);
//                    SumPtdRsub[0][ipt]->SetLineColor(ipt+1);
// //               }
//                mcSumPtdRsub[0][ipt]->DrawCopy("same eHIST");
//                SumPtdRsub[0][ipt]->DrawCopy("same PE");
//                
////                t1->AddEntry(SumPtdRsub[0][ipt],Form("%.f<p_{T}^{trk}<%.f (GeV/c)", trkpt[ipt],trkpt[ipt+1]), "PL");
//                if(ipt==0){
//                    t2->AddEntry(mcSumPtdRsub[0][ipt], "PYTHIA+HYDJET", "LF2");
//                    t2->AddEntry(SumPtdRsub[0][ipt], "PbPb, #sqrt{s_{NN}} = 2.76 TeV", "PL");
//                }
//                
//                
//                
//            }
//            drawText(Form("%d-%d%%",centr[0],centr[1]),0.55,0.88, 17);
//            drawCMS(0.13,0.85,pbpbLumi);
            
            drawText("anti-k_{T} ( R =0.3 ), PF Jets",0.55,0.65,17);
            drawText(Form("p_{T}^{jet}>%.f GeV/c", pt[0]),0.65,0.6,17);
            drawText(Form(" %.1f < |#eta|^{jet} < %.f",etalimit, etacut),0.65,0.55,17);
            t2->Draw("same");
            
            if(SavePlot)
                //    c1->Print(Form("%s/%sJetPt%.f_%.fSumTrkFraction%s.eps",plotsdir.Data(),type.Data(), pt[0],pt[0+1],level.Data()));
                c1->Print(Form("%s/ppHIJetPt%.f_%.fTrkBin%dLogYRawSumPtBkgSubRadius.pdf",plotsdir.Data(),pt[0],pt[0+1],nptbin));
            c1->Update();
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

void drawCMS(float px, float py, float nLumi) {
  //  TLatex *cms = new TLatex(px,py,"CMS Preliminary");
    TLatex *cms = new TLatex(px,py,"PbPb, #sqrt{s} = 2.76 TeV");
   cms->SetTextFont(63);
    cms->SetTextSize(17);
    cms->SetNDC();
//    cms->Draw();
    TLatex *lumi = new TLatex(px,py-0.05,Form("#intL dt = %.1f #mub^{-1}",nLumi));
    lumi->SetTextFont(63);
    lumi->SetTextSize(17);
    lumi->SetNDC();
    lumi->Draw();
}

void drawCMSpp(float px, float py, float nLumi) {
    TLatex *cms = new TLatex(px,py,"pp, #sqrt{s} = 2.76 TeV");
    cms->SetTextFont(63);
    cms->SetTextSize(16);
    cms->SetNDC();
//    cms->Draw();
    TLatex *lumi = new TLatex(px,py-0.05,Form("pp, #intL dt = %.1f pb^{-1}",nLumi));
//    TLatex *lumi = new TLatex(px,py-0.05,Form("#intL dt = %.1f nb^{-1}",nLumi));
    lumi->SetTextFont(63);
    lumi->SetTextSize(16);
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
    h->SetTitleSize(22);
    h->SetLabelSize(22,"Y");
    h->SetLabelSize(22,"X");
    h->SetTitleFont(43,"Y");
    h->SetTitleSize(22,"Y");
    h->SetTitleSize(22,"X");
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
