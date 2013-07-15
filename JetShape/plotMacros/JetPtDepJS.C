/*
 *  Compare pp jet shape with smeared and unsmeared pt resolution for jet pt
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
const double etacut = 2.0 ;
const double etalimit = 0.3 ;
const double dphicut = 7*TMath::TwoPi()/8. ;
const double leadingjetcut = 100. ;
const double subleadjetcut = 40. ;
const double trackcut = 1.0;
bool SavePlot = kTRUE ;

const bool IsMC=kFALSE ;
const bool SaveFile=kFALSE ;
bool DoOneSample=kFALSE ;

const int nfile = 6 ; //8 until 0.3, 13 until 0.5
//const int nfile = (conesize-deltacone/2.)/deltacone +1 ;
//const double rad[] = {0.05, 0.15, 0.25};
//const double ratio[] = {0.05/conesize, 0.15/conesize, 0.25/conesize};
//const double errR[] = {0.05, 0.05, 0.05};
double rad[nfile];
double ratio[nfile];
double errR[nfile]; 
//const double pt[]={100., 110., 120.,130., 140., 150.,200., 300.};
//const double pt[]={110., 120., 130., 160., 200., 300., 500.};
//const double pt[]={100., 300.};
const double pt[]={100., 120., 150., 200., 300.};
const int nptbin = 4 ;

//const double pt[]={100., 120., 140.,500.};
const double subpt[]={40., 50., 60., 70., 80., 100., 120};
const double radRebin[] = {0.05, 0.15, 0.25};

const int nrbin = 6 ;
double rbin[nrbin+1]={0.0,0.05,0.10, 0.15, 0.20,0.25,0.30};

const double kAj[]={0.0,0.1,0.4,1.0};

const Double_t pbpbLumi = 150 ;
const Double_t ppLumi = 231 ;

enum Option_t {kDiff, kInt, kDiffSub, kDiffSubRcp, kRcpRatio, kDiffRatio, kSmearRatio, kIntSub} ; 
enum Method_t {kCorr, kRaw} ;
enum Level_t {kReco, kGen} ;
//---------------------------------------------------------------------

void JetPtDepJS(Option_t opt=kDiff, Level_t lev =kReco, Method_t met=kCorr, int current = 0){
    for(int i = 0 ; i <nrbin; i++){
        //        if(opt==kInt) rad[i]=rbin[i+1] ;
        //        else         
        rad[i]=deltacone/2.+i*deltacone ;
        ratio[i]=rad[i]/conesize;
        errR[i]=deltacone/2. ;
    }
    TString Norm="NormJet" ;    
    TString plotsdir = "/Users/ymao/group/CMS/plots/finalJS/CorrectDR";
    TDatime d ;
    int date = d.GetDate();  
    TString type ;
    if(IsMC) type="MC";
    else type = "DATA" ;
    TString sample ;
    if(DoOneSample)sample="OneSample";
    else sample="MergedCSdiffer";
    TString bkg ;
    if(met==kCorr) bkg = "TrkEffCorrected";
    else bkg = "NoTrkEffCorr";
    TString level ;
    if(lev==kReco) level = "RECO" ;
    else level = "GEN" ;
    TString smf ="";
        // for HI centrality bin
    
    const int nbin = 4 ;
    const int centr[] ={0,10,30,50,100}; 
    
    TString effTab ;
    if(met==kRaw) effTab = "Trk";
    else effTab = "HistIterTrkCorrtest";  //"HistIterTrkCorrv14fXSecMBEff"; // "HistIterTrkCorrv14fXSec"; // "HistIterTrkCorrtest" ; //

    if(IsMC)
        //     const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/GenJet" ;
        const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/CorrectdR/NoMatchGen/NoTrkAlgoCut" ;
    else
        //    const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/FullStatPbPb/RebinJetPt" ;
        const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/CorrectdR/NoMatchGen/NoTrkAlgoCut" ;
    TString fileName ;
    TFile * f ;

    TString ppfileName ;
    TFile * ppf ;
    
    TString smearppfileName ;
    TFile * smearppf ;

    TH2F * ppIntJS[nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * ppbkgIntJS[nrbin]; //Integrated jet shapes (psi(r)) hist.
    
    TH2F * ppDiffJS[nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * ppbkgDiffJS[nrbin]; //Integrated jet shapes (psi(r)) hist.
    
    TH1F * ppIntegrated[nrbin][nptbin] ;
    TH1F * ppbkgIntegrated[nrbin][nptbin] ;
    
    TH1F * ppdifferential[nrbin][nptbin] ;
    TH1F * ppbkgdifferential[nrbin][nptbin] ;
    
    TGraphErrors * ppRadiusRho;
    TGraphErrors * ppRadiusPsi;
    
    TGraphErrors * ppbkgRadiusRho;
    TGraphErrors * ppbkgRadiusPsi;
    
    TGraphErrors * ppRhoNoBkg;
    TGraphErrors * ppPsiNoBkg ;
    
    TH2F * ppsmearIntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * ppsmearbkgIntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    
    TH2F * ppsmearDiffJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * ppsmearbkgDiffJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    
    TH1F * ppsmearIntegrated[nbin][nrbin][nptbin] ;
    TH1F * ppsmearbkgIntegrated[nbin][nrbin][nptbin] ;
    
    TH1F * ppsmeardifferential[nbin][nrbin][nptbin] ;
    TH1F * ppsmearbkgdifferential[nbin][nrbin][nptbin] ;
    
    TGraphErrors * ppsmearRadiusRho[nbin];
    TGraphErrors * ppsmearRadiusPsi[nbin];
    
    TGraphErrors * ppsmearbkgRadiusRho[nbin];
    TGraphErrors * ppsmearbkgRadiusPsi[nbin];
    
    TGraphErrors * ppsmearRhoNoBkg[nbin];
    TGraphErrors * ppsmearPsiNoBkg[nbin];

    TH2F * IntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * bkgIntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    
    TH2F * DiffJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * bkgDiffJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    
    TH1F * Integrated[nbin][nrbin][nptbin] ;
    TH1F * bkgIntegrated[nbin][nrbin][nptbin] ;
    
    TH1F * differential[nbin][nrbin][nptbin] ;
    TH1F * bkgdifferential[nbin][nrbin][nptbin] ;
    
    TGraphErrors * RadiusRho[nbin];
    TGraphErrors * RadiusPsi[nbin];
    
    TGraphErrors * bkgRadiusRho[nbin];
    TGraphErrors * bkgRadiusPsi[nbin];
    
    TGraphErrors * RhoNoBkg[nbin];
    TGraphErrors * PsiNoBkg[nbin] ;

    TGraphErrors * RebinRhoNoBkg[nbin];
    TGraphErrors * RebinSmearRhoNoBkg[nbin];

    TGraphErrors * RatioRhoNoBkg[nbin];
    TGraphErrors * RatioPsiNoBkg[nbin] ;
    TGraphErrors * smearRatioRhoNoBkg[nbin];
    TGraphErrors * smearRatioPsiNoBkg[nbin] ;
    TGraphErrors * RcpRatioRhoNoBkg[nbin];

    TH1::SetDefaultSumw2();
    
    double rho[nbin][nrbin][nptbin];
    double errrho[nbin][nrbin][nptbin];
    double psi[nbin][nrbin][nptbin];
    double errpsi[nbin][nrbin][nptbin];
    double bkgrho[nbin][nrbin][nptbin];
    double errbkgrho[nbin][nrbin][nptbin];
    double bkgpsi[nbin][nrbin][nptbin];
    double errbkgpsi[nbin][nrbin][nptbin];
    
    double rhosub[nbin][nrbin][nptbin];
    double errrhosub[nbin][nrbin][nptbin];
    double psisub[nbin][nrbin][nptbin];
    double errpsisub[nbin][nrbin][nptbin];

    double pprho[nrbin][nptbin];
    double errpprho[nrbin][nptbin];
    double pppsi[nrbin][nptbin];
    double errpppsi[nrbin][nptbin];
    
    double ppbkgrho[nrbin][nptbin];
    double errppbkgrho[nrbin][nptbin];
    double ppbkgpsi[nrbin][nptbin];
    double errppbkgpsi[nrbin][nptbin];
    
    double pprhosub[nrbin][nptbin];
    double errpprhosub[nrbin][nptbin];
    double pppsisub[nrbin][nptbin];
    double errpppsisub[nrbin][nptbin];
    
    double ppsmearrho[nbin][nrbin][nptbin];
    double errppsmearrho[nbin][nrbin][nptbin];
    double ppsmearpsi[nbin][nrbin][nptbin];
    double errppsmearpsi[nbin][nrbin][nptbin];
    
    double ppsmearbkgrho[nbin][nrbin][nptbin];
    double errppsmearbkgrho[nbin][nrbin][nptbin];
    double ppsmearbkgpsi[nbin][nrbin][nptbin];
    double errppsmearbkgpsi[nbin][nrbin][nptbin];
    
    double ppsmearrhosub[nbin][nrbin][nptbin];
    double errppsmearrhosub[nbin][nrbin][nptbin];
    double ppsmearpsisub[nbin][nrbin][nptbin];
    double errppsmearpsisub[nbin][nrbin][nptbin];

    double rhoratio[nbin][nrbin][nptbin];
    double errrhoratio[nbin][nrbin][nptbin];
    double psiratio[nbin][nrbin][nptbin];
    double errpsiratio[nbin][nrbin][nptbin];

    double rcprhoratio[nbin][nrbin][nptbin];
    double errrcprhoratio[nbin][nrbin][nptbin];

    double smearrhoratio[nbin][nrbin][nptbin];
    double errsmearrhoratio[nbin][nrbin][nptbin];
    double smearpsiratio[nbin][nrbin][nptbin];
    double errsmearpsiratio[nbin][nrbin][nptbin];

    double sumrho[nptbin];
    double scalepsi[nbin][nptbin];
    double scalerho[nbin][nptbin];
    double scalepppsi[nptbin];
    double scalepprho[nptbin] ;
    double scaleppsmearpsi[nbin][nptbin];
    double scaleppsmearrho[nbin][nptbin];

    for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
        sumrho[ipt]=0.;
        scalepppsi[ipt]=0. ;
        scalepprho[ipt]=0.;
    }

//    TCanvas  * c1 = new TCanvas("c1", "c1", kw, kh);
//    c1->SetFillColor(0);
//    c1->SetBorderSize(0);
//    c1->SetFrameBorderMode(0); 
//    gStyle->SetOptStat(0);    

    TString canv_name = "c1";
//    canv_name += parti;
//    canv_name += cent;
    switch(opt){
        case kDiff:
            const Double_t kw = 1000;
            const Double_t kh = 350;
            c1 = new TCanvas(canv_name," ",10,10,kw,kh);
            makeMultiPanelCanvas(c1,nbin,1,0.0,0.0,0.2,0.2,0.02); 
            break ;
        case kDiffSub:
            const Double_t kw = 1000;
            const Double_t kh = 550;
//            const Double_t kw = 1300;
//            const Double_t kh = 560;
            c1 = new TCanvas(canv_name," ",10,10,kw,kh);
           makeMultiPanelCanvas(c1,nbin,2,0.0,0.0,0.2,0.2,0.15); 
            break ;
        case kDiffSubRcp:
            const Double_t kw = 1000;
            const Double_t kh = 700;
            c1 = new TCanvas(canv_name," ",10,10,kw,kh);
           makeMultiPanelCanvas(c1,nbin-1,2,0.0,0.0,0.2,0.2,0.02); 
            break ;
        case kDiffRatio:
            const Double_t kw = 1000;
            const Double_t kh = 350;
            c1 = new TCanvas(canv_name," ",10,10,kw,kh);
            makeMultiPanelCanvas(c1,nbin,1,0.0,0.0,0.2,0.2,0.02); 
            break ;     
        case kRcpRatio:
            const Double_t kw = 1000;
            const Double_t kh = 450;
            c1 = new TCanvas(canv_name," ",10,10,kw,kh);
            makeMultiPanelCanvas(c1,nbin-1,1,0.0,0.0,0.2,0.2,0.02); 
            break ;    
        case kSmearRatio:
            const Double_t kw = 1000;
            const Double_t kh = 560;
            c1 = new TCanvas(canv_name," ",10,10,kw,kh);
            makeMultiPanelCanvas(c1,nbin,2,0.0,0.0,0.2,0.2,0.02); 
            break ;
        case kCompJSRatioTh:
            const Double_t kw = 560;
            const Double_t kh = 560;
            c1 = new TCanvas(canv_name," ",10,10,kw,kh);
          //  makeMultiPanelCanvas(c1,1,1,0.0,0.0,0.2,0.2,0.02); 
            break ;
        case kIntSub:
            const Double_t kw = 1000;
            const Double_t kh = 560;
            //            const Double_t kw = 1300;
            //            const Double_t kh = 560;
            c1 = new TCanvas(canv_name," ",10,10,kw,kh);
            makeMultiPanelCanvas(c1,nbin,2,0.0,0.0,0.2,0.2,0.02); 
            break ;

    }

    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);

    
    TLegend *t1=new TLegend(0.20,0.65,0.7,0.82);
 //   TLegend *t1=new TLegend(0.25,0.6,0.8,0.92);
    t1->SetFillColor(0);
    t1->SetBorderSize(0);
    t1->SetFillStyle(0);
    t1->SetTextFont(63);
    t1->SetTextSize(17);
    TLegend *t2=new TLegend(0.20,0.25,0.35,0.4);
    t2->SetFillColor(0);
    t2->SetBorderSize(0);
    t2->SetFillStyle(0);
    t2->SetTextFont(63);
    t2->SetTextSize(17);
    // open the pp data file
    if(IsMC){
        if(lev==kGen){
            if(DoOneSample)
                ppfileName =  Form("MCPP_RefJetPt100_SimTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1pthat80_mergedFile.root", trackcut); 
else 
    ppfileName =  Form("mergedCSdiff_MCPP_RefJetPt100_SimTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_mergedFile.root", trackcut); 

        }
        else {
            switch(met){
                case kCorr:
                    if(DoOneSample)
                        ppfileName =  Form("MCPP_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3pthat80_mergedFile.root", effTab.Data(), trackcut);   
                    else 
                    ppfileName =  Form("mergedCSdiff_MCPP_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_merged.root", effTab.Data(), trackcut);           
              //      ppfileName =  Form("merged_MCPP_Ak3PFIncJetPt100_HistCent4BinCorrTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_merged.root", trackcut);              
                    break ;
                case kRaw:
                    if(DoOneSample)
                      ppfileName =  Form("MCPP_Ak3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3pthat80_mergedFile.root", trackcut); 
                    else 
                      ppfileName =  Form("mergedCSdiff_MCPP_Ak3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_merged.root", trackcut); 
                    break ;
            }
        }
    }
    else {
        switch(met){
            case kCorr:
                
                //    ppfileName =  Form("DATAPP2011_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_pp_merged_full.root", effTab.Data(), trackcut);
                ppfileName =  Form("DATAPP_Ak3PFIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut);
                break ;
            case kRaw:
                ppfileName =  Form("DATAPP_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut);
                //                ppfileName =  Form("DATAPP2011_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_pp_merged_full.root", effTab.Data(), trackcut);
                break ;
        }
    }
    if(lev==kGen) ppf = TFile::Open(Form("/Users/ymao/group/CMS/anaOutputs/JetShape/CorrectdR/%s", ppfileName.Data()), "readonly");
   // if(lev==kGen) ppf = TFile::Open(Form("/Users/ymao/group/CMS/anaOutputs/JetShape/GenJet/%s", ppfileName.Data()), "readonly");
    else {
        if(met==kRaw) ppf = TFile::Open(Form("%s/%s", kHomeDir, ppfileName.Data()), "readonly");
        else ppf = TFile::Open(Form("%s/MyselfTrkEff/%s", kHomeDir, ppfileName.Data()), "readonly");
        
    }

    cout <<"reading pp file: " << Form("%s/%s", kHomeDir, ppfileName.Data()) <<endl ;    
    for(int ir =0 ; ir <nrbin; ir++){
        ppDiffJS[ir]=(TH2F*)ppf->Get(Form("differentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,0, 100));
        ppIntJS[ir]=(TH2F*)ppf->Get(Form("IntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,0, 100));
        ppbkgIntJS[ir]=(TH2F*)ppf->Get(Form("bkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,0, 100));
        ppbkgDiffJS[ir]=(TH2F*)ppf->Get(Form("bkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,0, 100));
        ppDiffJS[ir]->Sumw2();
        ppIntJS[ir]->Sumw2();
        ppbkgDiffJS[ir]->Sumw2();
        ppbkgIntJS[ir]->Sumw2();
        for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
            pprho[ir][ipt] = 0 ;
            errpprho[ir][ipt] = 0 ;
            pppsi[ir][ipt]=0 ;
            errpppsi[ir][ipt]=0 ;
            ppbkgrho[ir][ipt] = 0 ;
            errppbkgrho[ir][ipt] = 0 ;
            ppbkgpsi[ir][ipt]=0 ;
            errppbkgpsi[ir][ipt]=0 ;
            pprhosub[ir][ipt] = 0 ;
            errpprhosub[ir][ipt] = 0 ;
            pppsisub[ir][ipt]=0 ;
            errpppsisub[ir][ipt]=0 ;
            
            ppIntegrated[ir][ipt]=(TH1F*)ppIntJS[ir]->ProjectionY(Form("ppJetPt%.f_%.fPsidR%.f_%.f_Cen%.f-%.f%%",pt[ipt], pt[ipt+1],rbin[ir]*100,rbin[ir+1]*100,0, 100),ppIntJS[ir]->GetXaxis()->FindBin(pt[ipt]), ppIntJS[ir]->GetXaxis()->FindBin(pt[ipt+1])-1, "");
            
            ppbkgIntegrated[ir][ipt]=(TH1F*)ppbkgIntJS[ir]->ProjectionY(Form("ppJetPt%.f_%.fbkgPsidR%.f_%.f_Cen%.f-%.f%%",pt[ipt], pt[ipt+1],rbin[ir]*100,rbin[ir+1]*100,0, 100),ppbkgIntJS[ir]->GetXaxis()->FindBin(pt[ipt]), ppbkgIntJS[ir]->GetXaxis()->FindBin(pt[ipt+1])-1, "");
            ppdifferential[ir][ipt]=(TH1F*)ppDiffJS[ir]->ProjectionY(Form("ppJetPt%.f_%.fRhodR%.f_%.f_Cen%.f-%.f%%",pt[ipt], pt[ipt+1],rbin[ir]*100,rbin[ir+1]*100,0, 100),ppDiffJS[ir]->GetXaxis()->FindBin(pt[ipt]), ppDiffJS[ir]->GetXaxis()->FindBin(pt[ipt+1])-1, "");
            cout << pt[ipt]<<" < Jet pt < " << pt[ipt+1] <<"Njets =" << ppdifferential[ir][ipt]->GetEntries()<<endl;
            ppbkgdifferential[ir][ipt]=(TH1F*)ppbkgDiffJS[ir]->ProjectionY(Form("ppJetPt%.f_%.fbkgRhodR%.f_%.f_Cen%.f-%.f%%",pt[ipt], pt[ipt+1],rbin[ir]*100,rbin[ir+1]*100,0, 100),ppbkgDiffJS[ir]->GetXaxis()->FindBin(pt[ipt]), ppbkgDiffJS[ir]->GetXaxis()->FindBin(pt[ipt+1])-1, "");
            
            pppsi[ir][ipt]=ppIntegrated[ir][ipt]->GetMean(1);
            errpppsi[ir][ipt]=ppIntegrated[ir][ipt]->GetMeanError(1);
            pprho[ir][ipt]=ppdifferential[ir][ipt]->GetMean(1);
            errpprho[ir][ipt]=ppdifferential[ir][ipt]->GetMeanError(1);
            ppbkgpsi[ir][ipt]=ppbkgIntegrated[ir][ipt]->GetMean(1);
            errppbkgpsi[ir][ipt]=ppbkgIntegrated[ir][ipt]->GetMeanError(1);
            ppbkgrho[ir][ipt]=ppbkgdifferential[ir][ipt]->GetMean(1);
            errppbkgrho[ir][ipt]=ppbkgdifferential[ir][ipt]->GetMeanError(1);
            //nromalize to delta cone size to consistent for definition
            pprho[ir][ipt]/=deltacone ;
            ppbkgrho[ir][ipt]/=deltacone ;
            errpprho[ir][ipt]/=deltacone ;
            errppbkgrho[ir][ipt]/=deltacone ;
            //calculate the bkg subtracted JS
            pprhosub[ir][ipt]=pprho[ir][ipt]-ppbkgrho[ir][ipt];
            errpprhosub[ir][ipt]=TMath::Sqrt(TMath::Power(errpprho[ir][ipt],2)+TMath::Power(errppbkgrho[ir][ipt],2));
            pppsisub[ir][ipt]=pppsi[ir][ipt]-ppbkgpsi[ir][ipt];
            errpppsisub[ir][ipt]=TMath::Sqrt(TMath::Power(errpppsi[ir][ipt],2)+TMath::Power(errppbkgpsi[ir][ipt],2));
            cout <<"ir = "<<ir<<"psi =" << pppsisub[ir][ipt] <<"rho =" << pprhosub[ir][ipt]*deltacone<<endl ;
            //        pprhosub[ir]=pprho[ir];
            //        errpprhosub[ir]=TMath::Sqrt(TMath::Power(errpprho[ir],2));
            //        pppsisub[ir]=pppsi[ir]-ppbkgpsi[ir];
            //        errpppsisub[ir]=TMath::Sqrt(TMath::Power(errpppsi[ir],2)+TMath::Power(errppbkgpsi[ir],2));
            
     //       scalepprho[ipt]+=pprhosub[ir][ipt]*deltacone;
            //    scaleppsmearrho+=ppsmearrhosub[ir]*deltacone;
        }
    }
    for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
        scalepppsi[ipt]=0. ;
        scalepprho[ipt]=0.;
        scalepprho[ipt]=(pprhosub[0][ipt]+pprhosub[1][ipt]+pprhosub[2][ipt]+pprhosub[3][ipt]+pprhosub[4][ipt]+pprhosub[5][ipt])*deltacone;
        scalepppsi[ipt]=pppsisub[5][ipt];
        cout <<" pp scale =" << scalepprho[ipt] <<endl ;
    }
    
        //rescale bkg subtracted JS in order to get unity for psi
        for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
            sumrho[ipt]=0.;
            for(int ir =0 ; ir <nrbin; ir++){
            pppsisub[ir][ipt]/=scalepppsi[ipt] ;
            pprhosub[ir][ipt]/=scalepprho[ipt] ;
            errpppsisub[ir][ipt]/=scalepppsi[ipt] ;
            errpprhosub[ir][ipt]/=scalepprho[ipt] ;
            sumrho[ipt]+=pprhosub[ir][ipt];
            
//            for(int i =0 ; i <ir; i++){
//                pppsisub[ir][ipt]+=pprhosub[i][ipt]*deltacone ;
//                errpppsisub[ir][ipt]+=TMath::Power(errpprhosub[i][ipt]*deltacone,2);
//            }
//            pppsisub[ir][ipt]= 1- pppsisub[ir][ipt] ;
//            errpppsisub[ir][ipt]=TMath::Sqrt(errpppsisub[ir][ipt]);
            
        }
            cout <<" pp rho =" <<sumrho[ipt] <<" psi =" << pppsisub[5][ipt]<<endl ;
    }
    
        ppRadiusRho = new TGraphErrors(nptbin, pt, pprho[current], 0, errpprho[current]);
        ppRadiusRho->SetMarkerStyle(24);
        ppRadiusRho->SetMarkerColor(2);
        ppRadiusRho->SetMarkerSize(1.5);
        ppRadiusRho->SetLineColor(2);
        ppbkgRadiusRho = new TGraphErrors(nptbin, pt,  ppbkgrho[current], 0, errppbkgrho[current]);
        ppbkgRadiusRho->SetMarkerStyle(30);
        ppbkgRadiusRho->SetMarkerSize(1.5);
        ppbkgRadiusRho->SetMarkerColor(2);
        ppbkgRadiusRho->SetLineColor(2);
        
        ppRadiusPsi = new TGraphErrors(nptbin, pt,  pppsi[current], 0, errpppsi[current]);
        ppRadiusPsi->SetMarkerStyle(21);
        ppRadiusPsi->SetMarkerColor(4);
        ppRadiusPsi->SetLineColor(4);
        ppbkgRadiusPsi = new TGraphErrors(nptbin, pt,  ppbkgpsi[current], 0, errppbkgpsi[current]);
        ppbkgRadiusPsi->SetMarkerStyle(25);
        ppbkgRadiusPsi->SetMarkerColor(4);
        ppbkgRadiusPsi->SetLineColor(4);
        
        ppRhoNoBkg = new TGraphErrors(nptbin, pt,  pprhosub[current], 0, errpprhosub[current]);
        ppRhoNoBkg->SetMarkerStyle(25);
        ppRhoNoBkg->SetMarkerColor(2);
         ppRhoNoBkg->SetMarkerSize(1.5);
        ppRhoNoBkg->SetLineColor(2);
        ppPsiNoBkg = new TGraphErrors(nptbin, pt,  pppsisub[current], 0, errpppsisub[current]);
        ppPsiNoBkg->SetMarkerStyle(25);
        ppPsiNoBkg->SetMarkerColor(1);    
        ppPsiNoBkg->SetMarkerSize(1.5);
        ppPsiNoBkg->SetLineColor(1);        

    // open the HI data file
    if(IsMC){
        if(lev==kGen){
            if(DoOneSample)
                fileName =  Form("MCHI_RefJetPt100_SimTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1pthat80_Dijet80_HydjetDrum_v27_mergedV1.root", trackcut); 
else 
    fileName =  Form("mergedCSdiff_MCHI_RefJetPt100_SimTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1_Dijet_HydjetDrum_v27_mergedV1.root", trackcut); 
        }
        else {
            switch(met){
                case kCorr:
                    if(DoOneSample){
                        smearppfileName = Form("MCPP_Ak3PFwtIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1_trkbin3pthat80_mergedFile.root", effTab.Data(), trackcut);
                        fileName=Form("MCHI_AkPu3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin3_Dijet80_HydjetDrum_v27_mergedV1.root", effTab.Data(), trackcut, nbin);
                    }
                    else {   
                    fileName = Form("mergedCSdiff_MCHI_AkPu3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_Dijet_HydjetDrum_v27v28.root", effTab.Data(), trackcut);
                    smearppfileName=Form("mergedCSdiff_MCPP_Ak3PFwtIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_trkbin3_pp276Dijet_merged.root", effTab.Data(), trackcut, nbin);
                    }
                    break ;
                case kRaw:
                    if(DoOneSample){
                        smearppfileName = Form("MCPP_Ak3PFwtIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3pthat80_mergedFile.root", trackcut);
                        fileName=Form("MCHI_AkPu3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin3_Dijet80_HydjetDrum_v27_mergedV1.root", trackcut, nbin); 
                    }
                    else {
                    fileName=Form("mergedCSdiff_MCHI_AkPu3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_Dijet_HydjetDrum_v27v28.root", trackcut);
                    smearppfileName=Form("mergedCSdiff_MCPP_Ak3PFwtIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1_trkbin3_pp276Dijet_merged.root", trackcut);
                    }
                    break ;
            }
        }
    }
    else {
        switch(met){
            case kCorr:
                fileName=Form("DATAHI_AkPu3PFIncJetPt100_%s%.fEtaCut%.fLimit%.f_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin3_promptskim-hihighpt-hltjet80-pt90-v20.root", effTab.Data(), trackcut, etacut*10, etalimit*10, nbin);
                
                //       smearppfileName=Form("DATAPP_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin3_pp_merged_full.root", effTab.Data(), trackcut, nbin);
                //        smearppfileName=Form("DATAPP_Ak3PFwtIncJetYaxianSMParaPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin3_pp_merged_full.root", effTab.Data(), trackcut, nbin);
                
                smearppfileName=Form("DATAPP_Ak3PFwtIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut);
                //       smearppfileName=Form("mergedCSdiff_MCPP_Ak3PFwtIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_trkbin3_pp276Dijet_merged.root", effTab.Data(), trackcut, nbin);
                //      smearppfileName=Form("DATAPP_Ak3PFIncJetPt100_%s%.fEtaCut%.fLimit%.f_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82.root", effTab.Data(), trackcut, etacut*10, etalimit*10, nbin);
                break ;
            case kRaw:
                fileName=Form("DATAHI_AkPu3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_promptskim-hihighpt-hltjet80-pt90-v20.root", trackcut);
                smearppfileName=Form("DATAPP_Ak3PFwtIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_JetDB_forestv78.root", trackcut);
//                smearppfileName=Form("DATAPP_Ak3PFwtIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_pp_merged_full.root", trackcut);
                break ;
        }  
    }
    if(IsMC){
        if(lev==kGen) f = TFile::Open(Form("/Users/ymao/group/CMS/anaOutputs/JetShape/CorrectdR/%s", fileName.Data()), "readonly");
//        if(lev==kGen) f = TFile::Open(Form("/Users/ymao/group/CMS/anaOutputs/JetShape/GenJet/%s", fileName.Data()), "readonly");
        else {
            f = TFile::Open(Form("%s/%s", kHomeDir, fileName.Data()), "readonly");
            if(lev==kReco){
                if(met==kRaw) smearppf = TFile::Open(Form("%s/%s", kHomeDir, smearppfileName.Data()), "readonly");
                else smearppf = TFile::Open(Form("%s/MyselfTrkEff/%s", kHomeDir, smearppfileName.Data()), "readonly");
            }
        }
    }
    else {
        f = TFile::Open(Form("%s/%s", kHomeDir, fileName.Data()), "readonly");
        if(lev==kReco){
            if(met==kRaw) smearppf = TFile::Open(Form("%s/%s", kHomeDir, smearppfileName.Data()), "readonly");
            else smearppf = TFile::Open(Form("%s/MyselfTrkEff/%s", kHomeDir, smearppfileName.Data()), "readonly");
        }
    }
    cout <<"reading ppsmeared file: " << Form("%s/%s", kHomeDir, smearppfileName.Data()) <<endl ;
    cout <<"reading HI file: " << Form("%s/%s", kHomeDir, fileName.Data()) <<endl ;    

        for(int ibin = 0 ; ibin <nbin; ibin++){
            for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
                scalepsi[ibin][ipt]=0. ;
                scalerho[ibin][ipt]=0. ;
                scaleppsmearpsi[ibin][ipt]=0. ;
                scaleppsmearrho[ibin][ipt]=0. ;
            }
        }
    for(int ibin = 0 ; ibin <nbin; ibin++){
    //    cout <<" Centality =" << centr[ibin] << " - " << centr[ibin+1] <<"% "<<endl;
        for(int ir =0 ; ir <nrbin; ir++){
            IntJS[ibin][ir]=(TH2F*)f->Get(Form("IntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));
            bkgIntJS[ibin][ir]=(TH2F*)f->Get(Form("bkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));
            
            DiffJS[ibin][ir]=(TH2F*)f->Get(Form("differentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));
            bkgDiffJS[ibin][ir]=(TH2F*)f->Get(Form("bkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));
            
            DiffJS[ibin][ir]->Sumw2();
            IntJS[ibin][ir]->Sumw2();
            bkgIntJS[ibin][ir]->Sumw2();
            bkgDiffJS[ibin][ir]->Sumw2();
            
            for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
                rho[ibin][ir][ipt] = 0 ;
                errrho[ibin][ir][ipt] = 0 ;
                psi[ibin][ir][ipt]=0 ;
                errpsi[ibin][ir][ipt]=0 ;
                bkgrho[ibin][ir][ipt] = 0 ;
                errbkgrho[ibin][ir][ipt] = 0 ;
                bkgpsi[ibin][ir][ipt]=0 ;
                errbkgpsi[ibin][ir][ipt]=0 ;
                
                rhosub[ibin][ir][ipt] = 0 ;
                errrhosub[ibin][ir][ipt] = 0 ;
                psisub[ibin][ir][ipt]=0 ;
                errpsisub[ibin][ir][ipt]=0 ;
                
                rhoratio[ibin][ir][ipt] = 0 ;
                errrhoratio[ibin][ir][ipt] = 0 ;
                psiratio[ibin][ir][ipt]=0 ;
                errpsiratio[ibin][ir][ipt]=0 ;
                
                rcprhoratio[ibin][ir][ipt] = 0 ;
                errrcprhoratio[ibin][ir][ipt] = 0 ;
                
                Integrated[ibin][ir][ipt]=(TH1F*)IntJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fPsidR%.f_%.f_Cen%.f-%.f%%",pt[ipt], pt[ipt+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),IntJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt]), IntJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt+1])-1, "");
                
                bkgIntegrated[ibin][ir][ipt]=(TH1F*)bkgIntJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fbkgPsidR%.f_%.f_Cen%.f-%.f%%",pt[ipt], pt[ipt+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),bkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt]), bkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt+1])-1, "");
                differential[ibin][ir][ipt]=(TH1F*)DiffJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fRhodR%.f_%.f_Cen%.f-%.f%%",pt[ipt], pt[ipt+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),DiffJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt]), DiffJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt+1])-1, "");
                bkgdifferential[ibin][ir][ipt]=(TH1F*)bkgDiffJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fbkgRhodR%.f_%.f_Cen%.f-%.f%%",pt[ipt], pt[ipt+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),bkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt]), bkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt+1])-1, "");
                
                
           //     cout <<" Centality =" << centr[ibin] << " - " << centr[ibin+1] <<"% "<< pt[ipt]<<" < Jet pt < " << pt[ipt+1] <<"Njets =" << differential[ibin][ir][ipt]->GetEntries()<<endl;
                
                psi[ibin][ir][ipt]= Integrated[ibin][ir][ipt]->GetMean(1);
                errpsi[ibin][ir][ipt]= Integrated[ibin][ir][ipt]->GetMeanError(1);
                bkgpsi[ibin][ir][ipt]= bkgIntegrated[ibin][ir][ipt]->GetMean(1);
                errbkgpsi[ibin][ir][ipt]= bkgIntegrated[ibin][ir][ipt]->GetMeanError(1);
                rho[ibin][ir][ipt]= differential[ibin][ir][ipt]->GetMean(1);
                errrho[ibin][ir][ipt]= differential[ibin][ir][ipt]->GetMeanError(1);
                bkgrho[ibin][ir][ipt]= bkgdifferential[ibin][ir][ipt]->GetMean(1);
                errbkgrho[ibin][ir][ipt]= bkgdifferential[ibin][ir][ipt]->GetMeanError(1);
                
                //nromalize to delta cone size to consistent for definition
                rho[ibin][ir][ipt]/=deltacone ;
                bkgrho[ibin][ir][ipt]/=deltacone ;
                errrho[ibin][ir][ipt]/=deltacone ;
                errbkgrho[ibin][ir][ipt]/=deltacone ;
                
                //calculate the bkg subtracted JS
                rhosub[ibin][ir][ipt]=rho[ibin][ir][ipt]-bkgrho[ibin][ir][ipt];
                //   errrhosub[ibin][ir]=errrho[ibin][ir]+errbkgrho[ibin][ir];
                errrhosub[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errrho[ibin][ir][ipt],2)+TMath::Power(errbkgrho[ibin][ir][ipt],2));
                // cout <<"ir = "<<ir<<"rho =" << rhosub[ibin][ir][ipt]<<endl ;
                //     cout <<"ir = "<<ir<<"rho =" << rhosub[ibin][ir][ipt]<<"mean =" << differential[ibin][ir][ipt]->GetMean(1)<<" bkg mean =" << bkgdifferential[ibin][ir][ipt]->GetMean(1)<<endl ;
                
                //                rhosub[ibin][ir]=rho[ibin][ir];
                //             //   errrhosub[ibin][ir]=errrho[ibin][ir]+errbkgrho[ibin][ir];
                //                  errrhosub[ibin][ir]=TMath::Sqrt(TMath::Power(errrho[ibin][ir],2));
                //      errrhosub[ibin][ir]*=rhosub[ibin][ir];
                
                psisub[ibin][ir][ipt]=psi[ibin][ir][ipt]-bkgpsi[ibin][ir][ipt];
                errpsisub[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errpsi[ibin][ir][ipt],2)+TMath::Power(errbkgpsi[ibin][ir][ipt],2));
            }  // jet pt bin loop
            if(lev==kReco){
                ppsmearDiffJS[ibin][ir]=(TH2F*)smearppf->Get(Form("differentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));
                ppsmearIntJS[ibin][ir]=(TH2F*)smearppf->Get(Form("IntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));
                ppsmearbkgIntJS[ibin][ir]=(TH2F*)smearppf->Get(Form("bkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));
                ppsmearbkgDiffJS[ibin][ir]=(TH2F*)smearppf->Get(Form("bkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));
                ppsmearDiffJS[ibin][ir]->Sumw2();
                ppsmearIntJS[ibin][ir]->Sumw2();
                ppsmearbkgDiffJS[ibin][ir]->Sumw2();
                ppsmearbkgIntJS[ibin][ir]->Sumw2();
                for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
                    ppsmearrho[ibin][ir][ipt] = 0 ;
                    errppsmearrho[ibin][ir][ipt] = 0 ;
                    ppsmearpsi[ibin][ir][ipt]=0 ;
                    errppsmearpsi[ibin][ir][ipt]=0 ;
                    ppsmearbkgrho[ibin][ir][ipt] = 0 ;
                    errppsmearbkgrho[ibin][ir][ipt] = 0 ;
                    ppsmearbkgpsi[ibin][ir][ipt]=0 ;
                    errppsmearbkgpsi[ibin][ir][ipt]=0 ;
                    ppsmearrhosub[ibin][ir][ipt] = 0 ;
                    errppsmearrhosub[ibin][ir][ipt] = 0 ;
                    ppsmearpsisub[ibin][ir][ipt]=0 ;
                    errppsmearpsisub[ibin][ir][ipt]=0 ;
                    
                    ppsmearIntegrated[ibin][ir][ipt]=(TH1F*)ppsmearIntJS[ibin][ir]->ProjectionY(Form("ppsmearJetPt%.f_%.fPsidR%.f_%.f_Cen%.f-%.f%%",pt[ipt], pt[ipt+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),ppsmearIntJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt]), ppsmearIntJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt+1])-1, "");
                    
                    ppsmearbkgIntegrated[ibin][ir][ipt]=(TH1F*)ppsmearbkgIntJS[ibin][ir]->ProjectionY(Form("ppsmearJetPt%.f_%.fbkgPsidR%.f_%.f_Cen%.f-%.f%%",pt[ipt], pt[ipt+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),ppsmearbkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt]), ppsmearbkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt+1])-1, "");
                    
                    ppsmeardifferential[ibin][ir][ipt]=(TH1F*)ppsmearDiffJS[ibin][ir]->ProjectionY(Form("ppsmearJetPt%.f_%.fRhodR%.f_%.f_Cen%.f-%.f%%",pt[ipt], pt[ipt+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),ppsmearDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt]), ppsmearDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt+1])-1, "");
                    
                    ppsmearbkgdifferential[ibin][ir][ipt]=(TH1F*)ppsmearbkgDiffJS[ibin][ir]->ProjectionY(Form("ppsmearJetPt%.f_%.fbkgRhodR%.f_%.f_Cen%.f-%.f%%",pt[ipt], pt[ipt+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),ppsmearbkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt]), ppsmearbkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[ipt+1])-1, "");
                    
                    //     cout << pt[ipt]<<" < Jet pt < " << pt[ipt+1] <<"   smeared pp Njets =" << ppsmeardifferential[ibin][ir][ipt]->GetEntries()<<endl;
                    
                    ppsmearpsi[ibin][ir][ipt]=ppsmearIntegrated[ibin][ir][ipt]->GetMean(1);
                    errppsmearpsi[ibin][ir][ipt]=ppsmearIntegrated[ibin][ir][ipt]->GetMeanError(1);
                    ppsmearrho[ibin][ir][ipt]=ppsmeardifferential[ibin][ir][ipt]->GetMean(1);
                    errppsmearrho[ibin][ir][ipt]=ppsmeardifferential[ibin][ir][ipt]->GetMeanError(1);
                    ppsmearbkgpsi[ibin][ir][ipt]=ppsmearbkgIntegrated[ibin][ir][ipt]->GetMean(1);
                    errppsmearbkgpsi[ibin][ir][ipt]=ppsmearbkgIntegrated[ibin][ir][ipt]->GetMeanError(1);
                    ppsmearbkgrho[ibin][ir][ipt]=ppsmearbkgdifferential[ibin][ir][ipt]->GetMean(1);
                    errppsmearbkgrho[ibin][ir][ipt]=ppsmearbkgdifferential[ibin][ir][ipt]->GetMeanError(1);
                    ppsmearrho[ibin][ir][ipt]/=deltacone ;
                    ppsmearbkgrho[ibin][ir][ipt]/=deltacone ;
                    errppsmearrho[ibin][ir][ipt]/=deltacone ;
                    errppsmearbkgrho[ibin][ir][ipt]/=deltacone ;
                    
                    //                ppsmearrhosub[ibin][ir]=ppsmearrho[ibin][ir];
                    //                errppsmearrhosub[ibin][ir]=TMath::Sqrt(TMath::Power(errppsmearrho[ibin][ir],2));
                    ppsmearrhosub[ibin][ir][ipt]=ppsmearrho[ibin][ir][ipt]-ppsmearbkgrho[ibin][ir][ipt];
                    errppsmearrhosub[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errppsmearrho[ibin][ir][ipt],2)+TMath::Power(errppsmearbkgrho[ibin][ir][ipt],2));
                    ppsmearpsisub[ibin][ir][ipt]=ppsmearpsi[ibin][ir][ipt]-ppsmearbkgpsi[ibin][ir][ipt];
                    errppsmearpsisub[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errppsmearpsi[ibin][ir][ipt],2)+TMath::Power(errppsmearbkgpsi[ibin][ir][ipt],2));
                    cout <<"ir = "<<ir<<" smeared rho =" << ppsmearrhosub[ibin][ir][ipt]<<endl ;
                    
                    smearrhoratio[ibin][ir][ipt] = 0 ;
                    errsmearrhoratio[ibin][ir][ipt] = 0 ;
                    smearpsiratio[ibin][ir][ipt]=0 ;
                    errsmearpsiratio[ibin][ir][ipt]=0 ;
                }  //jet pt bin loop
                
            }  //for smeared JS, only at RECO level
            
            
        }  //radius loop

            for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
                scalerho[ibin][ipt]=(rhosub[ibin][0][ipt]+rhosub[ibin][1][ipt]+rhosub[ibin][2][ipt]+rhosub[ibin][3][ipt]+rhosub[ibin][4][ipt]+rhosub[ibin][5][ipt])*deltacone;
                scalepsi[ibin][ipt]=psisub[ibin][5][ipt];
                if(lev==kReco){
                   scaleppsmearrho[ibin][ipt]=(ppsmearrhosub[ibin][0][ipt]+ppsmearrhosub[ibin][1][ipt]+ppsmearrhosub[ibin][2][ipt]+ppsmearrhosub[ibin][3][ipt]+ppsmearrhosub[ibin][4][ipt]+ppsmearrhosub[ibin][5][ipt])*deltacone;
                    scaleppsmearpsi[ibin][ipt]=ppsmearpsisub[ibin][5][ipt];
                }
        cout <<"ibin = "<<ibin<<"scale =" << scalerho[ibin][ipt] <<"ppsmeared =" << scaleppsmearrho[ibin][ipt]<<endl ;
            }
        
            for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
                sumrho[ipt] = 0.0 ;
                for(int ir =0 ; ir <nrbin; ir++){
                //     if(ir>0){
                if(lev==kReco){
                    ppsmearpsisub[ibin][ir][ipt]/=scaleppsmearpsi[ibin][ipt] ;
                    ppsmearrhosub[ibin][ir][ipt]/=scaleppsmearrho[ibin][ipt] ;
                    errppsmearpsisub[ibin][ir][ipt]/=scaleppsmearpsi[ibin][ipt] ;
                    errppsmearrhosub[ibin][ir][ipt]/=scaleppsmearrho[ibin][ipt] ;
                }
                psisub[ibin][ir][ipt]/=scalepsi[ibin][ipt] ;
                rhosub[ibin][ir][ipt]/=scalerho[ibin][ipt] ;
                errpsisub[ibin][ir][ipt]/=scalepsi[ibin][ipt] ;
                errrhosub[ibin][ir][ipt]/=scalerho[ibin][ipt] ;
                //    }
                
                sumrho[ipt]+=rhosub[ibin][ir][ipt];
                
//                // try to calculate 1-psi instead of psi from rho
//                if(lev==kReco){
//                    for(int i =0 ; i <ir; i++){
//                        ppsmearpsisub[ibin][ir][ipt]+=ppsmearrhosub[ibin][i][ipt]*deltacone ;
//                        errppsmearpsisub[ibin][ir][ipt]+=TMath::Power(errppsmearrhosub[ibin][i][ipt]*deltacone,2);
//                    }
//                    ppsmearpsisub[ibin][ir][ipt]= 1- ppsmearpsisub[ibin][ir][ipt] ;
//                    errppsmearpsisub[ibin][ir][ipt]=TMath::Sqrt(errppsmearpsisub[ibin][ir][ipt]);
//                }
//                for(int i =0 ; i <ir; i++){
//                    psisub[ibin][ir][ipt]+=rhosub[ibin][i][ipt]*deltacone ;
//                    errpsisub[ibin][ir][ipt]+=TMath::Power(errrhosub[ibin][i][ipt]*deltacone,2);
//                }
//                psisub[ibin][ir][ipt]= 1- psisub[ibin][ir][ipt] ;
//                errpsisub[ibin][ir][ipt]=TMath::Sqrt(errpsisub[ibin][ir][ipt]);
                
                //           cout <<"ir = "<<ir<<"psi =" << psisub[ibin][ir] <<"rho =" << rhosub[ibin][ir]<<endl ;
                
                if(IsMC==kTRUE){
                    if(lev==kGen){
                        rhoratio[ibin][ir][ipt]=rhosub[ibin][ir][ipt]/pprhosub[ir][ipt];
                        errrhoratio[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errrhosub[ibin][ir][ipt]/rhosub[ibin][ir][ipt],2)+TMath::Power(errpprhosub[ir][ipt]/pprhosub[ir][ipt],2));
                        errrhoratio[ibin][ir][ipt]*=rhoratio[ibin][ir][ipt];
                        if(ir==1) cout <<"pprho ="<<pprhosub[ir][ipt] <<"  HI rho ="<<rhosub[ibin][ir][ipt]<<endl;
                        psiratio[ibin][ir][ipt]=psisub[ibin][ir][ipt]/pppsisub[ir][ipt];
                        errpsiratio[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errpsisub[ibin][ir][ipt]/psisub[ibin][ir][ipt],2)+TMath::Power(errpppsisub[ir][ipt]/pppsisub[ir][ipt],2));
                        errpsiratio[ibin][ir][ipt]*=psiratio[ibin][ir][ipt];
                    }
                    else {
                        rhoratio[ibin][ir][ipt]=rhosub[ibin][ir][ipt]/ppsmearrhosub[ibin][ir][ipt];
                        errrhoratio[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errrhosub[ibin][ir][ipt]/rhosub[ibin][ir][ipt],2)+TMath::Power(errppsmearrhosub[ibin][ir][ipt]/ppsmearrhosub[ibin][ir][ipt],2));
                        errrhoratio[ibin][ir][ipt]*=rhoratio[ibin][ir][ipt];
                        
                        psiratio[ibin][ir][ipt]=psisub[ibin][ir][ipt]/ppsmearpsisub[ibin][ir][ipt];
                        errpsiratio[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errpsisub[ibin][ir][ipt]/psisub[ibin][ir][ipt],2)+TMath::Power(errppsmearpsisub[ibin][ir][ipt]/ppsmearpsisub[ibin][ir][ipt],2));
                        errpsiratio[ibin][ir][ipt]*=psiratio[ibin][ir][ipt];
                        
                        smearrhoratio[ibin][ir][ipt]=ppsmearrhosub[ibin][ir][ipt]/pprhosub[ir][ipt];
                        //   errsmearrhoratio[ibin][ir]=errppsmearrhosub[ibin][ir];
                        errsmearrhoratio[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errppsmearrhosub[ibin][ir][ipt]/ppsmearrhosub[ibin][ir][ipt],2));
                        errsmearrhoratio[ibin][ir][ipt]*=smearrhoratio[ibin][ir][ipt];
                        //             errsmearrhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errppsmearrhosub[ibin][ir],2)+TMath::Power(errpprhosub[ir],2));
                        smearpsiratio[ibin][ir][ipt]=ppsmearpsisub[ibin][ir][ipt]/pppsisub[ir][ipt];
                        //   errsmearpsiratio[ibin][ir]=errppsmearpsisub[ibin][ir];
                        errsmearpsiratio[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errppsmearpsisub[ibin][ir][ipt]/ppsmearpsisub[ibin][ir][ipt],2));
                        errsmearpsiratio[ibin][ir][ipt]*=smearpsiratio[ibin][ir][ipt];
                    }
                }
                else {
                    rhoratio[ibin][ir][ipt]=rhosub[ibin][ir][ipt]/ppsmearrhosub[ibin][ir][ipt];
                    errrhoratio[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errrhosub[ibin][ir][ipt]/rhosub[ibin][ir][ipt],2)+TMath::Power(errppsmearrhosub[ibin][ir][ipt]/ppsmearrhosub[ibin][ir][ipt],2));
                    errrhoratio[ibin][ir][ipt]*=rhoratio[ibin][ir][ipt];
                    
                    psiratio[ibin][ir][ipt]=psisub[ibin][ir][ipt]/ppsmearpsisub[ibin][ir][ipt];
                    errpsiratio[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errpsisub[ibin][ir][ipt]/psisub[ibin][ir][ipt],2)+TMath::Power(errppsmearpsisub[ibin][ir][ipt]/ppsmearpsisub[ibin][ir][ipt],2));
                    errpsiratio[ibin][ir][ipt]*=psiratio[ibin][ir][ipt];
                    
                    smearrhoratio[ibin][ir][ipt]=ppsmearrhosub[ibin][ir][ipt]/pprhosub[ir][ipt];
                    //   errsmearrhoratio[ibin][ir]=errppsmearrhosub[ibin][ir];
                    errsmearrhoratio[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errppsmearrhosub[ibin][ir][ipt]/ppsmearrhosub[ibin][ir][ipt],2)+TMath::Power(errpprhosub[ir][ipt]/pprhosub[ir][ipt],2));
                    errsmearrhoratio[ibin][ir][ipt]*=smearrhoratio[ibin][ir][ipt];

                    smearpsiratio[ibin][ir][ipt]=ppsmearpsisub[ibin][ir][ipt]/pppsisub[ir][ipt];
                    //   errsmearpsiratio[ibin][ir]=errppsmearpsisub[ibin][ir];
                    errsmearpsiratio[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errppsmearpsisub[ibin][ir][ipt]/ppsmearpsisub[ibin][ir][ipt],2));
                    errsmearpsiratio[ibin][ir][ipt]*=smearpsiratio[ibin][ir][ipt];  
                    
                }
                
                //   cout <<"ir = "<<ir<<"rho ratio =" << rhoratio[ibin][ir] <<"psi ratio =" << psiratio[ibin][ir]<<endl ;  
            } //radius loop
                cout <<"sum rho =" <<sumrho[ipt] << "psi =" << psisub[ibin][5][ipt]<<endl ;
        }  //jet pt loop
    }
    for(int ibin = 0 ; ibin <nbin; ibin++){
        
        for(int ir =0 ; ir <nrbin; ir++){
            for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
                
                rcprhoratio[ibin][ir][ipt]=rhosub[ibin][ir][ipt]/rhosub[nbin-1][ir][ipt];
                errrcprhoratio[ibin][ir][ipt]=TMath::Sqrt(TMath::Power(errrhosub[ibin][ir][ipt]/rhosub[ibin][ir][ipt],2));
                errrcprhoratio[ibin][ir][ipt]*=rcprhoratio[ibin][ir][ipt];
                
                //  cout <<"ibin = "<<ibin <<"ir = "<<ir <<"rhoratio =" << rhoratio[ibin][ir]<<endl ;
            }
        }
    }
    TH1F * dummy = new TH1F("dummy", "dummy", 500, 0., 500.);
    dummy->SetTitle("") ;
    dummy->SetStats(kFALSE) ;
           dummy->SetAxisRange(80., 220., "X") ;

    dummy->GetXaxis()->SetTitle(" p_{T} (GeV/c) ");

    if(opt==kDiff|| opt==kDiffRatio || opt==kRcpRatio)fixedFontHist(dummy,1.1, 1.3);
    else fixedFontHist(dummy,1.8, 2.0);

    switch (opt) {
        case kDiff:
        //    if(nbin>1)c1->Divide(nbin/2, 2,0, 0);
            //  else c1->Divide(1, 2,0, 0);           
            for(int ipad = 1 ; ipad <=nbin; ipad++){
                c1->cd(ipad);
                c1->cd(ipad)->SetLogy();
                if(trackcut>2.)dummy->SetAxisRange(1.e-3, 20., "Y") ;
                else 
                dummy->SetAxisRange(1.e-4, 100., "Y") ;
                dummy->DrawCopy();
                
                RadiusRho[nbin-ipad] = new TGraphErrors(nptbin, pt, rho[nbin-ipad][current], 0, errrho[nbin-ipad][current]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                RadiusRho[nbin-ipad]->SetMarkerStyle(20);
                RadiusRho[nbin-ipad]->SetMarkerColor(1);
                RadiusRho[nbin-ipad]->SetMarkerSize(1.5);
                RadiusRho[nbin-ipad]->SetLineColor(1);
                RadiusRho[nbin-ipad]->Draw("same PE") ;
                
                bkgRadiusRho[nbin-ipad] = new TGraphErrors(nptbin, pt, bkgrho[nbin-ipad][current], 0, errbkgrho[nbin-ipad][current]);
                bkgRadiusRho[nbin-ipad]->SetMarkerStyle(29);
                bkgRadiusRho[nbin-ipad]->SetMarkerColor(1);
                bkgRadiusRho[nbin-ipad]->SetMarkerSize(1.5);
                bkgRadiusRho[nbin-ipad]->SetLineColor(1);
                bkgRadiusRho[nbin-ipad]->Draw("same PE") ;
                if(lev==kReco){
                ppsmearRadiusRho[nbin-ipad] = new TGraphErrors(nptbin, pt,ppsmearrho[nbin-ipad][current], 0, errppsmearrho[nbin-ipad][current]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                ppsmearRadiusRho[nbin-ipad]->SetMarkerStyle(24);
                ppsmearRadiusRho[nbin-ipad]->SetMarkerColor(2);
                ppsmearRadiusRho[nbin-ipad]->SetMarkerSize(1.5);
                ppsmearRadiusRho[nbin-ipad]->SetLineColor(2);
                ppsmearRadiusRho[nbin-ipad]->Draw("same PE") ;
                
                ppsmearbkgRadiusRho[nbin-ipad] = new TGraphErrors(nptbin, pt, ppsmearbkgrho[nbin-ipad][current], 0, errppsmearbkgrho[nbin-ipad][current]);
                ppsmearbkgRadiusRho[nbin-ipad]->SetMarkerStyle(30);
                ppsmearbkgRadiusRho[nbin-ipad]->SetMarkerSize(1.5);
                ppsmearbkgRadiusRho[nbin-ipad]->SetMarkerColor(2);
                ppsmearbkgRadiusRho[nbin-ipad]->SetLineColor(2);
                ppsmearbkgRadiusRho[nbin-ipad]->Draw("same PE") ;
                }
                else {
                    ppRadiusRho->Draw("same PE") ;
                    ppbkgRadiusRho->Draw("same PE") ;
                }
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[nbin-ipad+1]),0.2,0.92, 15);
                if(ipad==2){
                    if(IsMC){
                        t1->AddEntry(RadiusRho[nbin-1],"PYTHIA+HYDJET", "PL");
                        //      t1->AddEntry(ppRhoNoBkg,"pp unsmeard", "P");  
                        t1->AddEntry(bkgRadiusRho[nbin-1],"PYTHIA+HYDJET:bkg", "PL");
                        if(lev==kReco){
                        t1->AddEntry(ppsmearRadiusRho[nbin-1],"PYTHIA", "PL");   
                        //      t1->AddEntry(ppRhoNoBkg,"pp unsmeard", "P");  
                        t1->AddEntry(ppsmearbkgRadiusRho[nbin-1],"PYTHIA:bkg", "PL"); 
                        }
                        else {
                            t1->AddEntry(ppRadiusRho,"PYTHIA", "PL");   
                            //      t1->AddEntry(ppRhoNoBkg,"pp unsmeard", "P");  
                            t1->AddEntry(ppbkgRadiusRho,"PYTHIA:bkg", "PL"); 
 
                        }
                    }
                    else {
                        t1->AddEntry(RadiusRho[nbin-1],"PbPb #sqrt{s}=2.76 TeV", "PL");
                        //    t1->AddEntry(ppRhoNoBkg,"pp measured", "PL");  
                        t1->AddEntry(ppsmearRadiusRho[nbin-1],"pp reference", "PL");
                        //        t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"PYTHIA reference", "PL");
                    }
                    t1->Draw("same");
                }
               // if(ipad==1)drawText(Form("%s jet",level.Data()),0.55,0.7,17);
                if(nbin>1){
                    //       if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[ipt],pt[ipt+1]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.45,0.85,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.78,17);
                    if(ipad==3)drawText(Form("%s jet (anti-k_{T},R=0.3)",level.Data()),0.25,0.85,16);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.78,17);
                }
                else {
                    drawText(Form("%.2f<r<%.2f (GeV/c)", rbin[current],rbin[current+1]),0.6,0.75,17);
                }
                
            }
         //   drawCMS(0.20,0.85,pbpbLumi);

       //     t1->AddEntry(ppsmearRadiusRho[current],"pp Data", "P");
        //    t1->AddEntry(ppRadiusRho,"pp Raw", "P");  
         //   t1->AddEntry(ppsmearbkgRadiusRho[current],"Data Bkg", "P");  
          //  t1->AddEntry(ppbkgRadiusRho,"pp Bkg", "P");

          //  t1->Draw("same");
            //            if(nbin==1){
            //                c1->cd(2);
            //                dummy->SetAxisRange(0.5, 1.5, "Y");
            //                dummy->SetYTitle("Gen_matched/Gen_only");
            //                dummy->DrawCopy();
            //                TGraphErrors * ratio = new TGraphErrors(nfile, rad, rhoratio[0], errR, errrhoratio[0]);
            ////                TH1F * ratio =(TH1F*)RhoGen[0]->Clone("ratio");
            ////                ratio->Divide(RadiusRho[0]);
            //                ratio->Draw("same PE");
            //
            //            }
            if(SavePlot)c1->Print(Form("%s/%s%s%sJetPt%.f_%.fDiffJSwith%sBkg.gif",plotsdir.Data(),type.Data(), level.Data(),sample.Data(), pt[ipt],pt[ipt+1],bkg.Data()));  
            
            if(SaveFile){
                if(IsMC)
                    TFile * outf = new TFile(Form("%s%s%sIncJet%sDiffJSbkg.root", type.Data(), sample.Data(), level.Data(), bkg.Data()), "UPDATE");
                else
                    TFile * outf = new TFile(Form("%s%sIncJet%sDiffJSbkg.root", type.Data(), level.Data(), bkg.Data()), "UPDATE");
                ppRadiusRho->Write(Form("ppRhoJetPt%.f_%.f_Cen%d-%d", pt[ipt],pt[ipt+1], centr[0],centr[nbin]));
                ppbkgRadiusRho->Write(Form("ppbkgRhoJetPt%.f_%.f_Cen%d-%d", pt[ipt],pt[ipt+1], centr[0],centr[nbin]));
                for(int ibin = 0 ; ibin <nbin; ibin++){                    
                    RadiusRho[ibin]->Write(Form("PbPbRhoJetPt%.f_%.f_Cen%d-%d", pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                    bkgRadiusRho[ibin]->Write(Form("PbPbbkgRhoJetPt%.f_%.f_Cen%d-%d", pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                    if(lev==kReco){
                    ppsmearRadiusRho[ibin]->Write(Form("ppsmearRhoJetPt%.f_%.f_Cen%d-%d", pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                    ppsmearbkgRadiusRho[ibin]->Write(Form("ppsmearbkgRhoJetPt%.f_%.f_Cen%d-%d", pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]));
                    }
                }
                outf->Close();  
            }

            c1->Update();   
            break ;
        case kInt:
         //   if(nbin>1)c1->Divide(nbin/2, 2,0, 0);
            //   else c1->Divide(1, 2,0, 0);           
            for(int ipad = 1 ; ipad <=nbin; ipad++){
                c1->cd(ipad);
                //   c1->cd(ipad)->SetLogy();
                dummy->GetYaxis()->SetTitle("#psi (r)");
                //if(current==4)dummy->SetAxisRange(5.e-3, 20., "Y") ;
                // else 
                dummy->SetAxisRange(0., 1.5, "Y") ;
                dummy->DrawCopy();
                if(lev==kReco){
                ppsmearRadiusPsi[nbin-ipad] = new TGraphErrors(nptbin, pt, ppsmearpsi[nbin-ipad][current], 0, errppsmearpsi[nbin-ipad][current]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                ppsmearRadiusPsi[nbin-ipad]->SetMarkerStyle(20);
                ppsmearRadiusPsi[nbin-ipad]->SetMarkerColor(1);
                ppsmearRadiusPsi[nbin-ipad]->SetLineColor(1);
                ppsmearRadiusPsi[nbin-ipad]->Draw("same PE") ;
                
                ppsmearbkgRadiusPsi[nbin-ipad] = new TGraphErrors(nptbin, pt,ppsmearbkgpsi[nbin-ipad], errR, errppsmearbkgpsi[nbin-ipad]);
                ppsmearbkgRadiusPsi[nbin-ipad]->SetMarkerStyle(22);
                ppsmearbkgRadiusPsi[nbin-ipad]->SetMarkerColor(1);
                ppsmearbkgRadiusPsi[nbin-ipad]->SetLineColor(1);
                ppsmearbkgRadiusPsi[nbin-ipad]->Draw("same PE") ;
                }
                else {
                    ppRadiusPsi->Draw("same PE") ;
                    ppbkgRadiusPsi->Draw("same PE") ;  
                }
                RadiusPsi[nbin-ipad] = new TGraphErrors(nptbin, pt, psi[nbin-ipad][current], 0, errpsi[nbin-ipad][current]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                RadiusPsi[nbin-ipad]->SetMarkerStyle(24);
                RadiusPsi[nbin-ipad]->SetMarkerColor(2);
                RadiusPsi[nbin-ipad]->SetLineColor(2);
                RadiusPsi[nbin-ipad]->Draw("same PE") ;
                
                bkgRadiusPsi[nbin-ipad] = new TGraphErrors(nptbin, pt, bkgpsi[nbin-ipad][current], 0, errbkgpsi[nbin-ipad][current]);
                bkgRadiusPsi[nbin-ipad]->SetMarkerStyle(26);
                bkgRadiusPsi[nbin-ipad]->SetMarkerColor(2);
                bkgRadiusPsi[nbin-ipad]->SetLineColor(2);
                bkgRadiusPsi[nbin-ipad]->Draw("same PE") ;
                                
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[ipad]),0.2,0.65, 15);
                if(ipad==1)drawText("Ak PF, R=0.3",0.6,0.75,15);
                if(nbin>1){
                    if(ipad==2)drawText(Form("%.2f<r<%.2f (GeV/c)", rbin[current],rbin[current+1]),0.45,0.75,15);
                    if(ipad==nbin-1){
                        t2->AddEntry(RadiusPsi[0],"PbPb", "P");
                        t2->AddEntry(bkgRadiusPsi[0],"PbPb bkg", "P");  
                        if(lev==kReco){
                        t2->AddEntry(ppsmearRadiusPsi[0],"pp smeared", "P");  
                        t2->AddEntry(ppsmearbkgRadiusPsi[0],"smeared bkg", "P");
                        }
                        else {
                            t2->AddEntry(ppRadiusPsi,"pp ", "P");  
                            t2->AddEntry(ppbkgRadiusPsi,"pp bkg", "P");   
                        }
                        t2->Draw("same");
                    }
                }
                else {
                    drawText(Form("%.2f<r<%.2f (GeV/c)", rbin[current],rbin[current+1]),0.45,0.75,15);
                    t2->AddEntry(RadiusPsi[0],"PbPb", "P");
                    t2->AddEntry(ppsmearRadiusPsi[0],"pp smeared", "P");  
                    t2->AddEntry(bkgRadiusPsi[0],"PbPb bkg", "P");  
                    t2->AddEntry(ppsmearbkgRadiusPsi[0],"smeared bkg", "P");
                    t2->Draw("same");
                    
                }
                
            }
            drawCMS(0.20,0.85,pbpbLumi);
            
   //         t1->AddEntry(ppsmearRadiusPsi[current],"pp Data", "P");
            t1->AddEntry(ppRadiusPsi,"pp Raw", "P");  
   //         t1->AddEntry(ppsmearbkgRadiusPsi[current],"Data Bkg", "P");  
            t1->AddEntry(ppbkgRadiusPsi,"pp Bkg", "P");
            t1->Draw("same");
            if(SavePlot)c1->Print(Form("%s/%sJetPt%.f_%.fIntJSwith%Bkg.eps",plotsdir.Data(),type.Data(), pt[ipt],pt[ipt+1],bkg.Data()));               
            c1->Update();   
            
            break ;
        case kDiffSub:
        //    if(nbin>1)c1->Divide(nbin/2, 2,0, 0);
            //  else c1->Divide(1, 2,0, 0); 
            
//            c1->SetLogy();
//            dummy->DrawCopy();
            for(int ipad = 1 ; ipad <=nbin; ipad++){
                c1->cd(ipad);
                c1->cd(ipad)->SetLogy();
                dummy->GetYaxis()->SetTitle(Form("#rho (%.2f<r<%.2f)", rbin[current],rbin[current+1]));
              if(trackcut>2.)dummy->SetAxisRange(8.e-3, 50., "Y") ;
                else 
                dummy->SetAxisRange(7.e-2, 70., "Y") ;
                dummy->DrawCopy();
          //      ppRhoNoBkg->Draw("same PE") ;

            //    t1->AddEntry(ppsmearRhoNoBkg[nbin-ipad], Form("%d-%d%%",centr[nbin-ipad],centr[ipad]), "P");
                RhoNoBkg[nbin-ipad] = new TGraphErrors(nptbin, pt, rhosub[nbin-ipad][current], 0, errrhosub[nbin-ipad][current]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                RhoNoBkg[nbin-ipad]->SetMarkerStyle(20);
                RhoNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                RhoNoBkg[nbin-ipad]->SetMarkerColor(1);
                RhoNoBkg[nbin-ipad]->SetLineColor(1);
            //    if(!IsMC)drawSys(RhoNoBkg[nbin-ipad], sysprofile[nbin-ipad],deltacone/2.,kGray+2, 1001, 1);  
               RhoNoBkg[nbin-ipad]->Draw("same PE") ;
                cout <<"icentrality= " <<nbin-ipad <<endl ;
                RhoNoBkg[nbin-ipad]->Print();

                if(lev==kReco){
            //        ppsmearRhoNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, pprhosub, errR, errpprhosub);
                ppsmearRhoNoBkg[nbin-ipad] = new TGraphErrors(nptbin, pt,ppsmearrhosub[nbin-ipad][current], 0, errppsmearrhosub[nbin-ipad][current]);
                 //   ppRhoNoBkg = new TGraphErrors(nrbin, rad, pprhosub, errR, errpprhosub);
                ppsmearRhoNoBkg[nbin-ipad]->SetMarkerStyle(24);
                ppsmearRhoNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                ppsmearRhoNoBkg[nbin-ipad]->SetMarkerColor(1);
                ppsmearRhoNoBkg[nbin-ipad]->SetLineColor(1);
                ppsmearRhoNoBkg[nbin-ipad]->Draw("same PE") ;
                    cout << "pp smeared: " <<"icentrality= " <<centr[nbin-ipad] << " - "<<centr[nbin-ipad+1]<<endl ;
                    ppsmearRhoNoBkg[nbin-ipad]->Print();

                }
                else ppRhoNoBkg->Draw("same PE") ;
           //     ppRhoNoBkg->Draw("same PE") ;

                //     t1->AddEntry(RhoNoBkg[nbin-ipad], Form("%d-%d%%",centr[nbin-ipad],centr[ipad]), "P");
            //    if(ipad==1)drawCMS(0.5,0.85,pbpbLumi);
            //    if(ipad==1)drawCMS(0.5,0.85,ppLumi);
           //     if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[nbin-ipad+1]),0.65,0.8, 17);
                if(ipad==4){
                    if(IsMC){
                        t1->AddEntry(RhoNoBkg[nbin-1],"PYTHIA+HYDJET", "PL");
                      //  t1->AddEntry(ppRhoNoBkg,"PYTHIA", "P");
                        if(lev==kReco)t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"PYTHIA Ref", "PL"); 
                        else t1->AddEntry(ppRhoNoBkg,"PYTHIA", "PL"); 
                    }
                    else {
                    t1->AddEntry(RhoNoBkg[nbin-1],"PbPb", "PL");
                //        t1->AddEntry(ppRhoNoBkg,"PYTHIA", "PL");
//                     t1->AddEntry(ppRhoNoBkg,"pp measured", "PL");
         //               if(lev==kReco) t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"2011 pp reference", "PL");
                       if(lev==kReco) t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"2013 pp reference", "PL");
               //         if(lev==kReco)t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"PYTHIA reference", "PL");
//                       if(smf=="MySmear") t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"ak3PF pp:Yaxian", "PL");
//                        else t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"PYTHIA reference", "PL");
                  //      t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"TrkCorr in 4 CentBin", "");
                    }
                    t1->Draw("same");
                }
             //   if(ipad==1)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.7,17);
                if(nbin>1){
                    if(IsMC==kFALSE){
                    if(ipad==1) drawText(Form("CMS PbPb, #sqrt{s_{NN}} = 2.76 TeV, "),0.25,0.9,17);
                    if(ipad==2) drawText(Form(" L_{int} = %.f #mub^{-1}", pbpbLumi),0.01,0.9,17);
                    }
                    if(ipad==1){
                       drawText("anti-k_{T} ( R =0.3 ), PF Jets",0.25,0.76,17);
                        drawText(Form("%.2f<r<%.2f (GeV/c)", rbin[current],rbin[current+1]),0.3,0.68,17);
                   //     drawText(Form("p_{T}^{jet} > %.f GeV/c", pt[ipt]),0.45,0.68,17);
                       drawText(Form("%.1f < |#eta|^{jet} < %.f",etalimit, etacut),0.55,0.58,17);
                    }
                    if(ipad==2){
                        if(lev==kGen)drawText(Form("Generator Level"),0.3,0.76,17); 
                   //     drawText(Form("Tracks in cone (#Delta R < 0.3)"),0.1,0.76,17); 
                        drawText(Form("p_{T}^{track} >%.f GeV/c",trackcut),0.3,0.7,17); 
                    }
//                //    if(ipad==4)drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[ipt],pt[ipt+1]),0.25,0.8,17);
//                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", pt[ipt]),0.25,0.8,17);
//                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.7,17);
//                    if(ipad==3)drawText("Ak PF, R=0.3",0.25,0.85,17);
//            //        if(ipad==3)drawText(Form("%s jet(anti-k_{T},R=0.3)",level.Data()),0.2,0.85,17);
//                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.25,0.7,17);
                }
                else {
                    drawText(Form("%.2f<r<%.2f (GeV/c)", rbin[current],rbin[current+1]),0.6,0.75,17);
                }
            }

            for(int ibin = 0 ; ibin <nbin; ibin++){
            c1->cd(ibin+nbin+1);
                if(IsMC)dummy->SetAxisRange(0.82, 1.22, "Y") ;
                else dummy->SetAxisRange(0.3, 1.7, "Y") ;
            //    else dummy->SetAxisRange(0.8, 1.25, "Y") ;
          //      dummy->GetYaxis()->SetTitle("#rho(r)^{reference}/#rho(r)^{measured}");
        //        dummy->GetYaxis()->SetTitle(Form("#rho(r)/#rho(r)^{50-100%%}"));
         //      dummy->GetYaxis()->SetTitle("ratio");
                if(IsMC)
                    dummy->GetYaxis()->SetTitle("#rho(r)^{PYTHIA+HYDJET}/#rho(r)^{PYTHIA}");
                else dummy->GetYaxis()->SetTitle("#rho(r)^{PbPb}/#rho(r)^{pp}");
          //     else dummy->GetYaxis()->SetTitle("Ratio");
           //     else dummy->GetYaxis()->SetTitle("#rho(r)^{PbPb}/#rho(r)^{PYTHIA}");
                dummy->DrawCopy();
                //                dummy->DrawCopy();
                RatioRhoNoBkg[nbin-ibin-1] = new TGraphErrors(nptbin, pt, rhoratio[nbin-ibin-1][current], 0, errrhoratio[nbin-ibin-1][current]);
                //       RatioRhoNoBkg[ibin] = new TGraphErrors(3, radRebin, rhoratioRebin[ibin], errR, 0);
                RatioRhoNoBkg[nbin-ibin-1]->SetMarkerStyle(20);
                RatioRhoNoBkg[nbin-ibin-1]->SetMarkerColor(1);
                RatioRhoNoBkg[nbin-ibin-1]->SetLineColor(1);
                RatioRhoNoBkg[nbin-ibin-1]->SetMarkerSize(1.5);
                RatioRhoNoBkg[nbin-ibin-1]->Draw("same PE") ;
                if(lev==kReco){
                    smearRatioRhoNoBkg[nbin-ibin-1] = new TGraphErrors(nptbin, pt, smearrhoratio[nbin-ibin-1][current], 0, errsmearrhoratio[nbin-ibin-1][current]);
                    //      smearRatioRhoNoBkg[ibin] = new TGraph(nrbin, rad, smearrhoratio[ibin]);
                    smearRatioRhoNoBkg[nbin-ibin-1]->SetMarkerStyle(25);
                    smearRatioRhoNoBkg[nbin-ibin-1]->SetMarkerColor(2);
                    smearRatioRhoNoBkg[nbin-ibin-1]->SetLineColor(2);
                    smearRatioRhoNoBkg[nbin-ibin-1]->SetMarkerSize(1.5);
             //       smearRatioRhoNoBkg[nbin-ibin-1]->Draw("same PE") ;
                }

                cout <<"ibin = " <<nbin-ibin-1 <<endl ;
                RatioRhoNoBkg[nbin-ibin-1]->Print();
//                if(smf=="MySmear"){
//                    if(ibin==1)drawText("Pawan fitting Para.",0.25,0.85,17); 
//                    if(ibin==2 )drawText("Pawan fitting Para.",0.25,0.85,17);  
//                }
//                else {
//                    if(ibin==1)drawText("Yaxian fitting Para.",0.25,0.85,17); 
//                    if(ibin==2 )drawText("Yaxian fitting Para.",0.25,0.85,17);   
//                }
                regSun(80.,1.,300,1.,1, 1);
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.4,0.9, 17);
                if(nbin>1){
                    if(ibin==3){
                        if(IsMC){
                            if(lev==kReco)t2->AddEntry(RatioRhoNoBkg[nbin-1],"PYTHIA+HYDJET/PYTHIA_ref", "PL");
                            if(lev==kReco)t2->AddEntry(smearRatioRhoNoBkg[nbin-1],"Reference/Measured", "PL");
                        }
                        else {
                            t2->AddEntry(RatioRhoNoBkg[nbin-1],"PbPb/pp_measure", "PL");
                           if(lev==kReco) t2->AddEntry(smearRatioRhoNoBkg[nbin-1],"pp_ref/pp_measure", "PL");
//                            t2->AddEntry(RatioRhoNoBkg[nbin-1],"PbPb/PYTHIA", "PL");
//                            t2->AddEntry(smearRatioRhoNoBkg[nbin-1],"PYTHIA ref/PYTHIA", "PL");
                        }
          //              if(lev==kReco)t2->Draw("same");
                    }
                }
           }
            if(SavePlot)
           //   c1->Print(Form("%sJetPtThres%.fTrkCut%.fDiffJSAfter%sBkgSub.C",type.Data(), pt[ipt],trackcut, bkg.Data()));          
                c1->Print(Form("%s/%s%sRadiusBin%.f_%.f%sTrkCut%.fEtaCut%.f%sDiffJSPtDep%dCentAfter%sBkgSub2013pp.pdf",plotsdir.Data(),type.Data(), level.Data(), rbin[current]*100,rbin[current+1]*100, effTab.Data(),trackcut, etacut*10, smf.Data(),nbin, bkg.Data(), date));
//             c1->Print(Form("%s/%s%sJetPtThres%.f%sTrkCut%.f%sDiffJS%dCentAfter%sBkgSubFullStat.pdf",plotsdir.Data(),type.Data(), level.Data(), pt[ipt],effTab.Data(),trackcut, smf.Data(),nbin, bkg.Data(), date));               
         //   c1->Print(Form("%s/%s%sJetSkimCut100PtBin%.f_%.fTrkCut%.fDiffJSAfter%sBkgSub.gif",plotsdir.Data(),type.Data(),level.Data(), pt[ipt],pt[ipt+1],trackcut, bkg.Data()));               
            c1->Update();   
            if(SaveFile){
                TFile * outf = new TFile(Form("%s%sTrig80IncJet%sTrk%sDiffJSbkgSub.root", type.Data(), level.Data(), bkg.Data(), effTab.Data()), "UPDATE");
                ppRhoNoBkg->Write(Form("ppRhoERbkgsub_Cen%d-%d", centr[0],centr[nbin-1]));
                for(int ibin = 0 ; ibin <nbin; ibin++){                    
                    RhoNoBkg[ibin]->Write(Form("PbPbRhoERbkgsub_Cen%d-%d", centr[ibin],centr[ibin+1]));
                    RatioRhoNoBkg[ibin]->Write(Form("RatioRhoERbkgsub_Cen%d-%d", centr[ibin],centr[ibin+1]));
                    if(lev==kReco)ppsmearRhoNoBkg[ibin]->Write(Form("ppsmearRhoERbkgsub_Cen%d-%d", centr[ibin],centr[ibin+1]));
                }
                outf->Close();  
            }

            break ;
        case kDiffSubRcp:
            for(int ipad = 2 ; ipad <=nbin; ipad++){
                c1->cd(ipad-1);
                c1->cd(ipad-1)->SetLogy();
                dummy->GetYaxis()->SetTitle("#rho (r)");
                if(trackcut>2.)dummy->SetAxisRange(8.e-3, 50., "Y") ;
                else 
                    dummy->SetAxisRange(5.e-2, 50., "Y") ;
                dummy->DrawCopy();
                //       ppRhoNoBkg->Draw("same PE") ;
                
                //    t1->AddEntry(ppsmearRhoNoBkg[nbin-ipad], Form("%d-%d%%",centr[nbin-ipad],centr[ipad]), "P");
                RhoNoBkg[nbin-ipad] = new TGraphErrors(nptbin, pt, rhosub[nbin-ipad][current], 0, errrhosub[nbin-ipad][current]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                RhoNoBkg[nbin-ipad]->SetMarkerStyle(20);
                RhoNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                RhoNoBkg[nbin-ipad]->SetMarkerColor(1);
                RhoNoBkg[nbin-ipad]->SetLineColor(1);
                //    if(!IsMC)drawSys(RhoNoBkg[nbin-ipad], sysprofile[nbin-ipad],deltacone/2.,kGray+2, 1001, 1);  
                RhoNoBkg[nbin-ipad]->Draw("same PE") ;
                if(ipad==2)drawCMS(0.5,0.85,pbpbLumi);
                //    if(ipad==1)drawCMS(0.5,0.85,ppLumi);
                //     if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[nbin-ipad+1]),0.65,0.8, 17);
                if(ipad==3){
                    if(IsMC){
                        t1->AddEntry(RhoNoBkg[nbin-1],"PYTHIA+HYDJET", "PL");
                    }
                    else {
                        t1->AddEntry(RhoNoBkg[nbin-ipad],"PbPb #sqrt{s}=2.76 TeV", "PL");
                        //    t1->AddEntry(ppRhoNoBkg,"pp measured", "PL");  
                        //        t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"PYTHIA reference", "PL");
                    }
                    t1->Draw("same");
                }
                //   if(ipad==1)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.7,17);
                if(nbin>1){
                    //         if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[ipt],pt[ipt+1]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.25,0.8,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.7,17);
                    if(ipad==3)drawText("Ak PF, R=0.3",0.55,0.75,17);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.7,17);
                }
                else {
                    drawText(Form("%.2f<r<%.2f (GeV/c)", rbin[current],rbin[current+1]),0.6,0.75,17);
                }
            }
            for(int ibin = 1 ; ibin <nbin; ibin++){
                c1->cd(ibin+nbin-1);
                if(IsMC)dummy->SetAxisRange(0.5, 1.5, "Y") ;
                else dummy->SetAxisRange(0., 2.1, "Y") ;
                //    else dummy->SetAxisRange(0.8, 1.25, "Y") ;
                //      dummy->GetYaxis()->SetTitle("#rho(r)^{reference}/#rho(r)^{measured}");
                    dummy->GetYaxis()->SetTitle(Form("#rho(r)/#rho(r)^{50-100%%}"));
                //       else dummy->GetYaxis()->SetTitle("#rho(r)^{PbPb}/#rho(r)^{PYTHIA_reference}");
                dummy->DrawCopy();
                //                dummy->DrawCopy();

                RcpRatioRhoNoBkg[nbin-ibin-1] = new TGraphErrors(nptbin, pt, rcprhoratio[nbin-ibin-1][current], 0, errrcprhoratio[nbin-ibin-1][current]);
                //       RatioRhoNoBkg[ibin] = new TGraphErrors(3, radRebin, rhoratioRebin[ibin], errR, 0);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetMarkerStyle(20);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetMarkerColor(1);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetLineColor(1);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetMarkerSize(1.5);
                if(!IsMC)drawSys(RcpRatioRhoNoBkg[nbin-ibin-1], sysRcppercent[nbin-ibin-1],deltacone/2.,7, 1001, 1);  
                RcpRatioRhoNoBkg[nbin-ibin-1]->Draw("same PE") ;
                regSun(0.,1.,0.3,1.,1, 1);
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.25,0.9, 17);
            }            
            
            if(SavePlot)
                c1->Print(Form("%s/%sJetPtThres%.fTrkCut%.fDiffJSRcpAfter%sBkgSubFullStat.pdf",plotsdir.Data(),type.Data(), pt[ipt],trackcut, bkg.Data()));               
             //   c1->Print(Form("%s/%sJetPtThres%.fTrkCut%.fDiffJSRcpAfter%sBkgSub.gif",plotsdir.Data(),type.Data(), pt[ipt],trackcut, bkg.Data()));               
            c1->Update();   
            break ;
        case kRcpRatio:
            for(int ibin = 1 ; ibin <nbin; ibin++){
                c1->cd(ibin);
                if(IsMC)dummy->SetAxisRange(0.5, 1.5, "Y") ;
                else dummy->SetAxisRange(0., 2.1, "Y") ;

                dummy->GetYaxis()->SetTitle(Form("#rho(r)/#rho(r)^{50-100%%}"));
                //       else dummy->GetYaxis()->SetTitle("#rho(r)^{PbPb}/#rho(r)^{PYTHIA_reference}");
                dummy->DrawCopy();
                //                dummy->DrawCopy();
                
                RcpRatioRhoNoBkg[nbin-ibin-1] = new TGraphErrors(nptbin, pt, rcprhoratio[nbin-ibin-1][current], 0, errrcprhoratio[nbin-ibin-1][current]);
                //       RatioRhoNoBkg[ibin] = new TGraphErrors(3, radRebin, rhoratioRebin[ibin], errR, 0);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetMarkerStyle(20);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetMarkerColor(1);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetLineColor(1);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetMarkerSize(1.5);
                RcpRatioRhoNoBkg[nbin-ibin-1]->Draw("same PE") ;
                regSun(0.,1.,0.3,1.,1, 1);
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.25,0.9, 17);
                if(ibin==1)drawCMS(0.5,0.85,pbpbLumi);
                //    if(ipad==1)drawCMS(0.5,0.85,ppLumi);
                //     if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[nbin-ipad+1]),0.65,0.8, 17);
                if(ibin==2){
                    t1->AddEntry(RcpRatioRhoNoBkg[nbin-ibin],"PbPb/PbPb^{(50-100%)}", "PL");
                    t1->Draw("same");
                }
                //   if(ipad==1)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.7,17);
                if(nbin>1){
                    //         if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[ipt],pt[ipt+1]),0.25,0.8,17);
                    if(ibin==3)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.45,0.82,17);
                    if(ibin==3)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.75,17);
                    if(ibin==2)drawText("Ak PF, R=0.3",0.55,0.8,17);
                    if(ibin==2)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.75,17);
                }
                else {
                    drawText(Form("%.2f<r<%.2f (GeV/c)", rbin[current],rbin[current+1]),0.6,0.75,17);
                }

            }            
            
            if(SavePlot)
                //      c1->Print(Form("%s/%sJetPtThres%.f%.sTrkCut%.fIncDiffJS%sRatioSmearAfterBkgSub.eps",plotsdir.Data(),type.Data(), leadingjetcut,bkg.Data(),trackcut,Norm.Data()));               
                c1->Print(Form("%s/%sJetPtThres%.fTrkCut%.fDiffJSRcpRatio.pdf",plotsdir.Data(),type.Data(), pt[ipt],trackcut, bkg.Data()));               
            c1->Update();   
            break ;
        case kDiffRatio:
            for(int ibin = 0 ; ibin <nbin; ibin++){
                c1->cd(ibin+1);
                if(IsMC)dummy->SetAxisRange(0.5, 1.5, "Y") ;
                else dummy->SetAxisRange(0., 2.1, "Y") ;
                //    else dummy->SetAxisRange(0.8, 1.25, "Y") ;
                //      dummy->GetYaxis()->SetTitle("#rho(r)^{reference}/#rho(r)^{measured}");
                //        dummy->GetYaxis()->SetTitle(Form("#rho(r)/#rho(r)^{50-100%%}"));
                //      dummy->GetYaxis()->SetTitle("ratio");
                if(IsMC)
                    dummy->GetYaxis()->SetTitle("#rho(r)^{PYTHIA+HYDJET}/#rho(r)^{PYTHIA}");
                else dummy->GetYaxis()->SetTitle("#rho(r)^{PbPb}/#rho(r)^{pp_reference}");
                //       else dummy->GetYaxis()->SetTitle("#rho(r)^{PbPb}/#rho(r)^{PYTHIA_reference}");
                dummy->DrawCopy();

                RatioRhoNoBkg[nbin-ibin-1] = new TGraphErrors(nptbin, pt, rhoratio[nbin-ibin-1][current], 0, errrhoratio[nbin-ibin-1][current]);
                //       RatioRhoNoBkg[ibin] = new TGraphErrors(3, radRebin, rhoratioRebin[ibin], errR, 0);
                RatioRhoNoBkg[nbin-ibin-1]->SetMarkerStyle(20);
                RatioRhoNoBkg[nbin-ibin-1]->SetMarkerColor(1);
                RatioRhoNoBkg[nbin-ibin-1]->SetLineColor(1);
                RatioRhoNoBkg[nbin-ibin-1]->SetMarkerSize(1.5);
                if(!IsMC)drawSys(RatioRhoNoBkg[nbin-ibin-1], syspercent[nbin-ibin-1],deltacone/2.,7, 1001, 1);  
                //     if(!IsMC)drawSys(RatioRhoNoBkg[nbin-ibin-1], syspercent[nbin-ibin-1],deltacone/2.,kGray+2, 1001, 1);  
                RatioRhoNoBkg[nbin-ibin-1]->Draw("same PE") ;
                regSun(0.,1.,0.3,1.,1, 1);
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.25,0.25, 17);
                if(ibin==0)drawCMS(0.5,0.9,pbpbLumi);
                //    if(ipad==1)drawCMS(0.5,0.85,ppLumi);
                //     if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[nbin-ipad+1]),0.65,0.8, 17);
                if(ibin==0){
                        t1->AddEntry(RatioRhoNoBkg[nbin-1],"PbPb/pp:#sqrt{s}=2.76 TeV", "PL");
                    t1->Draw("same");
                }
                //   if(ipad==1)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.7,17);
                if(nbin>1){
                    //         if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[ipt],pt[ipt+1]),0.25,0.8,17);
                    if(ibin==1)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.25,0.9,17);
                    if(ibin==1)drawText(Form("|#eta|_{jet} < %.f", etacut),0.25,0.83,17);
                    if(ibin==1)drawText("Ak PF, R=0.3",0.25,0.75,17);
                    if(ibin==1)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.25,0.7,17);
                }
                else {
                    drawText(Form("%.2f<r<%.2f (GeV/c)", rbin[current],rbin[current+1]),0.6,0.75,17);
                }

            }            
            
            if(SavePlot)
                    c1->Print(Form("%sJetPtThres%.fTrkCut%.fDiffJSRatio.C",type.Data(), pt[ipt],trackcut, bkg.Data()));               
              //  c1->Print(Form("%s/%sJetPtThres%.fTrkCut%.fDiffJSRatio.gif",plotsdir.Data(),type.Data(), pt[ipt],trackcut, bkg.Data()));               
            c1->Update();   
            break ;
        case kSmearRatio:
            for(int ipad = 1 ; ipad <=nbin; ipad++){
                c1->cd(ipad);
                c1->cd(ipad)->SetLogy();
                dummy->GetYaxis()->SetTitle("#rho (r)");
                if(trackcut>2.)dummy->SetAxisRange(8.e-3, 50., "Y") ;
                else 
                    dummy->SetAxisRange(5.e-2, 50., "Y") ;
                dummy->DrawCopy();
                ppRhoNoBkg->Draw("same PE") ;
                
                //    t1->AddEntry(ppsmearRhoNoBkg[nbin-ipad], Form("%d-%d%%",centr[nbin-ipad],centr[ipad]), "P");
                if(lev==kReco){
                    ppsmearRhoNoBkg[nbin-ipad] = new TGraphErrors(nptbin, pt, ppsmearrhosub[nbin-ipad][current], 0, errppsmearrhosub[nbin-ipad][current]);
                    //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                    ppsmearRhoNoBkg[nbin-ipad]->SetMarkerStyle(24);
                    ppsmearRhoNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                    ppsmearRhoNoBkg[nbin-ipad]->SetMarkerColor(1);
                    ppsmearRhoNoBkg[nbin-ipad]->SetLineColor(1);
                    ppsmearRhoNoBkg[nbin-ipad]->Draw("same PE") ;
                }
                else ppRhoNoBkg->Draw("same PE") ;
                if(ipad==1)drawCMS(0.5,0.85,ppLumi);
                //     if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[nbin-ipad+1]),0.65,0.8, 17);
                if(ipad==2){
                    if(IsMC){
                        t1->AddEntry(ppRhoNoBkg,"PYTHIA", "P");  
                    //    t1->AddEntry(ppsmearRhoNoBkg[nbin-ipad],"PYTHIA reference", "PL"); 
                        t1->AddEntry(ppsmearRhoNoBkg[nbin-ipad],"PYTHIA smeared", "PL"); 
                    }
                    else {
                        t1->AddEntry(ppRhoNoBkg,"pp measured", "PL");  
                     //   t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"pp reference", "PL");
                        if(smf=="MySmear")t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"pp smeared Yaxian", "PL");
                        else t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"pp smeared Pawan", "PL");
                    }
                    t1->Draw("same");
                }
                if(nbin>1){
                    //      if(ipad==4)drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[ipt],pt[ipt+1]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.25,0.8,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.7,17);
                    if(ipad==3)drawText(Form("ak3PF jet,R=0.3",level.Data()),0.2,0.85,17);
//                    if(ipad==3)drawText(Form("%s jet(anti-k_{T},R=0.3)",level.Data()),0.2,0.85,17);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.7,17);
                }
                else {
                    drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[ipt],pt[ipt+1]),0.6,0.75,17);
                }
            }
            for(int ibin = 0 ; ibin <nbin; ibin++){
                c1->cd(ibin+nbin+1);
                dummy->SetAxisRange(0.9, 1.12, "Y") ;
      //          dummy->GetYaxis()->SetTitle("#rho(r)^{reference}/#rho(r)^{measured}");
                dummy->GetYaxis()->SetTitle("#rho(r)^{smearedOnly}/#rho(r)^{measured}");
                dummy->DrawCopy();
                if(lev==kReco){
                    smearRatioRhoNoBkg[nbin-ibin-1] = new TGraphErrors(nptbin, pt, smearrhoratio[nbin-ibin-1][current], 0, errsmearrhoratio[nbin-ibin-1][current]);
                    //      smearRatioRhoNoBkg[ibin] = new TGraph(nrbin, rad, smearrhoratio[ibin]);
                    smearRatioRhoNoBkg[nbin-ibin-1]->SetMarkerStyle(24);
                    smearRatioRhoNoBkg[nbin-ibin-1]->SetMarkerColor(1);
                    smearRatioRhoNoBkg[nbin-ibin-1]->SetLineColor(1);
                    smearRatioRhoNoBkg[nbin-ibin-1]->SetMarkerSize(1.5);
                    smearRatioRhoNoBkg[nbin-ibin-1]->Draw("same PE") ;
                }
//                if(smf=="MySmear"){
//                    if(ibin==1)drawText("Pawan fitting Para.",0.25,0.85,17); 
//                    if(ibin==2 )drawText("Pawan fitting Para.",0.25,0.85,17);  
//                }
//                else {
//                    if(ibin==1)drawText("Yaxian fitting Para.",0.25,0.85,17); 
//                    if(ibin==2 )drawText("Yaxian fitting Para.",0.25,0.85,17);   
//                }

                regSun(0.,1.,0.3,1.,1, 1);
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.25,0.9, 17);
            }            
            
            if(SavePlot)
                //  c1->Print(Form("%sJetPtThres%.fTrkCut%.fDiffJSAfter%sBkgSub.C",type.Data(), pt[ipt],trackcut, bkg.Data()));          
                c1->Print(Form("%s/%s%sJetPtThres%.fTrkCut%.f%sSmearedDiffJSAfter%sBkgSub.pdf",plotsdir.Data(),type.Data(), level.Data(), pt[ipt],trackcut, smf.Data(),bkg.Data()));               
            //  c1->Print(Form("%s/%s%sJetPtBin%.f_%.fTrkCut%.fDiffJSAfter%sBkgSubppRef.gif",plotsdir.Data(),type.Data(),level.Data(), pt[ipt],pt[ipt+1],trackcut, bkg.Data()));               
            c1->Update();   
            break ;
        case kIntSub:
            //    if(nbin>1)c1->Divide(nbin/2, 2,0, 0);
            //  else c1->Divide(1, 2,0, 0); 
            
            //            c1->SetLogy();
            //            dummy->DrawCopy();
            ppPsiNoBkg->Draw("same PE") ;
            
            for(int ipad = 1 ; ipad <=nbin; ipad++){
                c1->cd(ipad);
                c1->cd(ipad)->SetLogy();
                dummy->GetYaxis()->SetTitle("1-#psi (r)");
                if(trackcut>2.)dummy->SetAxisRange(8.e-4, 5., "Y") ;
                else 
                    dummy->SetAxisRange(5.e-3, 5., "Y") ;
                dummy->DrawCopy();
                //      ppRhoNoBkg->Draw("same PE") ;
                
                //    t1->AddEntry(ppsmearRhoNoBkg[nbin-ipad], Form("%d-%d%%",centr[nbin-ipad],centr[ipad]), "P");
                PsiNoBkg[nbin-ipad] = new TGraphErrors(nptbin, pt, psisub[nbin-ipad][current], 0, errpsisub[nbin-ipad][current]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                PsiNoBkg[nbin-ipad]->SetMarkerStyle(20);
                PsiNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                PsiNoBkg[nbin-ipad]->SetMarkerColor(1);
                PsiNoBkg[nbin-ipad]->SetLineColor(1);
             //   if(!IsMC)drawSys(RhoNoBkg[nbin-ipad], sysprofile[nbin-ipad],deltacone/2.,7, 1001, 1);  
                //    if(!IsMC)drawSys(RhoNoBkg[nbin-ipad], sysprofile[nbin-ipad],deltacone/2.,kGray+2, 1001, 1);  
                PsiNoBkg[nbin-ipad]->Draw("same PE") ;
                cout <<"icentrality= " <<nbin-ipad <<endl ;
                PsiNoBkg[nbin-ipad]->Print();
                
                if(lev==kReco){
                    ppsmearPsiNoBkg[nbin-ipad] = new TGraphErrors(nptbin, pt, ppsmearpsisub[nbin-ipad][current], 0, errppsmearpsisub[nbin-ipad][current]);
                    //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                    ppsmearPsiNoBkg[nbin-ipad]->SetMarkerStyle(24);
                    ppsmearPsiNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                    ppsmearPsiNoBkg[nbin-ipad]->SetMarkerColor(1);
                    ppsmearPsiNoBkg[nbin-ipad]->SetLineColor(1);
                    ppsmearPsiNoBkg[nbin-ipad]->Draw("same PE") ;
                    cout << "pp smeared: " <<"icentrality= " <<centr[nbin-ipad] << " - "<<centr[nbin-ipad+1]<<endl ;
                    ppsmearPsiNoBkg[nbin-ipad]->Print();
                    
                }
                else ppPsiNoBkg->Draw("same PE") ;
                //     t1->AddEntry(RhoNoBkg[nbin-ipad], Form("%d-%d%%",centr[nbin-ipad],centr[ipad]), "P");
                if(ipad==1)drawCMS(0.5,0.85,pbpbLumi);
                //    if(ipad==1)drawCMS(0.5,0.85,ppLumi);
                //     if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[nbin-ipad+1]),0.65,0.8, 17);
                if(ipad==2){
                    if(IsMC){
                        t1->AddEntry(PsiNoBkg[nbin-1],"PYTHIA+HYDJET", "PL");
                        //      t1->AddEntry(ppRhoNoBkg,"pp unsmeard", "P");  
                        if(lev==kReco)t1->AddEntry(ppsmearPsiNoBkg[nbin-1],"PYTHIA", "PL"); 
                        else t1->AddEntry(ppPsiNoBkg,"PYTHIA", "PL"); 
                    }
                    else {
                        t1->AddEntry(PsiNoBkg[nbin-1],"PbPb #sqrt{s}=2.76 TeV", "PL");
                        //    t1->AddEntry(ppRhoNoBkg,"pp measured", "PL");  
                        //         t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"pp reference", "PL");
                        if(smf=="MySmear") t1->AddEntry(ppsmearPsiNoBkg[nbin-1],"ak3PF pp:Yaxian", "PL");
                        else t1->AddEntry(ppsmearPsiNoBkg[nbin-1],"ak3PF pp reference", "PL");
                        //      t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"TrkCorr in 4 CentBin", "");
                    }
                    t1->Draw("same");
                }
                //   if(ipad==1)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.7,17);
                if(nbin>1){
                    //     if(ipad==4)drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[ipt],pt[ipt+1]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", pt[ipt]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.7,17);
                    if(ipad==3)drawText("Ak PF, R=0.3",0.25,0.85,17);
                    //        if(ipad==3)drawText(Form("%s jet(anti-k_{T},R=0.3)",level.Data()),0.2,0.85,17);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.25,0.7,17);
                }
                else {
                    drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[ipt],pt[ipt+1]),0.6,0.75,17);
                }
            }
            for(int ibin = 0 ; ibin <nbin; ibin++){
                c1->cd(ibin+nbin+1);
                if(IsMC)dummy->SetAxisRange(0.5, 1.5, "Y") ;
                else dummy->SetAxisRange(0., 2.1, "Y") ;
                //    else dummy->SetAxisRange(0.8, 1.25, "Y") ;
                //      dummy->GetYaxis()->SetTitle("#rho(r)^{reference}/#rho(r)^{measured}");
                //        dummy->GetYaxis()->SetTitle(Form("#rho(r)/#rho(r)^{50-100%%}"));
                //      dummy->GetYaxis()->SetTitle("ratio");
                if(IsMC)
                    dummy->GetYaxis()->SetTitle("(1-#psi(r)^{PYTHIA+HYDJET})/(1-#psi(r)^{PYTHIA})");
                else dummy->GetYaxis()->SetTitle("(1-#psi(r)^{PbPb})/(1-#psi(r)^{pp})");
                //     else dummy->GetYaxis()->SetTitle("#rho(r)^{PbPb}/#rho(r)^{PYTHIA_reference}");
                dummy->DrawCopy();
                RatioPsiNoBkg[nbin-ibin-1] = new TGraphErrors(nptbin, pt, psiratio[nbin-ibin-1][current], 0, errpsiratio[nbin-ibin-1][current]);
                //       RatioRhoNoBkg[ibin] = new TGraphErrors(3, radRebin, rhoratioRebin[ibin], errR, 0);
                RatioPsiNoBkg[nbin-ibin-1]->SetMarkerStyle(20);
                RatioPsiNoBkg[nbin-ibin-1]->SetMarkerColor(1);
                RatioPsiNoBkg[nbin-ibin-1]->SetLineColor(1);
                RatioPsiNoBkg[nbin-ibin-1]->SetMarkerSize(1.5);
             //   if(!IsMC)drawSys(RatioRhoNoBkg[nbin-ibin-1], syspercent[nbin-ibin-1],deltacone/2.,7, 1001, 1);  
                //     if(!IsMC)drawSys(RatioRhoNoBkg[nbin-ibin-1], syspercent[nbin-ibin-1],deltacone/2.,kGray+2, 1001, 1);  
                RatioPsiNoBkg[nbin-ibin-1]->Draw("same PE") ;
                cout <<"ibin = " <<nbin-ibin-1 <<endl ;
                RatioPsiNoBkg[nbin-ibin-1]->Print();
                //                if(smf=="MySmear"){
                //                    if(ibin==1)drawText("Pawan fitting Para.",0.25,0.85,17); 
                //                    if(ibin==2 )drawText("Pawan fitting Para.",0.25,0.85,17);  
                //                }
                //                else {
                //                    if(ibin==1)drawText("Yaxian fitting Para.",0.25,0.85,17); 
                //                    if(ibin==2 )drawText("Yaxian fitting Para.",0.25,0.85,17);   
                //                }
                regSun(0.,1.,0.3,1.,1, 1);
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.25,0.9, 17);
            }            
            
            if(SavePlot)
                //  c1->Print(Form("%sJetPtThres%.fTrkCut%.fDiffJSAfter%sBkgSub.C",type.Data(), pt[ipt],trackcut, bkg.Data()));          
                c1->Print(Form("%s/%s%sJetPtThres%.fEff4CentBinTrkCut%.f%sIntJS%dCentAfter%sBkgSubFullStat.pdf",plotsdir.Data(),type.Data(), level.Data(), pt[ipt],trackcut, smf.Data(),nbin, bkg.Data()));               
            //   c1->Print(Form("%s/%s%sJetPtBin%.f_%.fTrkCut%.fDiffJSAfter%sBkgSub.gif",plotsdir.Data(),type.Data(),level.Data(), pt[ipt],pt[ipt+1],trackcut, bkg.Data()));               
            c1->Update();   
            break ;
        default :
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
    TLatex *cms = new TLatex(px,py,"CMS Preliminary");
    cms->SetTextFont(63);
    cms->SetTextSize(17);
    cms->SetNDC();
    cms->Draw();
    TLatex *lumi = new TLatex(px,py-0.10,Form("#intL dt = %.1f #mub^{-1}",nLumi));
    lumi->SetTextFont(63);
    lumi->SetTextSize(15);
    lumi->SetNDC();
    if(!IsMC)lumi->Draw();
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
    h->SetTitleFont(44);
    h->SetTitleSize(18);
    h->SetLabelSize(20,"Y");
    h->SetLabelSize(20,"X");
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
