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

const bool IsMC=kTRUE ;
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
const int nptbin = 6 ;
//const double pt[]={100., 105., 110., 120.,130., 140., 150.,160.,200., 500.};
//const double pt[]={110., 120., 130., 160., 200., 300., 500.};
const double pt[]={100., 300.};
//const double pt[]={100., 110., 150., 200., 300.};

//const double pt[]={100., 120., 140.,500.};
const double subpt[]={40., 50., 60., 70., 80., 100., 120};
const double radRebin[] = {0.05, 0.15, 0.25};

const int nrbin = 6 ;
double rbin[nrbin+1]={0.0,0.05,0.10, 0.15, 0.20,0.25,0.30};

const double kAj[]={0.0,0.1,0.4,1.0};

const Double_t pbpbLumi = 150 ;
const Double_t ppLumi = 5.3 ;

enum Option_t {kDiff, kInt, kDiffSub, kDiffSubRcp, kRcpRatio, kDiffRatio, kSmearRatio, kCompJSRatioTh, kIntSub} ; 
enum Method_t {kCorr, kRaw} ;
enum Level_t {kReco, kGen} ;
//---------------------------------------------------------------------

void InclusiveJS(Option_t opt=kDiff, Level_t lev =kReco, Method_t met=kCorr, int current = 0){
    for(int i = 0 ; i <nrbin; i++){
        //        if(opt==kInt) rad[i]=rbin[i+1] ;
        //        else         
        rad[i]=deltacone/2.+i*deltacone ;
        ratio[i]=rad[i]/conesize;
        errR[i]=deltacone/2. ;
    }
    TString Norm="NormJet" ;    
    TString trkSel ="Sim" ;
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
//    const int nbin = 6 ;
//    const int centr[] ={0,5,10,30,50,70,90}; 
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
    
    TH1F * ppIntegrated[nrbin] ;
    TH1F * ppbkgIntegrated[nrbin] ;
    
    TH1F * ppdifferential[nrbin] ;
    TH1F * ppbkgdifferential[nrbin] ;
    
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
    
    TH1F * ppsmearIntegrated[nbin][nrbin] ;
    TH1F * ppsmearbkgIntegrated[nbin][nrbin] ;
    
    TH1F * ppsmeardifferential[nbin][nrbin] ;
    TH1F * ppsmearbkgdifferential[nbin][nrbin] ;
    
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
    
    TH1F * Integrated[nbin][nrbin] ;
    TH1F * bkgIntegrated[nbin][nrbin] ;
    
    TH1F * differential[nbin][nrbin] ;
    TH1F * bkgdifferential[nbin][nrbin] ;
    
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
    
    double rho[nbin][nrbin];
    double errrho[nbin][nrbin];
    double psi[nbin][nrbin];
    double errpsi[nbin][nrbin];
    double bkgrho[nbin][nrbin];
    double errbkgrho[nbin][nrbin];
    double bkgpsi[nbin][nrbin];
    double errbkgpsi[nbin][nrbin];
    
    double rhosub[nbin][nrbin];
    double errrhosub[nbin][nrbin];
    double psisub[nbin][nrbin];
    double errpsisub[nbin][nrbin];

    double pprho[nrbin];
    double errpprho[nrbin];
    double pppsi[nrbin];
    double errpppsi[nrbin];
    
    double ppbkgrho[nrbin];
    double errppbkgrho[nrbin];
    double ppbkgpsi[nrbin];
    double errppbkgpsi[nrbin];
    
    double pprhosub[nrbin];
    double errpprhosub[nrbin];
    double pppsisub[nrbin];
    double errpppsisub[nrbin];
    
    double ppsmearrho[nbin][nrbin];
    double errppsmearrho[nbin][nrbin];
    double ppsmearpsi[nbin][nrbin];
    double errppsmearpsi[nbin][nrbin];
    
    double ppsmearbkgrho[nbin][nrbin];
    double errppsmearbkgrho[nbin][nrbin];
    double ppsmearbkgpsi[nbin][nrbin];
    double errppsmearbkgpsi[nbin][nrbin];
    
    double ppsmearrhosub[nbin][nrbin];
    double errppsmearrhosub[nbin][nrbin];
    double ppsmearpsisub[nbin][nrbin];
    double errppsmearpsisub[nbin][nrbin];

    double rhoratio[nbin][nrbin];
    double errrhoratio[nbin][nrbin];
    double psiratio[nbin][nrbin];
    double errpsiratio[nbin][nrbin];

    double rcprhoratio[nbin][nrbin];
    double errrcprhoratio[nbin][nrbin];

    double smearrhoratio[nbin][nrbin];
    double errsmearrhoratio[nbin][nrbin];
    double smearpsiratio[nbin][nrbin];
    double errsmearpsiratio[nbin][nrbin];

    double sumrho = 0.;
    double scalepsi[nbin];
    double scalerho[nbin];
    double scalepppsi = 0.0;
    double scalepprho = 0.0;
    double scaleppsmearpsi[nbin];
    double scaleppsmearrho[nbin];

    double rhoRebin[nbin][3];
    double errrhoRebin[nbin][3];
    double smearrhoRebin[nbin][3];
    double errsmearrhoRebin[nbin][3];
    double rhoratioRebin[nbin][3];

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
    TLegend *t2=new TLegend(0.20,0.65,0.7,0.82);
//    TLegend *t2=new TLegend(0.20,0.25,0.35,0.4);
    t2->SetFillColor(0);
    t2->SetBorderSize(0);
    t2->SetFillStyle(0);
    t2->SetTextFont(63);
    t2->SetTextSize(17);
    // open the pp data file
    if(IsMC){
        if(lev==kGen){
            if(DoOneSample)
                ppfileName =  Form("MCPP_RefJetPt100_%sTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1pthat80_JEC_ppHiIterativeTrack_P01_prod24_v84_merged_forest_0.root", trkSel.Data(), trackcut);
else 
    ppfileName =  Form("mergedCSdiff_MCPP_RefJetPt100_%sTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_ppHiIterativeTrack_P01_prod24_v84_merged_forest_0.root", trkSel.Data(),trackcut); 

        }
        else {
            switch(met){
                case kCorr:
                    if(DoOneSample)
                        ppfileName =  Form("MCPP_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3pthat80_mergedFile.root", effTab.Data(), trackcut);   
                    else 
                //    ppfileName =  Form("mergedCSdiff_MCPP_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_merged.root", effTab.Data(), trackcut);
                    ppfileName =  Form("mergedCSdiff_MCPP2013_Ak3PFIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_ppHiIterativeTrack_merged.root", effTab.Data(), trackcut);              
                    break ;
                case kRaw:
                    if(DoOneSample)
                      ppfileName =  Form("MCPP_Ak3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3pthat80_mergedFile.root", trackcut); 
                    else 
          //              ppfileName =  Form("mergedCSdiff_MCPP_Ak3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_merged.root", trackcut);
                    ppfileName =  Form("mergedCSdiff_MCPP2013_Ak3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_ppHiIterativeTrack_merged.root", trackcut); 
                    break ;
            }
        }
    }
    else {
        switch(met){
            case kCorr:

          //     ppfileName =  Form("DATAPP2011_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_pp_merged_full.root", effTab.Data(), trackcut);
                ppfileName =  Form("DATAPP_Ak3PFIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut);
                break ;
            case kRaw:
                ppfileName =  Form("DATAPP_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut);
//                ppfileName =  Form("DATAPP2011_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_pp_merged_full.root", effTab.Data(), trackcut);
                break ;
        }  
    }
    if(lev==kGen) ppf = TFile::Open(Form("%s/MyselfTrkEff/AfterCWR/%s", kHomeDir, ppfileName.Data()), "readonly");
   // if(lev==kGen) ppf = TFile::Open(Form("/Users/ymao/group/CMS/anaOutputs/JetShape/GenJet/%s", ppfileName.Data()), "readonly");
    else {
        if(met==kRaw) ppf = TFile::Open(Form("%s/%s", kHomeDir, ppfileName.Data()), "readonly");
        else ppf = TFile::Open(Form("%s/MyselfTrkEff/%s", kHomeDir, ppfileName.Data()), "readonly");

    }
 //   ppf = TFile::Open(Form("/Users/ymao/group/CMS/anaOutputs/JetShape/GenJet/%s", ppfileName.Data()), "readonly");
 //   ppf = TFile::Open(Form("/Users/ymao/group/CMS/anaOutputs/JetShape/FullStatPbPb/RebinJetPt/TrkCorr4CentBin/%s", ppfileName.Data()), "readonly");
    cout <<"reading pp file: " << Form("%s/%s", kHomeDir, ppfileName.Data()) <<endl ;    
    for(int ir =0 ; ir <nrbin; ir++){  
        pprho[ir] = 0 ;
        errpprho[ir] = 0 ;
        pppsi[ir]=0 ;
        errpppsi[ir]=0 ;  
        ppbkgrho[ir] = 0 ;
        errppbkgrho[ir] = 0 ;
        ppbkgpsi[ir]=0 ;
        errppbkgpsi[ir]=0 ;  
        pprhosub[ir] = 0 ;
        errpprhosub[ir] = 0 ;
        pppsisub[ir]=0 ;
        errpppsisub[ir]=0 ;  
        ppDiffJS[ir]=(TH2F*)ppf->Get(Form("differentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,0, 100));  
        ppIntJS[ir]=(TH2F*)ppf->Get(Form("IntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,0, 100));  
        ppbkgIntJS[ir]=(TH2F*)ppf->Get(Form("bkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,0, 100));  
        ppbkgDiffJS[ir]=(TH2F*)ppf->Get(Form("bkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,0, 100)); 
        ppIntJS[ir]->Sumw2();
        ppIntegrated[ir]=(TH1F*)ppIntJS[ir]->ProjectionY(Form("ppJetPt%.f_%.fPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,0, 100),ppIntJS[ir]->GetXaxis()->FindBin(pt[current]), ppIntJS[ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
        
        ppbkgIntJS[ir]->Sumw2();
        ppbkgIntegrated[ir]=(TH1F*)ppbkgIntJS[ir]->ProjectionY(Form("ppJetPt%.f_%.fbkgPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,0, 100),ppbkgIntJS[ir]->GetXaxis()->FindBin(pt[current]), ppbkgIntJS[ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
        ppDiffJS[ir]->Sumw2();
        ppdifferential[ir]=(TH1F*)ppDiffJS[ir]->ProjectionY(Form("ppJetPt%.f_%.fRhodR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,0, 100),ppDiffJS[ir]->GetXaxis()->FindBin(pt[current]), ppDiffJS[ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
        cout << pt[current]<<" < Jet pt < " << pt[current+1] <<"Njets =" << ppdifferential[ir]->GetEntries()<<endl;
        ppbkgDiffJS[ir]->Sumw2();
        ppbkgdifferential[ir]=(TH1F*)ppbkgDiffJS[ir]->ProjectionY(Form("ppJetPt%.f_%.fbkgRhodR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,0, 100),ppbkgDiffJS[ir]->GetXaxis()->FindBin(pt[current]), ppbkgDiffJS[ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
        
        pppsi[ir]=ppIntegrated[ir]->GetMean(1);
        errpppsi[ir]=ppIntegrated[ir]->GetMeanError(1);
        pprho[ir]=ppdifferential[ir]->GetMean(1);
        errpprho[ir]=ppdifferential[ir]->GetMeanError(1);
        ppbkgpsi[ir]=ppbkgIntegrated[ir]->GetMean(1);
        errppbkgpsi[ir]=ppbkgIntegrated[ir]->GetMeanError(1);
        ppbkgrho[ir]=ppbkgdifferential[ir]->GetMean(1);
        errppbkgrho[ir]=ppbkgdifferential[ir]->GetMeanError(1);
        //nromalize to delta cone size to consistent for definition
        pprho[ir]/=deltacone ;
        ppbkgrho[ir]/=deltacone ;
        errpprho[ir]/=deltacone ;
        errppbkgrho[ir]/=deltacone ;
        //calculate the bkg subtracted JS
        pprhosub[ir]=pprho[ir]-ppbkgrho[ir];
        errpprhosub[ir]=TMath::Sqrt(TMath::Power(errpprho[ir],2)+TMath::Power(errppbkgrho[ir],2));
//        pprhosub[ir]=pprho[ir];
//        errpprhosub[ir]=TMath::Sqrt(TMath::Power(errpprho[ir],2));
//        pppsisub[ir]=pppsi[ir]-ppbkgpsi[ir];
//        errpppsisub[ir]=TMath::Sqrt(TMath::Power(errpppsi[ir],2)+TMath::Power(errppbkgpsi[ir],2));

        scalepprho+=pprhosub[ir]*deltacone;
        //    scaleppsmearrho+=ppsmearrhosub[ir]*deltacone;
    } 
        // scalepprho-=1. ;
    //    scalepppsi=pppsisub[nrbin-1]; 
 //       scaleppsmearpsi=ppsmearpsisub[nrbin-1]; 
        
        //rescale bkg subtracted JS in order to get unity for psi
        for(int ir =0 ; ir <nrbin; ir++){ 
        //    pppsisub[ir]/=scalepppsi ;           
            pprhosub[ir]/=scalepprho ;
        //    errpppsisub[ir]/=scalepppsi ;           
            errpprhosub[ir]/=scalepprho ;
            sumrho+=pprhosub[ir];
            
            for(int i =0 ; i <ir; i++){
                pppsisub[ir]+=pprhosub[i]*deltacone ;   
                errpppsisub[ir]+=TMath::Power(errpprhosub[i]*deltacone,2);
            }
            pppsisub[ir]= 1- pppsisub[ir] ;
            errpppsisub[ir]=TMath::Sqrt(errpppsisub[ir]);
            
            cout <<"ir = "<<ir<<"psi =" << pppsisub[ir] <<"rho =" << pprhosub[ir]<<endl ;           
        }
        cout <<" pp rho =" <<sumrho<<endl ;
        
        ppRadiusRho = new TGraphErrors(nrbin, rad, pprho, errR, errpprho);
        ppRadiusRho->SetMarkerStyle(24);
        ppRadiusRho->SetMarkerColor(2);
        ppRadiusRho->SetMarkerSize(1.5);
        ppRadiusRho->SetLineColor(2);
        ppbkgRadiusRho = new TGraphErrors(nrbin, rad, ppbkgrho, errR, errppbkgrho);
        ppbkgRadiusRho->SetMarkerStyle(30);
        ppbkgRadiusRho->SetMarkerSize(1.5);
        ppbkgRadiusRho->SetMarkerColor(2);
        ppbkgRadiusRho->SetLineColor(2);
        
        ppRadiusPsi = new TGraphErrors(nrbin, rad, pppsi, errR, errpppsi);
        ppRadiusPsi->SetMarkerStyle(21);
        ppRadiusPsi->SetMarkerColor(4);
        ppRadiusPsi->SetLineColor(4);
        ppbkgRadiusPsi = new TGraphErrors(nrbin, rad, ppbkgpsi, errR, errppbkgpsi);
        ppbkgRadiusPsi->SetMarkerStyle(25);
        ppbkgRadiusPsi->SetMarkerColor(4);
        ppbkgRadiusPsi->SetLineColor(4);
        
        ppRhoNoBkg = new TGraphErrors(nrbin, rad, pprhosub, errR, errpprhosub);
        ppRhoNoBkg->SetMarkerStyle(25);
        ppRhoNoBkg->SetMarkerColor(2);
         ppRhoNoBkg->SetMarkerSize(1.5);
        ppRhoNoBkg->SetLineColor(2);
        ppPsiNoBkg = new TGraphErrors(nrbin, rad, pppsisub, errR, errpppsisub);
        ppPsiNoBkg->SetMarkerStyle(25);
        ppPsiNoBkg->SetMarkerColor(1);    
        ppPsiNoBkg->SetMarkerSize(1.5);
        ppPsiNoBkg->SetLineColor(1);        

    // open the HI data file
    if(IsMC){
        if(lev==kGen){
            if(DoOneSample)
                fileName =  Form("MCHI_RefJetPt100_%sTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1pthat80_Dijet80_HydjetDrum_v27_mergedV1.root", trkSel.Data(),trackcut); 
else 
    fileName =  Form("mergedCSdiff_MCHI_RefJetPt100_%sTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1_Dijet_HydjetDrum_v27_mergedV1.root", trkSel.Data(),trackcut); 
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
                    smearppfileName=Form("mergedCSdiff_MCPP2013_Ak3PFwtIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_trkbin3_pp276Dijet_ppHiIterativeTrack_merged.root", effTab.Data(), trackcut, nbin);
                    }
                    break ;
                case kRaw:
                    if(DoOneSample){
                        smearppfileName = Form("MCPP_Ak3PFwtIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3pthat80_mergedFile.root", trackcut);
                        fileName=Form("MCHI_AkPu3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin3_Dijet80_HydjetDrum_v27_mergedV1.root", trackcut, nbin); 
                    }
                    else {
                    fileName=Form("mergedCSdiff_MCHI_AkPu3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_Dijet_HydjetDrum_v27v28.root", effTab.Data(), trackcut);
                    smearppfileName=Form("mergedCSdiff_MCPP2013_Ak3PFwtIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_trkbin3_pp276Dijet_ppHiIterativeTrack_merged.root", effTab.Data(), trackcut, nbin);
                    }
                    break ;
            }
        }
    }
    else {
        switch(met){
            case kCorr:
                fileName=Form("DATAHI_AkPu3PFIncJetPt100_%s%.fEtaCut%.fLimit%.f_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin3_promptskim-hihighpt-hltjet80-pt90-v20.root", effTab.Data(), trackcut, etacut*10, etalimit*10, nbin);
           //     smearppfileName=Form("DATAPP_Ak3PFIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut);

                smearppfileName=Form("DATAPP_Ak3PFwtIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut);
         //       smearppfileName=Form("DATAPP2011_Ak3PFwtIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_pp_merged_full.root", effTab.Data(), trackcut);

                break ;
            case kRaw:
                fileName=Form("DATAHI_AkPu3PFIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_promptskim-hihighpt-hltjet80-pt90-v20.root", trackcut);
                smearppfileName=Form("DATAPP_Ak3PFwtIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut);
//                smearppfileName=Form("DATAPP_Ak3PFwtIncJetPt100_Trk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_pp_merged_full.root", trackcut);
                break ;
        }  
    }
    if(IsMC){
        if(lev==kGen) f = TFile::Open(Form("%s/MyselfTrkEff/AfterCWR/%s", kHomeDir, fileName.Data()), "readonly");
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
            scalepsi[ibin]=0. ;
            scalerho[ibin]=0. ;
            scaleppsmearpsi[ibin]=0. ;
            scaleppsmearrho[ibin]=0. ;
            sumrho = 0.0 ;
            for(int ir =0 ; ir <nrbin; ir++){
                if(lev==kReco){
                ppsmearrho[ibin][ir] = 0 ;
                errppsmearrho[ibin][ir] = 0 ;
                ppsmearpsi[ibin][ir]=0 ;
                errppsmearpsi[ibin][ir]=0 ;  
                ppsmearbkgrho[ibin][ir] = 0 ;
                errppsmearbkgrho[ibin][ir] = 0 ;
                ppsmearbkgpsi[ibin][ir]=0 ;
                errppsmearbkgpsi[ibin][ir]=0 ;  
                ppsmearrhosub[ibin][ir] = 0 ;
                errppsmearrhosub[ibin][ir] = 0 ;
                ppsmearpsisub[ibin][ir]=0 ;
                errppsmearpsisub[ibin][ir]=0 ; 
                ppsmearDiffJS[ibin][ir]=(TH2F*)smearppf->Get(Form("differentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
                ppsmearIntJS[ibin][ir]=(TH2F*)smearppf->Get(Form("IntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
                ppsmearbkgIntJS[ibin][ir]=(TH2F*)smearppf->Get(Form("bkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
                ppsmearbkgDiffJS[ibin][ir]=(TH2F*)smearppf->Get(Form("bkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));          
                ppsmearIntJS[ibin][ir]->Sumw2();
                ppsmearIntegrated[ibin][ir]=(TH1F*)ppsmearIntJS[ibin][ir]->ProjectionY(Form("ppsmearJetPt%.f_%.fPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),ppsmearIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), ppsmearIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
                
                ppsmearbkgIntJS[ibin][ir]->Sumw2();
                ppsmearbkgIntegrated[ibin][ir]=(TH1F*)ppsmearbkgIntJS[ibin][ir]->ProjectionY(Form("ppsmearJetPt%.f_%.fbkgPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),ppsmearbkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), ppsmearbkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
                
                ppsmearDiffJS[ibin][ir]->Sumw2();
                ppsmeardifferential[ibin][ir]=(TH1F*)ppsmearDiffJS[ibin][ir]->ProjectionY(Form("ppsmearJetPt%.f_%.fRhodR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),ppsmearDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), ppsmearDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
                
                ppsmearbkgDiffJS[ibin][ir]->Sumw2();
                ppsmearbkgdifferential[ibin][ir]=(TH1F*)ppsmearbkgDiffJS[ibin][ir]->ProjectionY(Form("ppsmearJetPt%.f_%.fbkgRhodR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),ppsmearbkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), ppsmearbkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
                    
                    cout <<" Centality =" << centr[ibin] << " - " << centr[ibin+1] <<"%% "<< pt[current]<<" < Jet pt < " << pt[current+1] <<"  ppNjets =" << ppsmeardifferential[ibin][ir]->GetEntries()<<endl;

                ppsmearpsi[ibin][ir]=ppsmearIntegrated[ibin][ir]->GetMean(1);
                errppsmearpsi[ibin][ir]=ppsmearIntegrated[ibin][ir]->GetMeanError(1);
                ppsmearrho[ibin][ir]=ppsmeardifferential[ibin][ir]->GetMean(1);
                errppsmearrho[ibin][ir]=ppsmeardifferential[ibin][ir]->GetMeanError(1);
                ppsmearbkgpsi[ibin][ir]=ppsmearbkgIntegrated[ibin][ir]->GetMean(1);
                errppsmearbkgpsi[ibin][ir]=ppsmearbkgIntegrated[ibin][ir]->GetMeanError(1);
                ppsmearbkgrho[ibin][ir]=ppsmearbkgdifferential[ibin][ir]->GetMean(1);
                errppsmearbkgrho[ibin][ir]=ppsmearbkgdifferential[ibin][ir]->GetMeanError(1);
                ppsmearrho[ibin][ir]/=deltacone ;
                ppsmearbkgrho[ibin][ir]/=deltacone ;
                errppsmearrho[ibin][ir]/=deltacone ;
                errppsmearbkgrho[ibin][ir]/=deltacone ;

//                ppsmearrhosub[ibin][ir]=ppsmearrho[ibin][ir];
//                errppsmearrhosub[ibin][ir]=TMath::Sqrt(TMath::Power(errppsmearrho[ibin][ir],2));
                ppsmearrhosub[ibin][ir]=ppsmearrho[ibin][ir]-ppsmearbkgrho[ibin][ir];
                errppsmearrhosub[ibin][ir]=TMath::Sqrt(TMath::Power(errppsmearrho[ibin][ir],2)+TMath::Power(errppsmearbkgrho[ibin][ir],2));
//                ppsmearpsisub[ibin][ir]=ppsmearpsi[ibin][ir]-ppsmearbkgpsi[ibin][ir];
//                errppsmearpsisub[ibin][ir]=TMath::Sqrt(TMath::Power(errppsmearpsi[ibin][ir],2)+TMath::Power(errppsmearbkgpsi[ibin][ir],2));
                    
                    smearrhoratio[ibin][ir] = 0 ;
                    errsmearrhoratio[ibin][ir] = 0 ;
                    smearpsiratio[ibin][ir]=0 ;
                    errsmearpsiratio[ibin][ir]=0 ; 
                }
                rho[ibin][ir] = 0 ;
                errrho[ibin][ir] = 0 ;
                psi[ibin][ir]=0 ;
                errpsi[ibin][ir]=0 ;  
                bkgrho[ibin][ir] = 0 ;
                errbkgrho[ibin][ir] = 0 ;
                bkgpsi[ibin][ir]=0 ;
                errbkgpsi[ibin][ir]=0 ;  

                rhosub[ibin][ir] = 0 ;
                errrhosub[ibin][ir] = 0 ;
                psisub[ibin][ir]=0 ;
                errpsisub[ibin][ir]=0 ;  

                rhoratio[ibin][ir] = 0 ;
                errrhoratio[ibin][ir] = 0 ;
                psiratio[ibin][ir]=0 ;
                errpsiratio[ibin][ir]=0 ;  

                rcprhoratio[ibin][ir] = 0 ;
                errrcprhoratio[ibin][ir] = 0 ;
                IntJS[ibin][ir]=(TH2F*)f->Get(Form("IntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
                bkgIntJS[ibin][ir]=(TH2F*)f->Get(Form("bkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
                
                DiffJS[ibin][ir]=(TH2F*)f->Get(Form("differentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
                bkgDiffJS[ibin][ir]=(TH2F*)f->Get(Form("bkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
                         
                IntJS[ibin][ir]->Sumw2();
                Integrated[ibin][ir]=(TH1F*)IntJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),IntJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), IntJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
                
                bkgIntJS[ibin][ir]->Sumw2();
                bkgIntegrated[ibin][ir]=(TH1F*)bkgIntJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fbkgPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),bkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), bkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
                DiffJS[ibin][ir]->Sumw2();
                differential[ibin][ir]=(TH1F*)DiffJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fRhodR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),DiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), DiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
                bkgDiffJS[ibin][ir]->Sumw2();                
                bkgdifferential[ibin][ir]=(TH1F*)bkgDiffJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fbkgRhodR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),bkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), bkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
                

                cout <<" Centality =" << centr[ibin] << " - " << centr[ibin+1] <<"%% "<< pt[current]<<" < Jet pt < " << pt[current+1] <<"Njets =" << differential[ibin][ir]->GetEntries()<<endl;

                psi[ibin][ir]= Integrated[ibin][ir]->GetMean(1);
                errpsi[ibin][ir]= Integrated[ibin][ir]->GetMeanError(1);
                bkgpsi[ibin][ir]= bkgIntegrated[ibin][ir]->GetMean(1);
                errbkgpsi[ibin][ir]= bkgIntegrated[ibin][ir]->GetMeanError(1);
                rho[ibin][ir]= differential[ibin][ir]->GetMean(1);
                errrho[ibin][ir]= differential[ibin][ir]->GetMeanError(1);
                bkgrho[ibin][ir]= bkgdifferential[ibin][ir]->GetMean(1);
                errbkgrho[ibin][ir]= bkgdifferential[ibin][ir]->GetMeanError(1);

                //nromalize to delta cone size to consistent for definition
                rho[ibin][ir]/=deltacone ;
                bkgrho[ibin][ir]/=deltacone ;
                errrho[ibin][ir]/=deltacone ;
                errbkgrho[ibin][ir]/=deltacone ;

                //calculate the bkg subtracted JS
                rhosub[ibin][ir]=rho[ibin][ir]-bkgrho[ibin][ir];
                //   errrhosub[ibin][ir]=errrho[ibin][ir]+errbkgrho[ibin][ir];
                errrhosub[ibin][ir]=TMath::Sqrt(TMath::Power(errrho[ibin][ir],2)+TMath::Power(errbkgrho[ibin][ir],2));
                
//                rhosub[ibin][ir]=rho[ibin][ir];
//             //   errrhosub[ibin][ir]=errrho[ibin][ir]+errbkgrho[ibin][ir];
//                  errrhosub[ibin][ir]=TMath::Sqrt(TMath::Power(errrho[ibin][ir],2));
             //      errrhosub[ibin][ir]*=rhosub[ibin][ir];
                
//                psisub[ibin][ir]=psi[ibin][ir]-bkgpsi[ibin][ir];
//                  errpsisub[ibin][ir]=TMath::Sqrt(TMath::Power(errpsi[ibin][ir],2)+TMath::Power(errbkgpsi[ibin][ir],2));
           //     if(ir>0){
                if(lev==kReco)scaleppsmearrho[ibin]+=ppsmearrhosub[ibin][ir]*deltacone;
                scalerho[ibin]+=rhosub[ibin][ir]*deltacone;                    
                } 
          //  } 
//            if(lev==kReco) scaleppsmearpsi[ibin]=ppsmearpsisub[ibin][nrbin-1]; 
//            
//            scalepsi[ibin]=psisub[ibin][nrbin-1];
       //     cout <<"ibin = "<<ibin<<"scale =" << scalerho[ibin] <<"rho =" << scalepsi[ibin]<<endl ;
            
            for(int ir =0 ; ir <nrbin; ir++){ 
           //     if(ir>0){
                if(lev==kReco){
          //      ppsmearpsisub[ibin][ir]/=scaleppsmearpsi[ibin] ;           
                ppsmearrhosub[ibin][ir]/=scaleppsmearrho[ibin] ;
           //     errppsmearpsisub[ibin][ir]/=scaleppsmearpsi[ibin] ;           
                errppsmearrhosub[ibin][ir]/=scaleppsmearrho[ibin] ;
                }
            //    psisub[ibin][ir]/=scalepsi[ibin] ;           
                rhosub[ibin][ir]/=scalerho[ibin] ;
            //    errpsisub[ibin][ir]/=scalepsi[ibin] ;           
                errrhosub[ibin][ir]/=scalerho[ibin] ;
        //    }

                sumrho+=rhosub[ibin][ir];
                
                // try to calculate 1-psi instead of psi from rho
                if(lev==kReco){
                    for(int i =0 ; i <ir; i++){
                      ppsmearpsisub[ibin][ir]+=ppsmearrhosub[ibin][i]*deltacone ;   
                      errppsmearpsisub[ibin][ir]+=TMath::Power(errppsmearrhosub[ibin][i]*deltacone,2);
                    }
                    ppsmearpsisub[ibin][ir]= 1- ppsmearpsisub[ibin][ir] ;
                    errppsmearpsisub[ibin][ir]=TMath::Sqrt(errppsmearpsisub[ibin][ir]);
                }
                for(int i =0 ; i <ir; i++){
                    psisub[ibin][ir]+=rhosub[ibin][i]*deltacone ;   
                    errpsisub[ibin][ir]+=TMath::Power(errrhosub[ibin][i]*deltacone,2);
                }
                psisub[ibin][ir]= 1- psisub[ibin][ir] ;
                errpsisub[ibin][ir]=TMath::Sqrt(errpsisub[ibin][ir]);
                
     //           cout <<"ir = "<<ir<<"psi =" << psisub[ibin][ir] <<"rho =" << rhosub[ibin][ir]<<endl ;  
                
                if(IsMC==kTRUE){
                    if(lev==kGen){
                        rhoratio[ibin][ir]=rhosub[ibin][ir]/pprhosub[ir];
                        errrhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errrhosub[ibin][ir]/rhosub[ibin][ir],2)+TMath::Power(errpprhosub[ir]/pprhosub[ir],2));
                        errrhoratio[ibin][ir]*=rhoratio[ibin][ir];
                        if(ir==1) cout <<"pprho ="<<pprhosub[ir] <<"  HI rho ="<<rhosub[ibin][ir]<<endl;
                        psiratio[ibin][ir]=psisub[ibin][ir]/pppsisub[ir];
                        errpsiratio[ibin][ir]=TMath::Sqrt(TMath::Power(errpsisub[ibin][ir]/psisub[ibin][ir],2)+TMath::Power(errpppsisub[ir]/pppsisub[ir],2));
                        errpsiratio[ibin][ir]*=psiratio[ibin][ir];
                    }
                    else {
                        rhoratio[ibin][ir]=rhosub[ibin][ir]/ppsmearrhosub[ibin][ir];
                        errrhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errrhosub[ibin][ir]/rhosub[ibin][ir],2)+TMath::Power(errppsmearrhosub[ibin][ir]/ppsmearrhosub[ibin][ir],2));
                        errrhoratio[ibin][ir]*=rhoratio[ibin][ir];
                        
                        psiratio[ibin][ir]=psisub[ibin][ir]/ppsmearpsisub[ibin][ir];
                        errpsiratio[ibin][ir]=TMath::Sqrt(TMath::Power(errpsisub[ibin][ir]/psisub[ibin][ir],2)+TMath::Power(errppsmearpsisub[ibin][ir]/ppsmearpsisub[ibin][ir],2));
                        errpsiratio[ibin][ir]*=psiratio[ibin][ir];
//                        rhoratio[ibin][ir]=rhosub[ibin][ir]/pprhosub[ir];
//                        errrhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errrhosub[ibin][ir]/rhosub[ibin][ir],2)+TMath::Power(errpprhosub[ir]/pprhosub[ir],2));
//                        errrhoratio[ibin][ir]*=rhoratio[ibin][ir];
//                        if(ir==1) cout <<"pprho ="<<pprhosub[ir] <<"  HI rho ="<<rhosub[ibin][ir]<<endl;
//                        psiratio[ibin][ir]=psisub[ibin][ir]/pppsisub[ir];
//                        errpsiratio[ibin][ir]=TMath::Sqrt(TMath::Power(errpsisub[ibin][ir]/psisub[ibin][ir],2)+TMath::Power(errpppsisub[ir]/pppsisub[ir],2));
//                        errpsiratio[ibin][ir]*=psiratio[ibin][ir];
                        
                        smearrhoratio[ibin][ir]=ppsmearrhosub[ibin][ir]/pprhosub[ir];
                        //   errsmearrhoratio[ibin][ir]=errppsmearrhosub[ibin][ir];
                        errsmearrhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errppsmearrhosub[ibin][ir]/ppsmearrhosub[ibin][ir],2));
                        errsmearrhoratio[ibin][ir]*=smearrhoratio[ibin][ir];
                        //             errsmearrhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errppsmearrhosub[ibin][ir],2)+TMath::Power(errpprhosub[ir],2));
                        smearpsiratio[ibin][ir]=ppsmearpsisub[ibin][ir]/pppsisub[ir];
                        //   errsmearpsiratio[ibin][ir]=errppsmearpsisub[ibin][ir];
                        errsmearpsiratio[ibin][ir]=TMath::Sqrt(TMath::Power(errppsmearpsisub[ibin][ir]/ppsmearpsisub[ibin][ir],2));
                        errsmearpsiratio[ibin][ir]*=smearpsiratio[ibin][ir];
                    }
                }
                else {
//                    rhoratio[ibin][ir]=rhosub[ibin][ir]/pprhosub[ir];
//                    errrhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errrhosub[ibin][ir]/rhosub[ibin][ir],2)+TMath::Power(errpprhosub[ir]/pprhosub[ir],2));
//                    errrhoratio[ibin][ir]*=rhoratio[ibin][ir];
//                    if(ir==1) cout <<"pprho ="<<pprhosub[ir] <<"  HI rho ="<<rhosub[ibin][ir]<<endl;
//                    psiratio[ibin][ir]=psisub[ibin][ir]/pppsisub[ir];
//                    errpsiratio[ibin][ir]=TMath::Sqrt(TMath::Power(errpsisub[ibin][ir]/psisub[ibin][ir],2)+TMath::Power(errpppsisub[ir]/pppsisub[ir],2));
//                    errpsiratio[ibin][ir]*=psiratio[ibin][ir];

                    rhoratio[ibin][ir]=rhosub[ibin][ir]/ppsmearrhosub[ibin][ir];
                    errrhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errrhosub[ibin][ir]/rhosub[ibin][ir],2)+TMath::Power(errppsmearrhosub[ibin][ir]/ppsmearrhosub[ibin][ir],2));
                    errrhoratio[ibin][ir]*=rhoratio[ibin][ir];
                    
                    psiratio[ibin][ir]=psisub[ibin][ir]/ppsmearpsisub[ibin][ir];
                    errpsiratio[ibin][ir]=TMath::Sqrt(TMath::Power(errpsisub[ibin][ir]/psisub[ibin][ir],2)+TMath::Power(errppsmearpsisub[ibin][ir]/ppsmearpsisub[ibin][ir],2));
                    errpsiratio[ibin][ir]*=psiratio[ibin][ir];
                    
                    smearrhoratio[ibin][ir]=ppsmearrhosub[ibin][ir]/pprhosub[ir];
                    //   errsmearrhoratio[ibin][ir]=errppsmearrhosub[ibin][ir];
                    errsmearrhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errppsmearrhosub[ibin][ir]/ppsmearrhosub[ibin][ir],2));
                    errsmearrhoratio[ibin][ir]*=smearrhoratio[ibin][ir];
                    //             errsmearrhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errppsmearrhosub[ibin][ir],2)+TMath::Power(errpprhosub[ir],2));
                    smearpsiratio[ibin][ir]=ppsmearpsisub[ibin][ir]/pppsisub[ir];
                    //   errsmearpsiratio[ibin][ir]=errppsmearpsisub[ibin][ir];
                    errsmearpsiratio[ibin][ir]=TMath::Sqrt(TMath::Power(errppsmearpsisub[ibin][ir]/ppsmearpsisub[ibin][ir],2));
                    errsmearpsiratio[ibin][ir]*=smearpsiratio[ibin][ir];  
                    
                }

             //   cout <<"ir = "<<ir<<"rho ratio =" << rhoratio[ibin][ir] <<"psi ratio =" << psiratio[ibin][ir]<<endl ;  
            }
            cout <<"sum rho =" <<sumrho<<endl ;
    }
    for(int ibin = 0 ; ibin <nbin; ibin++){  
        
        for(int ir =0 ; ir <nrbin; ir++){
            rcprhoratio[ibin][ir]=rhosub[ibin][ir]/rhosub[nbin-1][ir];
            errrcprhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errrhosub[ibin][ir]/rhosub[ibin][ir],2));
            errrcprhoratio[ibin][ir]*=rcprhoratio[ibin][ir];
            
          //  cout <<"ibin = "<<ibin <<"ir = "<<ir <<"rhoratio =" << rhoratio[ibin][ir]<<endl ;  
        }
    }
    TH1F * dummy = new TH1F("dummy", "dummy", 100, -0.05, 1.);
    dummy->SetTitle("") ;
    dummy->SetStats(kFALSE) ;
       if(opt==kCompJSRatioTh) dummy->SetAxisRange(0., 1., "X") ;
       else   dummy->SetAxisRange(-0.01, 0.32, "X") ;

    if(opt==kCompJSRatioTh) dummy->GetXaxis()->SetTitle("r/R");
    else  dummy->GetXaxis()->SetTitle(" r ");

    if(opt==kDiff|| opt==kDiffRatio || opt==kRcpRatio || opt==kCompJSRatioTh)fixedFontHist(dummy,1.1, 1.3);
    else fixedFontHist(dummy,1.8, 2.23);

    //for PbPb profile systematics uncertainty
    //QM period uncertainty (priminary)
//    double sysprofile[nbin][nrbin]={{1.70552, 3.93072, 6.16875, 8.40937, 10.651, 12.893},
//        {1.60929, 3.54131, 5.50656, 7.47889, 9.45386, 11.4301},
//        {1.07522, 2.48205, 3.93005, 5.38611, 6.8451, 8.30547},
//        {1.21473, 2.79387, 4.40546, 6.02352, 7.66711, 9.26547}};
//2012 end of year
//    double sysprofile[nbin][nrbin]={{2.61643, 3.38668, 5.88245, 7.30573, 9.15482, 10.9573},
//        {2.51537, 4.00112, 5.78009, 7.08917, 8.66978, 11.3077},
//        {1.80397, 3.40805, 3.68262, 4.40762, 5.77387, 6.03651},
//        {1.98379, 3.65005, 4.25504, 4.86365, 5.15183, 5.64737}};
    
//    double sysprofile[nbin][nrbin]={{6.32645, 6.39154, 7.02905, 8.57264, 9.27228, 12.4112},
//        {6.98647, 6.9366, 6.80919, 7.77162, 10.2519, 12.4892},
//        {6.53023, 6.73427, 6.98098, 6.37705, 6.69387, 6.72948},
//        {7.1653, 7.05692, 7.44621, 7.35506, 7.78266, 7.71648}};
    double sysprofile[nbin][nrbin]={{7.27875, 8.01062, 9.68922, 10.9142, 11.5135, 13.9534},
        {7.26943, 8.28372, 9.64631, 10.2867, 12.2556, 13.9916},
        {7.11105, 7.71839, 7.90216, 7.92854, 8.84386, 8.08016},
        {7.1653, 7.79587, 8.09171, 8.38418, 8.08607, 7.18608}};
    //for PbPb/pp ratio sys.
    //QM period uncertainty (priminary)
//    double syspercent[nbin][nrbin]={{1.40033, 4.04487, 4.18342, 5.09741, 11.2329, 7.64036},
//        {1.31925, 3.70484, 3.98242, 5.10289, 10.3913, 7.52044},
//        {0.812068, 2.08369, 3.289, 4.56697, 6.55107, 7.14174},
//        {0.975659, 2.53639, 3.43957, 4.66036, 7.77082, 7.21841}};
//    double syspercent[nbin][nrbin]={{1.52935, 1.32517, 3.68347, 5.27916, 7.37934, 9.3645},
//        {1.39669, 1.43364, 3.81572, 5.47206, 7.40952, 10.2376},
//        {1.06686, 1.15969, 2.81642, 3.29993, 4.85675, 5.19809},
//        {0.931882, 1.13045, 2.58427, 2.75043, 3.00904, 3.75146}};

//    double syspercent[nbin][nrbin]={{1.20498, 2.68224, 5.112, 7.14828, 8.66799, 11.1667},
//        {1.22724, 2.85424, 5.18561, 7.352, 8.94337, 10.992},
//        {1.19124, 1.45967, 3.79709, 3.1998, 5.34887, 4.71617},
//        {1.18007, 1.34272, 2.48922, 2.67718, 2.80794, 3.2557}};
//    //! when using 2011 pp reference
//    double syspercent[nbin][nrbin]={{1.26314, 2.7714, 5.15745, 7.18162, 8.69646, 11.1731},
//        {1.24532, 2.90773, 5.23464, 7.38403, 8.96785, 10.9962},
//        {1.22754, 1.58031, 3.86354, 3.26574, 5.39212, 4.74478},
//        {1.18007, 1.40265, 2.58602, 2.74553, 2.79803, 3.08416}};

    //! when using 2013 pp reference
    double syspercent[nbin][nrbin]={{1.72241, 3.00862, 5.28871, 7.27646, 8.77494, 11.2343},
        {1.70939, 3.13465, 5.36401, 7.4763, 9.04397, 11.0584},
        {1.69647, 1.96686, 4.0371, 3.46932, 5.5178, 4.88714},
        {1.66245, 1.82719, 2.83878, 2.98481, 3.03317, 3.29897}};
    
    //for Rcp ratio sys
//    //QM period uncertainty (priminary)
//    double sysRcppercent[nbin][nrbin]={{1.42016, 2.08759, 2.75502, 3.42245, 4.08988, 4.75731},
//        {0.604485, 1.43505, 2.265615, 3.09618, 3.926744, 4.75731},
//        {0.604485, 1.63792, 2.67136, 3.70479, 4.73822, 5.77166},
//        {0, 0, 0, 0, 0, 0}};
//    double sysRcppercent[nbin][nrbin]={{1.21728, 1.764, 2.73071, 4.1941, 4.74961, 4.96771},
//        {1.02829, 0.928683, 2.86656, 4.33202, 4.58184, 6.22735},
//        {0.327669, 0.39081, 0.89254, 1.33783, 2.93244, 3.78168},
//        {0, 0, 0, 0, 0, 0}};
    
//update with tracking efficiency uncertainty
    double sysRcppercent[nbin][nrbin]={{1.09671, 2.62119, 4.5079, 6.44626, 8.40059, 10.3619},
        {1.07768, 2.65708, 4.63095, 6.6572, 8.6993, 10.7482},
        {1.00336, 1.27668, 1.80713, 2.4319, 3.09439, 3.77478},
        {0, 0, 0, 0, 0, 0}};

    switch (opt) {
        case kDiff:
        //    if(nbin>1)c1->Divide(nbin/2, 2,0, 0);
            //  else c1->Divide(1, 2,0, 0);           
            for(int ipad = 1 ; ipad <=nbin; ipad++){
                c1->cd(ipad);
                c1->cd(ipad)->SetLogy();
                dummy->GetYaxis()->SetTitle("#rho (r)");
                if(trackcut>2.)dummy->SetAxisRange(1.e-3, 20., "Y") ;
                else 
                dummy->SetAxisRange(1.e-4, 100., "Y") ;
                dummy->DrawCopy();
                
                RadiusRho[nbin-ipad] = new TGraphErrors(nfile, rad, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                RadiusRho[nbin-ipad]->SetMarkerStyle(20);
                RadiusRho[nbin-ipad]->SetMarkerColor(1);
                RadiusRho[nbin-ipad]->SetMarkerSize(1.5);
                RadiusRho[nbin-ipad]->SetLineColor(1);
                RadiusRho[nbin-ipad]->Draw("same PE") ;
                
                bkgRadiusRho[nbin-ipad] = new TGraphErrors(nrbin, rad, bkgrho[nbin-ipad], errR, errbkgrho[nbin-ipad]);
                bkgRadiusRho[nbin-ipad]->SetMarkerStyle(29);
                bkgRadiusRho[nbin-ipad]->SetMarkerColor(1);
                bkgRadiusRho[nbin-ipad]->SetMarkerSize(1.5);
                bkgRadiusRho[nbin-ipad]->SetLineColor(1);
                bkgRadiusRho[nbin-ipad]->Draw("same PE") ;
                if(lev==kReco){
                ppsmearRadiusRho[nbin-ipad] = new TGraphErrors(nfile, rad, ppsmearrho[nbin-ipad], errR, errppsmearrho[nbin-ipad]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                ppsmearRadiusRho[nbin-ipad]->SetMarkerStyle(24);
                ppsmearRadiusRho[nbin-ipad]->SetMarkerColor(2);
                ppsmearRadiusRho[nbin-ipad]->SetMarkerSize(1.5);
                ppsmearRadiusRho[nbin-ipad]->SetLineColor(2);
                ppsmearRadiusRho[nbin-ipad]->Draw("same PE") ;
                
                ppsmearbkgRadiusRho[nbin-ipad] = new TGraphErrors(nrbin, rad, ppsmearbkgrho[nbin-ipad], errR, errppsmearbkgrho[nbin-ipad]);
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
                    //       if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.45,0.85,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.78,17);
                    if(ipad==3)drawText(Form("%s jet (anti-k_{T},R=0.3)",level.Data()),0.25,0.85,16);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.78,17);
                }
                else {
                    drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.6,0.75,17);
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
            if(SavePlot)c1->Print(Form("%s/%s%s%sJetPt%.f_%.fDiffJSwith%sBkg.gif",plotsdir.Data(),type.Data(), level.Data(),sample.Data(), pt[current],pt[current+1],bkg.Data()));  
            
            if(SaveFile){
                if(IsMC)
                    TFile * outf = new TFile(Form("%s%s%sIncJet%sDiffJSbkg.root", type.Data(), sample.Data(), level.Data(), bkg.Data()), "UPDATE");
                else
                    TFile * outf = new TFile(Form("%s%sIncJet%sDiffJSbkg.root", type.Data(), level.Data(), bkg.Data()), "UPDATE");
                ppRadiusRho->Write(Form("ppRhoJetPt%.f_%.f_Cen%d-%d", pt[current],pt[current+1], centr[0],centr[nbin]));
                ppbkgRadiusRho->Write(Form("ppbkgRhoJetPt%.f_%.f_Cen%d-%d", pt[current],pt[current+1], centr[0],centr[nbin]));
                for(int ibin = 0 ; ibin <nbin; ibin++){                    
                    RadiusRho[ibin]->Write(Form("PbPbRhoJetPt%.f_%.f_Cen%d-%d", pt[current],pt[current+1], centr[ibin],centr[ibin+1]));
                    bkgRadiusRho[ibin]->Write(Form("PbPbbkgRhoJetPt%.f_%.f_Cen%d-%d", pt[current],pt[current+1], centr[ibin],centr[ibin+1]));
                    if(lev==kReco){
                    ppsmearRadiusRho[ibin]->Write(Form("ppsmearRhoJetPt%.f_%.f_Cen%d-%d", pt[current],pt[current+1], centr[ibin],centr[ibin+1]));
                    ppsmearbkgRadiusRho[ibin]->Write(Form("ppsmearbkgRhoJetPt%.f_%.f_Cen%d-%d", pt[current],pt[current+1], centr[ibin],centr[ibin+1]));
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
                ppsmearRadiusPsi[nbin-ipad] = new TGraphErrors(nfile, rad, ppsmearpsi[nbin-ipad], errR, errppsmearpsi[nbin-ipad]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                ppsmearRadiusPsi[nbin-ipad]->SetMarkerStyle(20);
                ppsmearRadiusPsi[nbin-ipad]->SetMarkerColor(1);
                ppsmearRadiusPsi[nbin-ipad]->SetLineColor(1);
                ppsmearRadiusPsi[nbin-ipad]->Draw("same PE") ;
                
                ppsmearbkgRadiusPsi[nbin-ipad] = new TGraphErrors(nrbin, rad, ppsmearbkgpsi[nbin-ipad], errR, errppsmearbkgpsi[nbin-ipad]);
                ppsmearbkgRadiusPsi[nbin-ipad]->SetMarkerStyle(22);
                ppsmearbkgRadiusPsi[nbin-ipad]->SetMarkerColor(1);
                ppsmearbkgRadiusPsi[nbin-ipad]->SetLineColor(1);
                ppsmearbkgRadiusPsi[nbin-ipad]->Draw("same PE") ;
                }
                else {
                    ppRadiusPsi->Draw("same PE") ;
                    ppbkgRadiusPsi->Draw("same PE") ;  
                }
                RadiusPsi[nbin-ipad] = new TGraphErrors(nfile, rad, psi[nbin-ipad], errR, errpsi[nbin-ipad]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                RadiusPsi[nbin-ipad]->SetMarkerStyle(24);
                RadiusPsi[nbin-ipad]->SetMarkerColor(2);
                RadiusPsi[nbin-ipad]->SetLineColor(2);
                RadiusPsi[nbin-ipad]->Draw("same PE") ;
                
                bkgRadiusPsi[nbin-ipad] = new TGraphErrors(nrbin, rad, bkgpsi[nbin-ipad], errR, errbkgpsi[nbin-ipad]);
                bkgRadiusPsi[nbin-ipad]->SetMarkerStyle(26);
                bkgRadiusPsi[nbin-ipad]->SetMarkerColor(2);
                bkgRadiusPsi[nbin-ipad]->SetLineColor(2);
                bkgRadiusPsi[nbin-ipad]->Draw("same PE") ;
                                
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[ipad]),0.2,0.65, 15);
                if(ipad==1)drawText("Ak PF, R=0.3",0.6,0.75,15);
                if(nbin>1){
                    if(ipad==2)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.45,0.75,15);
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
                    drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.45,0.75,15);
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
            if(SavePlot)c1->Print(Form("%s/%sJetPt%.f_%.fIntJSwith%Bkg.eps",plotsdir.Data(),type.Data(), pt[current],pt[current+1],bkg.Data()));               
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
                dummy->GetYaxis()->SetTitle("#rho (r)");
              if(trackcut>2.)dummy->SetAxisRange(8.e-3, 50., "Y") ;
                else 
                dummy->SetAxisRange(7.e-2, 70., "Y") ;
                dummy->DrawCopy();
          //      ppRhoNoBkg->Draw("same PE") ;

            //    t1->AddEntry(ppsmearRhoNoBkg[nbin-ipad], Form("%d-%d%%",centr[nbin-ipad],centr[ipad]), "P");
                RhoNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, rhosub[nbin-ipad], errR, errrhosub[nbin-ipad]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                RhoNoBkg[nbin-ipad]->SetMarkerStyle(20);
                RhoNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                RhoNoBkg[nbin-ipad]->SetMarkerColor(1);
                RhoNoBkg[nbin-ipad]->SetLineColor(1);
               if(!IsMC)drawSys(RhoNoBkg[nbin-ipad], sysprofile[nbin-ipad],deltacone/2.,7, 1001, 1);
            //    if(!IsMC)drawSys(RhoNoBkg[nbin-ipad], sysprofile[nbin-ipad],deltacone/2.,kGray+2, 1001, 1);  
               RhoNoBkg[nbin-ipad]->Draw("same PE") ;
                cout <<"icentrality= " <<nbin-ipad <<endl ;
                RhoNoBkg[nbin-ipad]->Print();

                if(lev==kReco){
            //        ppsmearRhoNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, pprhosub, errR, errpprhosub);
                ppsmearRhoNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, ppsmearrhosub[nbin-ipad], errR, errppsmearrhosub[nbin-ipad]);
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
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[nbin-ipad+1]),0.4,0.08, 17);
                if(ipad==2){
                    if(IsMC){
                        t1->AddEntry(RhoNoBkg[nbin-1],"PYTHIA+HYDJET", "PL");
                      //  t1->AddEntry(ppRhoNoBkg,"PYTHIA", "P");
                        if(lev==kReco)t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"PYTHIA reference", "PL"); 
                        else t1->AddEntry(ppRhoNoBkg,"PYTHIA", "PL"); 
                    }
                    else {
                        t1->AddEntry(RhoNoBkg[nbin-1],Form("PbPb"), "PL");
//                   t1->AddEntry(RhoNoBkg[nbin-1],Form("PbPb:L_{int} = %.f #mub^{-1}", pbpbLumi), "PL");
                //        t1->AddEntry(ppRhoNoBkg,"PYTHIA", "PL");
//                     t1->AddEntry(ppRhoNoBkg,"pp measured", "PL");
        //                if(lev==kReco) t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"pp 2013 (wrong JEC)", "PL");
                        if(lev==kReco) t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"pp reference", "PL");
               //         if(lev==kReco)t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"PYTHIA reference", "PL");
//                       if(smf=="MySmear") t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"ak3PF pp:Yaxian", "PL");
//                        else t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"PYTHIA reference", "PL");
                  //      t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"TrkCorr in 4 CentBin", "");
                    }
                    t1->Draw("same");
                }
//                if(ipad==3 ){
//                    if(lev==kReco && IsMC==kFALSE) t2->AddEntry(ppsmearRhoNoBkg[nbin-1],"pp reference", "PL");
//                    t2->Draw("same");
//                 //   drawCMS(0.5,0.85,pbpbLumi);
//                //    drawText(Form(" L_{int} = %.f #mub^{-1}", pbpbLumi),0.01,0.9,17);
//                }
             //   if(ipad==1)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.7,17);
                if(nbin>1){
                    if(IsMC==kFALSE){
                    if(ipad==1) drawText(Form("CMS,  #sqrt{s_{NN}} = 2.76 TeV "),0.25,0.9,17);
                    if(ipad==2) drawCMSpp(0.25, 0.95, ppLumi);
                    if(ipad==3) drawCMS(0.2, 0.95, pbpbLumi);
                //    if(ipad==3) drawCMSpp(0.3,0.65,ppLumi);
                //    if(ipad==4) drawCMS(0.3,0.65,pbpbLumi);
                    }
                    if(ipad==1){
                       drawText("anti-k_{T} ( R =0.3 ), PF Jets",0.25,0.76,17);
                   //     drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.3,0.68,17);
                        drawText(Form("p_{T}^{jet} > %.f GeV/c", pt[current]),0.45,0.68,17);
                       drawText(Form("%.1f < |#eta|^{jet} < %.f",etalimit, etacut),0.55,0.58,17);
                        drawText(Form("p_{T}^{track} >%.f GeV/c",trackcut),0.55,0.48,17);
                    }
//                    if(ipad==2){
//                        if(lev==kGen)drawText(Form("Generator Level"),0.3,0.76,17); 
//                   //     drawText(Form("Tracks in cone (#Delta R < 0.3)"),0.1,0.76,17); 
//                 //       drawText(Form("p_{T}^{track} >%.f GeV/c",trackcut),0.3,0.7,17);
//                    }
//                //    if(ipad==4)drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.8,17);
//                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", pt[current]),0.25,0.8,17);
//                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.7,17);
//                    if(ipad==3)drawText("Ak PF, R=0.3",0.25,0.85,17);
//            //        if(ipad==3)drawText(Form("%s jet(anti-k_{T},R=0.3)",level.Data()),0.2,0.85,17);
//                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.25,0.7,17);
                }
                else {
                    drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.6,0.75,17);
                }
            }

            for(int ibin = 0 ; ibin <nbin; ibin++){
            c1->cd(ibin+nbin+1);
                if(IsMC)dummy->SetAxisRange(0.3, 1.7, "Y") ;
//                if(IsMC)dummy->SetAxisRange(0.82, 1.22, "Y") ;
                else dummy->SetAxisRange(0.3, 1.7, "Y") ;
            //    else dummy->SetAxisRange(0.8, 1.25, "Y") ;
          //      dummy->GetYaxis()->SetTitle("#rho(r)^{reference}/#rho(r)^{measured}");
        //        dummy->GetYaxis()->SetTitle(Form("#rho(r)/#rho(r)^{50-100%%}"));
         //      dummy->GetYaxis()->SetTitle("ratio");
                if(IsMC)
   //                 dummy->GetYaxis()->SetTitle("PYTHIA+HYDJET/PYTHIA");
                    dummy->GetYaxis()->SetTitle("#rho(r)^{PYTHIA+HYDJET}/#rho(r)^{PYTHIA}");
                else
         //           dummy->GetYaxis()->SetTitle("PbPb/pp_{reference}");
                    dummy->GetYaxis()->SetTitle("#rho(r)^{PbPb}/#rho(r)^{pp}");
          //     else dummy->GetYaxis()->SetTitle("Ratio");
           //     else dummy->GetYaxis()->SetTitle("#rho(r)^{PbPb}/#rho(r)^{PYTHIA}");
                dummy->DrawCopy();
                //                dummy->DrawCopy();
                RatioRhoNoBkg[nbin-ibin-1] = new TGraphErrors(nrbin, rad, rhoratio[nbin-ibin-1], errR, errrhoratio[nbin-ibin-1]);
                //       RatioRhoNoBkg[ibin] = new TGraphErrors(3, radRebin, rhoratioRebin[ibin], errR, 0);
                RatioRhoNoBkg[nbin-ibin-1]->SetMarkerStyle(20);
                RatioRhoNoBkg[nbin-ibin-1]->SetMarkerColor(1);
                RatioRhoNoBkg[nbin-ibin-1]->SetLineColor(1);
                RatioRhoNoBkg[nbin-ibin-1]->SetMarkerSize(1.5);
               if(!IsMC)drawSys(RatioRhoNoBkg[nbin-ibin-1], syspercent[nbin-ibin-1],deltacone/2.,7, 1001, 1);
           //     if(!IsMC)drawSys(RatioRhoNoBkg[nbin-ibin-1], syspercent[nbin-ibin-1],deltacone/2.,kGray+2, 1001, 1);  
                RatioRhoNoBkg[nbin-ibin-1]->Draw("same PE") ;
                if(lev==kReco){
                    smearRatioRhoNoBkg[nbin-ibin-1] = new TGraphErrors(nrbin, rad, smearrhoratio[nbin-ibin-1], errR, errsmearrhoratio[nbin-ibin-1]);
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
                regSun(-0.01,1.,0.31,1.,1, 1);
         //       if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.4,0.9, 17);
//                if(nbin>1){
//                    if(ibin==3){
//                        if(IsMC){
//                            if(lev==kReco)t2->AddEntry(RatioRhoNoBkg[nbin-1],"PYTHIA+HYDJET/PYTHIA_ref", "PL");
//                            if(lev==kReco)t2->AddEntry(smearRatioRhoNoBkg[nbin-1],"Reference/Measured", "PL");
//                        }
//                        else {
//                            t2->AddEntry(RatioRhoNoBkg[nbin-1],"PbPb/pp_measure", "PL");
//                           if(lev==kReco) t2->AddEntry(smearRatioRhoNoBkg[nbin-1],"pp_ref/pp_measure", "PL");
////                            t2->AddEntry(RatioRhoNoBkg[nbin-1],"PbPb/PYTHIA", "PL");
////                            t2->AddEntry(smearRatioRhoNoBkg[nbin-1],"PYTHIA ref/PYTHIA", "PL");
//                        }
//          //              if(lev==kReco)t2->Draw("same");
//                    }
//                }
           }
            if(SavePlot)
           //   c1->Print(Form("%sJetPtThres%.fTrkCut%.fDiffJSAfter%sBkgSub.C",type.Data(), pt[current],trackcut, bkg.Data()));          
                c1->Print(Form("%s/%s%sJetPtBin%.f_%.f%sTrkCut%.fEtaCut%.f%sDiffJS%dCentAfter%sBkgSub2013ppSmall.pdf",plotsdir.Data(),type.Data(), level.Data(), pt[current],pt[current+1], effTab.Data(),trackcut, etacut*10, smf.Data(),nbin, bkg.Data(), date));
//             c1->Print(Form("%s/%s%sJetPtThres%.f%sTrkCut%.f%sDiffJS%dCentAfter%sBkgSubFullStat.pdf",plotsdir.Data(),type.Data(), level.Data(), pt[current],effTab.Data(),trackcut, smf.Data(),nbin, bkg.Data(), date));               
         //   c1->Print(Form("%s/%s%sJetSkimCut100PtBin%.f_%.fTrkCut%.fDiffJSAfter%sBkgSub.gif",plotsdir.Data(),type.Data(),level.Data(), pt[current],pt[current+1],trackcut, bkg.Data()));               
            c1->Update();   
            if(SaveFile){
                TFile * outf = new TFile(Form("%s%sIncJet%sTrk%sDiffJSbkgSub.root", type.Data(), level.Data(), bkg.Data(), effTab.Data()), "UPDATE");
                ppRhoNoBkg->Write(Form("ppRhoERbkgsubJetPt%.f_%.f_Cen%d-%d", pt[current],pt[current+1],centr[0],centr[nbin-1]));
                for(int ibin = 0 ; ibin <nbin; ibin++){                    
                    RhoNoBkg[ibin]->Write(Form("PbPbRhoERbkgsubJetPt%.f_%.f_Cen%d-%d", pt[current],pt[current+1], centr[ibin],centr[ibin+1]));
                    RatioRhoNoBkg[ibin]->Write(Form("RatioRhoERbkgsubJetPt%.f_%.f_Cen%d-%d", pt[current],pt[current+1], centr[ibin],centr[ibin+1]));
                    if(lev==kReco)ppsmearRhoNoBkg[ibin]->Write(Form("ppsmearRhoERbkgsubJetPt%.f_%.f_Cen%d-%d", pt[current],pt[current+1], centr[ibin],centr[ibin+1]));
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
                RhoNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, rhosub[nbin-ipad], errR, errrhosub[nbin-ipad]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                RhoNoBkg[nbin-ipad]->SetMarkerStyle(20);
                RhoNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                RhoNoBkg[nbin-ipad]->SetMarkerColor(1);
                RhoNoBkg[nbin-ipad]->SetLineColor(1);
                if(!IsMC)drawSys(RhoNoBkg[nbin-ipad], sysprofile[nbin-ipad],deltacone/2.,7, 1001, 1);  
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
                    //         if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.25,0.8,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.7,17);
                    if(ipad==3)drawText("Ak PF, R=0.3",0.55,0.9,17);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.82,17);
                }
                else {
                    drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.6,0.75,17);
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

                RcpRatioRhoNoBkg[nbin-ibin-1] = new TGraphErrors(nrbin, rad, rcprhoratio[nbin-ibin-1], errR, errrcprhoratio[nbin-ibin-1]);
                //       RatioRhoNoBkg[ibin] = new TGraphErrors(3, radRebin, rhoratioRebin[ibin], errR, 0);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetMarkerStyle(20);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetMarkerColor(1);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetLineColor(1);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetMarkerSize(1.5);
                if(!IsMC)drawSys(RcpRatioRhoNoBkg[nbin-ibin-1], sysRcppercent[nbin-ibin-1],deltacone/2.,7, 1001, 1);  
                RcpRatioRhoNoBkg[nbin-ibin-1]->Draw("same PE") ;
                RcpRatioRhoNoBkg[nbin-ibin-1]->Print();
                regSun(0.,1.,0.3,1.,1, 1);
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.25,0.9, 17);
            }            
            
            if(SavePlot)
                c1->Print(Form("%s/%sJetPtThres%.fTrkCut%.fDiffJSRcpAfter%sBkgSubFullStat.pdf",plotsdir.Data(),type.Data(), pt[current],trackcut, bkg.Data()));               
             //   c1->Print(Form("%s/%sJetPtThres%.fTrkCut%.fDiffJSRcpAfter%sBkgSub.gif",plotsdir.Data(),type.Data(), pt[current],trackcut, bkg.Data()));               
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
                
                RcpRatioRhoNoBkg[nbin-ibin-1] = new TGraphErrors(nrbin, rad, rcprhoratio[nbin-ibin-1], errR, errrcprhoratio[nbin-ibin-1]);
                //       RatioRhoNoBkg[ibin] = new TGraphErrors(3, radRebin, rhoratioRebin[ibin], errR, 0);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetMarkerStyle(20);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetMarkerColor(1);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetLineColor(1);
                RcpRatioRhoNoBkg[nbin-ibin-1]->SetMarkerSize(1.5);
                if(!IsMC)drawSys(RcpRatioRhoNoBkg[nbin-ibin-1], sysRcppercent[nbin-ibin-1],deltacone/2.,7, 1001, 1);  
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
                    //         if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.8,17);
                    if(ibin==3)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.45,0.82,17);
                    if(ibin==3)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.75,17);
                    if(ibin==2)drawText("Ak PF, R=0.3",0.55,0.8,17);
                    if(ibin==2)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.75,17);
                }
                else {
                    drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.6,0.75,17);
                }

            }            
            
            if(SavePlot)
                //      c1->Print(Form("%s/%sJetPtThres%.f%.sTrkCut%.fIncDiffJS%sRatioSmearAfterBkgSub.eps",plotsdir.Data(),type.Data(), leadingjetcut,bkg.Data(),trackcut,Norm.Data()));               
                c1->Print(Form("%s/%sJetPtThres%.fTrkCut%.fDiffJSRcpRatio.pdf",plotsdir.Data(),type.Data(), pt[current],trackcut, bkg.Data()));               
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

                RatioRhoNoBkg[nbin-ibin-1] = new TGraphErrors(nrbin, rad, rhoratio[nbin-ibin-1], errR, errrhoratio[nbin-ibin-1]);
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
                    //         if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.8,17);
                    if(ibin==1)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.25,0.9,17);
                    if(ibin==1)drawText(Form("|#eta|_{jet} < %.f", etacut),0.25,0.83,17);
                    if(ibin==1)drawText("Ak PF, R=0.3",0.25,0.75,17);
                    if(ibin==1)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.25,0.7,17);
                }
                else {
                    drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.6,0.75,17);
                }

            }            
            
            if(SavePlot)
               //     c1->Print(Form("%sJetPtThres%.fTrkCut%.fDiffJSRatio.C",type.Data(), pt[current],trackcut, bkg.Data()));
                c1->Print(Form("%s/%sJetPtThres%.fTrkCut%.fDiffJSRatio.gif",plotsdir.Data(),type.Data(), pt[current],trackcut, bkg.Data()));               
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
                    ppsmearRhoNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, ppsmearrhosub[nbin-ipad], errR, errppsmearrhosub[nbin-ipad]);
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
                        t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"pp smeared", "PL");
//                        t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"pp reference", "PL");
//                        if(smf=="MySmear")t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"pp smeared Yaxian", "PL");
//                        else t1->AddEntry(ppsmearRhoNoBkg[nbin-1],"pp smeared Pawan", "PL");
                    }
                    t1->Draw("same");
                }
                if(nbin>1){
                    //      if(ipad==4)drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.25,0.8,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.7,17);
                    if(ipad==3)drawText(Form("ak3PF jet,R=0.3",level.Data()),0.2,0.85,17);
//                    if(ipad==3)drawText(Form("%s jet(anti-k_{T},R=0.3)",level.Data()),0.2,0.85,17);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.7,17);
                }
                else {
                    drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.6,0.75,17);
                }
            }
            for(int ibin = 0 ; ibin <nbin; ibin++){
                c1->cd(ibin+nbin+1);
                dummy->SetAxisRange(0.9, 1.12, "Y") ;
        //        dummy->GetYaxis()->SetTitle("#rho(r)^{reference}/#rho(r)^{measured}");
                dummy->GetYaxis()->SetTitle("#rho(r)^{smearedOnly}/#rho(r)^{measured}");
                dummy->DrawCopy();
                if(lev==kReco){
                    smearRatioRhoNoBkg[nbin-ibin-1] = new TGraphErrors(nrbin, rad, smearrhoratio[nbin-ibin-1], errR, errsmearrhoratio[nbin-ibin-1]);
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
                //  c1->Print(Form("%sJetPtThres%.fTrkCut%.fDiffJSAfter%sBkgSub.C",type.Data(), pt[current],trackcut, bkg.Data()));          
                c1->Print(Form("%s/%s%sJetPtThres%.fTrkCut%.f%sSmearedOnlyDiffJSAfter%sBkgSub.pdf",plotsdir.Data(),type.Data(), level.Data(), pt[current],trackcut, smf.Data(),bkg.Data()));               
            //  c1->Print(Form("%s/%s%sJetPtBin%.f_%.fTrkCut%.fDiffJSAfter%sBkgSubppRef.gif",plotsdir.Data(),type.Data(),level.Data(), pt[current],pt[current+1],trackcut, bkg.Data()));               
            c1->Update();   
            break ;
        case kCompJSRatioTh:
                c1->cd(1);
                if(IsMC)dummy->SetAxisRange(0.5, 1.5, "Y") ;
                else dummy->SetAxisRange(0., 2.1, "Y") ;
                if(IsMC)
                    dummy->GetYaxis()->SetTitle("#rho(r)^{PYTHIA+HYDJET}/#rho(r)^{PYTHIA}");
                else dummy->GetYaxis()->SetTitle("#rho(r)^{PbPb}/#rho(r)^{pp}");
                //       else dummy->GetYaxis()->SetTitle("#rho(r)^{PbPb}/#rho(r)^{PYTHIA_reference}");
                dummy->DrawCopy();
                
                RatioRhoNoBkg[0] = new TGraphErrors(nrbin, ratio, rhoratio[0], errR, errrhoratio[0]);
                //       RatioRhoNoBkg[ibin] = new TGraphErrors(3, radRebin, rhoratioRebin[ibin], errR, 0);
                RatioRhoNoBkg[0]->SetMarkerStyle(20);
                RatioRhoNoBkg[0]->SetMarkerColor(1);
                RatioRhoNoBkg[0]->SetLineColor(1);
                RatioRhoNoBkg[0]->SetMarkerSize(1.5);
                if(!IsMC)drawSys(RatioRhoNoBkg[0], syspercent[0],deltacone/2.,7, 1001, 1);  
                //     if(!IsMC)drawSys(RatioRhoNoBkg[nbin-ibin-1], syspercent[nbin-ibin-1],deltacone/2.,kGray+2, 1001, 1);  
                RatioRhoNoBkg[0]->Draw("same PE") ;
                regSun(0.,1.,1.0,1.,1, 1);
                drawCMS(0.7,0.3,pbpbLumi);

            const int Nbin = 39 ;
            double x[Nbin], y50[Nbin], y100[Nbin], y200[Nbin];
//            TH1D * ThJS100 = new TH1D("ThJS100", "ThJS100", 40, 0., 1.);
//            ThJS100->Sumw2();
//            TH1D * ThJS200 = new TH1D("ThJS200", "ThJS200", 40, 0., 1.);
//            ThJS200->Sumw2();
            ifstream indata50;
            ifstream indata100;
            ifstream indata200;
            indata50.open("JSRatioIvanE50.txt");
            indata100.open("JSRatioIvanE100.txt");
            indata200.open("JSRatioIvanE200.txt");       
             
            for(Int_t i = 0 ; i <Nbin; i++){
                indata50 >>x[i];
                indata50 >>y50[i];
                indata100 >>x[i];
                indata100 >>y100[i];
                indata200 >>x[i];
                indata200 >>y200[i]; 
//                ThJS100->Fill(x100[i], y100[i]);
//                ThJS200->Fill(x200[i], y200[i]);
                
            }
            indata50.close();
            indata100.close();
            indata200.close();
            TGraph * ThJS50 = new TGraph(Nbin, x, y50);
            TGraph * ThJS100 = new TGraph(Nbin, x, y100);
            TGraph * ThJS200 = new TGraph(Nbin, x, y200);

            ThJS50->SetLineColor(6);
            ThJS50->Draw("same L");
            ThJS100->SetLineColor(2);
            ThJS100->Draw("same L");
            ThJS200->SetLineColor(4);
            ThJS200->Draw("same L");
            
            t1->AddEntry(RatioRhoNoBkg[0],Form("CMS #sqrt{s}=2.76 TeV, %d-%d%%:p_{T}^{jet}>%.f GeV/c", centr[0],centr[1],leadingjetcut), "PL");
            t1->AddEntry(ThJS50,"Vitev, #sqrt{s}=5.5 TeV, E_{T}=50GeV, b=3fm", "L");
            t1->AddEntry(ThJS100,"Vitev, #sqrt{s}=5.5 TeV, E_{T}=100GeV, b=3fm", "L");
            t1->AddEntry(ThJS200,"Vitev, #sqrt{s}=5.5 TeV, E_{T}=200GeV, b=3fm", "L");            
            t1->Draw("same");
            drawText(Form("I.Vitev,et.al, JHEP11 (2008)093"),0.2,0.2,17);

//                if(nbin>1){
//                    //         if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.8,17);
//                    if(ibin==1)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.25,0.9,17);
//                    if(ibin==1)drawText(Form("|#eta|_{jet} < %.f", etacut),0.25,0.83,17);
//                    if(ibin==1)drawText("Ak PF, R=0.3",0.25,0.75,17);
//                    if(ibin==1)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.25,0.7,17);
//                }
//                else {
//                    drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.6,0.75,17);
//                }
//                
          
            
            if(SavePlot)
             //   c1->Print(Form("%sJetPtThres%.fTrkCut%.fDiffJSRatioWithTheory.gif",type.Data(), pt[current],trackcut, bkg.Data()));               
              c1->Print(Form("%s/%sJetPtThres%.fTrkCut%.fDiffJSRatioWithTheory.pdf",plotsdir.Data(),type.Data(), pt[current],trackcut, bkg.Data()));               
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
                PsiNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, psisub[nbin-ipad], errR, errpsisub[nbin-ipad]);
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
                    ppsmearPsiNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, ppsmearpsisub[nbin-ipad], errR, errppsmearpsisub[nbin-ipad]);
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
                    //     if(ipad==4)drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", pt[current]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.7,17);
                    if(ipad==3)drawText("Ak PF, R=0.3",0.25,0.85,17);
                    //        if(ipad==3)drawText(Form("%s jet(anti-k_{T},R=0.3)",level.Data()),0.2,0.85,17);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.25,0.7,17);
                }
                else {
                    drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.6,0.75,17);
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
                RatioPsiNoBkg[nbin-ibin-1] = new TGraphErrors(nrbin, rad, psiratio[nbin-ibin-1], errR, errpsiratio[nbin-ibin-1]);
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
                //  c1->Print(Form("%sJetPtThres%.fTrkCut%.fDiffJSAfter%sBkgSub.C",type.Data(), pt[current],trackcut, bkg.Data()));          
                c1->Print(Form("%s/%s%sJetPtThres%.fEff4CentBinTrkCut%.f%sIntJS%dCentAfter%sBkgSubFullStat.pdf",plotsdir.Data(),type.Data(), level.Data(), pt[current],trackcut, smf.Data(),nbin, bkg.Data()));               
            //   c1->Print(Form("%s/%s%sJetPtBin%.f_%.fTrkCut%.fDiffJSAfter%sBkgSub.gif",plotsdir.Data(),type.Data(),level.Data(), pt[current],pt[current+1],trackcut, bkg.Data()));               
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
//    cms->Draw();
    TLatex *lumi = new TLatex(px,py-0.05,Form("PbPb, #intL dt = %.1f #mub^{-1}",nLumi));
    lumi->SetTextFont(63);
    lumi->SetTextSize(17);
    lumi->SetNDC();
    if(!IsMC)lumi->Draw();
}

void drawCMSpp(float px, float py, float nLumi) {
   TLatex *cms = new TLatex(px,py,"CMS Preliminary");
    cms->SetTextFont(63);
    cms->SetTextSize(17);
    cms->SetNDC();
//    cms->Draw();
    TLatex *lumi = new TLatex(px,py-0.05,Form("pp, #intL dt = %.1f pb^{-1}",nLumi));
//    TLatex *lumi = new TLatex(px,py-0.05,Form("#intL dt = %.1f nb^{-1}",nLumi));
    lumi->SetTextFont(63);
    lumi->SetTextSize(17);
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
    h->SetTitleSize(22);
    h->SetLabelSize(22,"Y");
    h->SetLabelSize(22,"X");
    h->SetTitleFont(44,"Y");
    h->SetTitleSize(20,"Y");
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
