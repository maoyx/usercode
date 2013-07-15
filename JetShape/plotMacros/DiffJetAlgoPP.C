/*
 *  Compare JS with different statistics  
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
#include <iostream.h>
#include <stdlib.h>
#include <fstream.h>
#include <iomanip>
using namespace std;

//define the kinematics cuts for the analysis
const double conesize = 0.3 ;
const double deltacone = 0.05 ;
const double etacut = 2.0 ;
const double dphicut = 7*TMath::TwoPi()/8. ;
const double leadingjetcut = 100. ;
const double subleadjetcut = 40. ;
const double trackcut = 1.0;
bool SavePlot = kTRUE ;
const bool SaveFile = kFALSE ;

const int nfile = 6 ; //8 until 0.3, 13 until 0.5
//const int nfile = (conesize-deltacone/2.)/deltacone +1 ;
//const double rad[] = {0.05, 0.15, 0.25};
//const double ratio[] = {0.05/conesize, 0.15/conesize, 0.25/conesize};
//const double errR[] = {0.05, 0.05, 0.05};
double rad[nfile];
double ratio[nfile];
double errR[nfile]; 
const int nptbin = 6 ;
//const double pt[]={100., 120.,160., 300., 500.};
//const double pt[]={100., 120., 140., 160., 200., 300., 500.};
const double pt[]={100., 300.};
//const double pt[]={100., 120., 150., 200., 300.};

const int nrbin = 6 ;
double rbin[nrbin+1]={0.0,0.05,0.10, 0.15, 0.20,0.25,0.30};

const Double_t kw = 1000;
const Double_t kh = 560;

const Double_t pbpbLumi = 150 ;
const Double_t ppLumi = 5.3 ;
//const Double_t ppLumi = 231 ;

enum Option_t {kDiff, kInt, kDiffSub, kIntSub} ; 
enum Display_t {kMea, kSmear, kRewt} ;
//---------------------------------------------------------------------

void DiffJetAlgoPP(Option_t opt=kDiff, Display_t data = kMea, int current=0){
    for(int i = 0 ; i <nrbin; i++){
        //        if(opt==kInt) rad[i]=rbin[i+1] ;
        //        else         
        rad[i]=deltacone/2.+i*deltacone ;
        ratio[i]=rad[i]/conesize;
        errR[i]=deltacone/2. ;
    }
    TString plotsdir = "/Users/ymao/group/CMS/plots/finalJS/CorrectDR";
    TDatime d ;
    int date = d.GetDate();  
    TString type ;
    if(data==kMea) type="Measure";
    else if(data==kSmear)type = "Smeared" ;
    else type ="Rewted";
    
//    TString bkg ;
//    if(met==kER) bkg = "EtaReflect";
//    else bkg = "RandomCone";
//    const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/FullStatPbPb/RebinJetPt" ;      
    const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/CorrectdR/NoMatchGen/NoTrkAlgoCut/MyselfTrkEff" ;      
    TString fileName ;
    TFile * f ;
    TString ak3fileName ;
    TFile * ak3f ;
    TString effTab = "HistIterTrkCorrtest";  //"HistIterTrkCorrv14fXSecMBEff"; // "HistIterTrkCorrv14fXSec"; // "HistIterTrkCorrtest" ; //
    
    if(data==kMea){
        ////if it is pp, no centrality bins, only one
        const int nbin = 1 ;
        const int centr[] ={0,100};  
    }
    else {
        // for HI centrality bin
        const int nbin = 4 ;
        const int centr[] ={0, 10,30,50,100}; 
    }
    
    TH2F * IntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * bkgIntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    
    TH2F * DiffJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * bkgDiffJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    
    //    TProfile * differential[nbin][nrbin] ;
    //    TProfile * Integrated[nbin][nrbin] ;
    
    TH1F * Integrated[nbin][nrbin] ;
    TH1F * bkgIntegrated[nbin][nrbin] ;
    
    TH1F * differential[nbin][nrbin] ;
    TH1F * bkgdifferential[nbin][nrbin] ;
    
    TGraphErrors * RadiusRho[nbin];
    TGraphErrors * RadiusPsi[nbin];
    
    TGraphErrors * bkgRadiusRho[nbin];
    TGraphErrors * bkgRadiusPsi[nbin];
    
    TH2F * ak3IntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * ak3bkgIntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * ak3DiffJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * ak3bkgDiffJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    
    TH1F * ak3Integrated[nbin][nrbin] ;
    TH1F * ak3bkgIntegrated[nbin][nrbin] ;
    
    TH1F * ak3diff[nbin][nrbin] ;
    TH1F * ak3bkgdiff[nbin][nrbin] ;
    
    TGraphErrors * ak3Rho[nbin];
    TGraphErrors * ak3Psi[nbin] ;
    TGraphErrors * bkgak3Rho[nbin];
    TGraphErrors * bkgak3Psi[nbin] ;

    TGraphErrors * RhoNoBkg[nbin];
    TGraphErrors * PsiNoBkg[nbin] ;
    TGraphErrors * ak3RhoNoBkg[nbin];
    TGraphErrors * ak3PsiNoBkg[nbin] ;

    TGraphErrors * RatioRho[nbin];
    TGraphErrors * RatioPsi[nbin] ;

    TH1::SetDefaultSumw2();
    
    double rho[nbin][nrbin];
    double errrho[nbin][nrbin];
    double psi[nbin][nrbin];
    double errpsi[nbin][nrbin];
    
    double ak3rho[nbin][nrbin];
    double errak3rho[nbin][nrbin];
    double ak3psi[nbin][nrbin];
    double errak3psi[nbin][nrbin];
    
    double bkgrho[nbin][nrbin];
    double errbkgrho[nbin][nrbin];
    double bkgpsi[nbin][nrbin];
    double errbkgpsi[nbin][nrbin];
    
    double ak3bkgrho[nbin][nrbin];
    double errak3bkgrho[nbin][nrbin];
    double ak3bkgpsi[nbin][nrbin];
    double errak3bkgpsi[nbin][nrbin];
    
    double rhosub[nbin][nrbin];
    double errrhosub[nbin][nrbin];
    double psisub[nbin][nrbin];
    double errpsisub[nbin][nrbin];
    
    double ak3rhosub[nbin][nrbin];
    double errak3rhosub[nbin][nrbin];
    double ak3psisub[nbin][nrbin];
    double errak3psisub[nbin][nrbin];


    double rhoratio[nbin][nrbin];
    double errrhoratio[nbin][nrbin];
    double psiratio[nbin][nrbin];
    double errpsiratio[nbin][nrbin];
    
    double scalepsi[nbin];
    double scaleak3psi[nbin];
    double scalerho[nbin];
    double scaleak3rho[nbin];
    
    TString canv_name = "c1";
    //    canv_name += parti;
    //    canv_name += cent;
    //    makeMultiPanelCanvas(c1,4,1,0.0,0.0,0.2,0.2,0.02);
    if(data==kMea){
        c1 = new TCanvas(canv_name," ",800, 800);
        if(opt==kDiffSub) c1->Divide(1, 2, 0, 0);  
    }
    else {
        c1 = new TCanvas(canv_name," ",10,10,kw,kh);
        if(opt==kDiffSub)
            makeMultiPanelCanvas(c1,4,2,0.0,0.0,0.2,0.2,0.02);
        else 
            makeMultiPanelCanvas(c1,4,1,0.0,0.0,0.2,0.2,0.02);   
    }
    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);
    
    
    if(data==kMea)
        TLegend *t1=new TLegend(0.65,0.8,0.8,0.92);
    else 
        TLegend *t1=new TLegend(0.35,0.75,0.8,0.93);
    t1->SetFillColor(0);
    t1->SetBorderSize(0);
    t1->SetFillStyle(0);
    t1->SetTextFont(63);
    t1->SetTextSize(17);
    TLegend *t2=new TLegend(0.20,0.45,0.35,0.6);
    t2->SetFillColor(0);
    t2->SetBorderSize(0);
    t2->SetFillStyle(0);
    t2->SetTextFont(63);
    t2->SetTextSize(17);
    // open the data file
    switch(data){
        case kMea:
//            fileName=Form("DATAPP2011_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_pp_merged_full.root", effTab.Data(), trackcut);
//         //  fileName = Form("DATAPP_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82.root", effTab.Data(), trackcut);
//            ak3fileName = Form("DATAPP_Ak3PFIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut);
            fileName=Form("DATAPP_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut);
            ak3fileName = Form("DATAPP_AkPu3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut);
//            fileName=Form("mergedCSdiff_MCPP2013_Ak3PFIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_ppHiIterativeTrack_merged.root", effTab.Data(), trackcut);
//            ak3fileName = Form("mergedCSdiff_MCPP2013_ShiftVertexAk3PFIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_ppHiIterativeTrack_merged.root", effTab.Data(), trackcut); 
            break ;
        case kSmear:
            fileName=Form("DATAPP_PFIncJetPt100_HistCorrTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1_pp_merged_full.root", trackcut); 
            ak3fileName = Form("DATAPP_Ak3PFIncAk3PFJetPt100_HistCorrTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1_pp_merged_full.root", trackcut);  
            break ;
        case kRewt:
            fileName=Form("DATAPP_Ak3PFwtIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut, nbin); 
            ak3fileName = Form("DATAPP_Ak3PFwtIncJetYaxianParaPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut, nbin);
            
//            fileName=Form("DATAPP2011_Ak3PFwtIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_pp_merged_full.root", effTab.Data(), trackcut);
//            //  fileName = Form("DATAPP_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82.root", effTab.Data(), trackcut);
//            ak3fileName = Form("DATAPP_Ak3PFwtIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_PP2013_HiForest_PromptReco_JsonPP_Jet80_HIReco_forestv84_v2.root", effTab.Data(), trackcut);
            break ;
    }
    
 //   f = TFile::Open(Form("/Users/ymao/group/CMS/anaOutputs/JetShape/SameTrkBinFrank/%s", fileName.Data()), "readonly");
    f = TFile::Open(Form("%s/%s", kHomeDir, fileName.Data()), "readonly");
    ak3f = TFile::Open(Form("%s/%s", kHomeDir, ak3fileName.Data()), "readonly");
  //  ak3f = TFile::Open(Form("/Users/ymao/group/CMS/anaOutputs/JetShape/FullStatPbPb/%s", ak3fileName.Data()), "readonly");
    for(int ibin = 0 ; ibin <nbin; ibin++){   
        scaleak3psi[ibin]=0. ;
        scalepsi[ibin]=0. ;
        scaleak3rho[ibin]=0. ;
        scalerho[ibin]=0. ;
        double sumrho = 0.0 ;
        for(int ir =0 ; ir <nrbin; ir++){  
            rho[ibin][ir] = 0 ;
            errrho[ibin][ir] = 0 ;
            psi[ibin][ir]=0 ;
            errpsi[ibin][ir]=0 ;  
            ak3rho[ibin][ir] = 0 ;
            errak3rho[ibin][ir] = 0 ;
            ak3psi[ibin][ir]=0 ;
            errak3psi[ibin][ir]=0 ;  
            bkgrho[ibin][ir] = 0 ;
            errbkgrho[ibin][ir] = 0 ;
            bkgpsi[ibin][ir]=0 ;
            errbkgpsi[ibin][ir]=0 ;  
            ak3bkgrho[ibin][ir] = 0 ;
            errak3bkgrho[ibin][ir] = 0 ;
            ak3bkgpsi[ibin][ir]=0 ;
            errak3bkgpsi[ibin][ir]=0 ;  
            rhosub[ibin][ir] = 0 ;
            errrhosub[ibin][ir] = 0 ;
            psisub[ibin][ir]=0 ;
            errpsisub[ibin][ir]=0 ;  
            ak3rhosub[ibin][ir] = 0 ;
            errak3rhosub[ibin][ir] = 0 ;
            ak3psisub[ibin][ir]=0 ;
            errak3psisub[ibin][ir]=0 ;  
            
            rhoratio[ibin][ir] = 0 ;
            errrhoratio[ibin][ir] = 0 ;
            psiratio[ibin][ir] = 0 ;
            errpsiratio[ibin][ir] = 0 ;
            
            ak3IntJS[ibin][ir]=(TH2F*)ak3f->Get(Form("IntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            IntJS[ibin][ir]=(TH2F*)f->Get(Form("IntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            ak3bkgIntJS[ibin][ir]=(TH2F*)ak3f->Get(Form("bkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            bkgIntJS[ibin][ir]=(TH2F*)f->Get(Form("bkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            
            ak3DiffJS[ibin][ir]=(TH2F*)ak3f->Get(Form("differentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            DiffJS[ibin][ir]=(TH2F*)f->Get(Form("differentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            ak3bkgDiffJS[ibin][ir]=(TH2F*)ak3f->Get(Form("bkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            bkgDiffJS[ibin][ir]=(TH2F*)f->Get(Form("bkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            
            ak3IntJS[ibin][ir]->Sumw2();
            ak3Integrated[ibin][ir]=(TH1F*)ak3IntJS[ibin][ir]->ProjectionY(Form("qJetPt%.f_%.fPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),ak3IntJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), ak3IntJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            
            IntJS[ibin][ir]->Sumw2();
            Integrated[ibin][ir]=(TH1F*)IntJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),IntJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), IntJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            
            ak3bkgIntJS[ibin][ir]->Sumw2();
            ak3bkgIntegrated[ibin][ir]=(TH1F*)ak3bkgIntJS[ibin][ir]->ProjectionY(Form("qJetPt%.f_%.fbkgPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),ak3bkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), ak3bkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            bkgIntJS[ibin][ir]->Sumw2();
            bkgIntegrated[ibin][ir]=(TH1F*)bkgIntJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fbkgPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),bkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), bkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            
            ak3DiffJS[ibin][ir]->Sumw2();
            ak3diff[ibin][ir]=(TH1F*)ak3DiffJS[ibin][ir]->ProjectionY(Form("qJetPt%.f_%.fRhodR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),ak3DiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), ak3DiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            
            DiffJS[ibin][ir]->Sumw2();
            differential[ibin][ir]=(TH1F*)DiffJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fRhodR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),DiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), DiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            
            ak3bkgDiffJS[ibin][ir]->Sumw2();
            ak3bkgdiff[ibin][ir]=(TH1F*)ak3bkgDiffJS[ibin][ir]->ProjectionY(Form("qJetPt%.f_%.fbkgRhodR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),ak3bkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), ak3bkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            bkgDiffJS[ibin][ir]->Sumw2();
            bkgdifferential[ibin][ir]=(TH1F*)bkgDiffJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fbkgRhodR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),bkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), bkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            
            psi[ibin][ir]=Integrated[ibin][ir]->GetMean(1);
            errpsi[ibin][ir]=Integrated[ibin][ir]->GetMeanError(1);
            bkgpsi[ibin][ir]=bkgIntegrated[ibin][ir]->GetMean(1);
            errbkgpsi[ibin][ir]=bkgIntegrated[ibin][ir]->GetMeanError(1);
            rho[ibin][ir]=differential[ibin][ir]->GetMean(1);
            errrho[ibin][ir]=differential[ibin][ir]->GetMeanError(1);
            bkgrho[ibin][ir]=bkgdifferential[ibin][ir]->GetMean(1);
            errbkgrho[ibin][ir]=bkgdifferential[ibin][ir]->GetMeanError(1);
            
            ak3psi[ibin][ir]=ak3Integrated[ibin][ir]->GetMean(1);
            errak3psi[ibin][ir]=ak3Integrated[ibin][ir]->GetMeanError(1);
            ak3bkgpsi[ibin][ir]=ak3bkgIntegrated[ibin][ir]->GetMean(1);
            errak3bkgpsi[ibin][ir]=ak3bkgIntegrated[ibin][ir]->GetMeanError(1);
            ak3rho[ibin][ir]=ak3diff[ibin][ir]->GetMean(1);
            errak3rho[ibin][ir]=ak3diff[ibin][ir]->GetMeanError(1);
            ak3bkgrho[ibin][ir]=ak3bkgdiff[ibin][ir]->GetMean(1);
            errak3bkgrho[ibin][ir]=ak3bkgdiff[ibin][ir]->GetMeanError(1);

            //nromalize to delta cone size to consistent for definition
            rho[ibin][ir]/=deltacone ;
            ak3rho[ibin][ir]/=deltacone ;
            bkgrho[ibin][ir]/=deltacone ;
            ak3bkgrho[ibin][ir]/=deltacone ;
            errrho[ibin][ir]/=deltacone ;
            errak3rho[ibin][ir]/=deltacone ;
            errbkgrho[ibin][ir]/=deltacone ;
            errak3bkgrho[ibin][ir]/=deltacone ;
//            if(ir==0){
//                rho[ibin][ir]=rho[ibin][ir]*1.1;
//                errrho[ibin][ir]=errrho[ibin][ir]*1.1;
//            }
            
            //calculate the bkg subtracted JS
            rhosub[ibin][ir]=rho[ibin][ir]-bkgrho[ibin][ir];
            errrhosub[ibin][ir]=TMath::Sqrt(TMath::Power(errrho[ibin][ir],2)+TMath::Power(errbkgrho[ibin][ir],2));
//            errrhosub[ibin][ir]=errrho[ibin][ir]+errbkgrho[ibin][ir];
            psisub[ibin][ir]=psi[ibin][ir]-bkgpsi[ibin][ir];
            errpsisub[ibin][ir]=TMath::Sqrt(TMath::Power(errpsi[ibin][ir],2)+TMath::Power(errbkgpsi[ibin][ir],2));
            
            ak3rhosub[ibin][ir]=ak3rho[ibin][ir]-ak3bkgrho[ibin][ir];
            errak3rhosub[ibin][ir]=TMath::Sqrt(TMath::Power(errak3rho[ibin][ir],2)+TMath::Power(errak3bkgrho[ibin][ir],2));
            ak3psisub[ibin][ir]=ak3psi[ibin][ir]-ak3bkgpsi[ibin][ir];
            errak3psisub[ibin][ir]=TMath::Sqrt(TMath::Power(errak3psi[ibin][ir],2)+TMath::Power(errak3bkgpsi[ibin][ir],2));

            scaleak3rho[ibin]+=ak3rhosub[ibin][ir]*deltacone;
            scalerho[ibin]+=rhosub[ibin][ir]*deltacone;                    
            
            //                    sumrho+=ak3Rhosub[ibin][ir]*deltacone;
            //                    cout <<"ir = "<<ir<<"psi =" << psisub[ibin][ir] <<"rho =" << rhosub[ibin][ir]<<endl ;
            
        }                
        scalepsi[ibin]=psisub[ibin][nrbin-1];
        scaleak3psi[ibin]=ak3psisub[ibin][nrbin-1];
        cout <<"ibin = "<<ibin<<"scale =" << scalerho[ibin] <<"rho =" << scaleak3rho[ibin]<<endl ;
        
        for(int ir =0 ; ir <nrbin; ir++){ 
            psisub[ibin][ir]/=scalepsi[ibin] ;           
            rhosub[ibin][ir]/=scalerho[ibin] ;
            errpsisub[ibin][ir]/=scalepsi[ibin] ;           
            errrhosub[ibin][ir]/=scalerho[ibin] ;
            
            ak3rhosub[ibin][ir]/=scaleak3rho[ibin] ;
            ak3psisub[ibin][ir]/=scaleak3psi[ibin] ;
            errak3rhosub[ibin][ir]/=scaleak3rho[ibin] ;
            errak3psisub[ibin][ir]/=scaleak3psi[ibin] ;

            sumrho+=rhosub[ibin][ir];
            
            rhoratio[ibin][ir]=ak3rhosub[ibin][ir]/rhosub[ibin][ir];
//            rhoratio[ibin][ir]=rhosub[ibin][ir]/ak3rhosub[ibin][ir];
            errrhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errrhosub[ibin][ir]/rhosub[ibin][ir],2)+TMath::Power(errak3rhosub[ibin][ir]/ak3rhosub[ibin][ir],2));
            errrhoratio[ibin][ir]*=rhoratio[ibin][ir];

            psiratio[ibin][ir]=ak3psisub[ibin][ir]/psisub[ibin][ir];
//            psiratio[ibin][ir]=psisub[ibin][ir]/ak3psisub[ibin][ir];
            errpsiratio[ibin][ir]=TMath::Sqrt(TMath::Power(errpsisub[ibin][ir]/psisub[ibin][ir],2)+TMath::Power(errak3psisub[ibin][ir]/ak3psisub[ibin][ir],2));
            errpsiratio[ibin][ir]*=psiratio[ibin][ir];

          //  cout <<"ir = "<<ir<<"psi =" << psiratio[ibin][ir] <<"rho =" << rhoratio[ibin][ir]<<endl ;
            cout <<"ibin = " <<ibin <<"ir = "<<ir<< "ratio =" <<TMath::Abs(rhoratio[ibin][ir]-1)*100<<endl;


        }
        cout <<"sum rho =" <<sumrho<<endl ;
    }
        
    TH1F * dummy = new TH1F("dummy", "dummy", 100, 0., 1.);
    dummy->SetTitle("") ;
    dummy->SetStats(kFALSE) ;
    dummy->SetAxisRange(0., conesize, "X") ;
    //    dummy->SetAxisRange(0., 1., "X") ;
    dummy->GetXaxis()->SetTitle("radius (r)");
    fixedFontHist(dummy,1.5, 2.0);
    
    switch (opt) {
        case kDiff:
        //    if(nbin>1)c1->Divide(nbin/2, 2,0, 0);
            //  else c1->Divide(1, 2,0, 0);           
            for(int ipad = 1 ; ipad <=nbin; ipad++){
                c1->cd(ipad);
                c1->cd(ipad)->SetLogy();
                dummy->GetYaxis()->SetTitle("#rho (r)");
                //if(current==4)dummy->SetAxisRange(5.e-3, 20., "Y") ;
                // else 
                dummy->SetAxisRange(1.e-4, 50., "Y") ;         
                dummy->DrawCopy();
                //   ppRadiusRho->Draw("same P") ;
                
                RadiusRho[nbin-ipad] = new TGraphErrors(nfile, rad, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                RadiusRho[nbin-ipad]->SetMarkerStyle(20);
                RadiusRho[nbin-ipad]->SetMarkerColor(1);
                RadiusRho[nbin-ipad]->SetLineColor(1);
                RadiusRho[nbin-ipad]->Draw("same PE") ;
                
                ak3Rho[nbin-ipad] = new TGraphErrors(nfile, rad, ak3rho[nbin-ipad], errR, errak3rho[nbin-ipad]);
                //ak3Rho[nbin-ipad]->Print();
                ak3Rho[nbin-ipad]->SetMarkerStyle(24);
                ak3Rho[nbin-ipad]->SetMarkerColor(2);
                ak3Rho[nbin-ipad]->SetLineColor(2);                   
                ak3Rho[nbin-ipad]->Draw("same PE") ;

                bkgRadiusRho[nbin-ipad] = new TGraphErrors(nrbin, rad, bkgrho[nbin-ipad], errR, errbkgrho[nbin-ipad]);
                bkgRadiusRho[nbin-ipad]->SetMarkerStyle(22);
                bkgRadiusRho[nbin-ipad]->SetMarkerColor(1);
                bkgRadiusRho[nbin-ipad]->SetLineColor(1);
                bkgRadiusRho[nbin-ipad]->Draw("same PE") ;
                
                bkgak3Rho[nbin-ipad] = new TGraphErrors(nrbin, rad, ak3bkgrho[nbin-ipad], errR, errak3bkgrho[nbin-ipad]);
                //ak3Rho[nbin-ipad]->Print();
                bkgak3Rho[nbin-ipad]->SetMarkerStyle(26);
                bkgak3Rho[nbin-ipad]->SetMarkerColor(2);
                bkgak3Rho[nbin-ipad]->SetLineColor(2);
                bkgak3Rho[nbin-ipad]->Draw("same PE") ;

                if(ipad==1)drawCMSmc(0.45,0.9,type.Data());
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[nbin-ipad+1]),0.2,0.92, 15);
                if(nbin>1){
                    if(ipad==2){
                        t1->AddEntry(RadiusRho[nbin-ipad],"akPu3PF", "PL");
                        t1->AddEntry(ak3Rho[nbin-ipad],"ak3PF", "PL");  
                        t1->AddEntry(bkgRadiusRho[nbin-ipad],"akPu3PF bkg", "PL");  
                        t1->AddEntry(bkgak3Rho[nbin-ipad],"ak3PF bkg", "PL");
                        t1->Draw("same");
                    }
                    //       if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.45,0.85,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.78,17);
                    if(ipad==3)drawText(Form("jet (anti-k_{T},R=0.3)"),0.25,0.85,16);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.78,17);
                }
                else {
                    t1->AddEntry(RadiusRho[nbin-1],"akPu3PF.", "PL");
                    t1->AddEntry(ak3Rho[nbin-1],"akPu3PF", "PL");  
                    t1->AddEntry(bkgRadiusRho[nbin-1],"akPu3PF bkg", "PL");  
                    t1->AddEntry(bkgak3Rho[nbin-1],"akPu3PF bkg", "PL");

                    drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.65,0.75,17);
                    drawText(Form("|#eta|_{jet} < %.f", etacut),0.65,0.71,17);
                    drawText(Form("jet (anti-k_{T},R=0.3)"),0.65,0.25,17);
                    drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.65,0.3,17);
                }
                
            }
            
            if(SavePlot)c1->Print(Form("%s/%sDiffJetAlgoPt%.f_%.fTrkCut%.fDiffJSwithBkg_%d.gif",plotsdir.Data(),type.Data(), pt[current],pt[current+1],trackcut, date));               
            c1->Update();   
            break ;
        case kInt:
           // if(nbin>1)c1->Divide(nbin/2, 2,0, 0);
            //   else c1->Divide(1, 2,0, 0);           
            for(int ipad = 1 ; ipad <=nbin; ipad++){
                c1->cd(ipad);
                c1->cd(ipad)->SetLogy();
                dummy->GetYaxis()->SetTitle("#psi (r)");
                //if(current==4)dummy->SetAxisRange(5.e-3, 20., "Y") ;
                // else 
                dummy->SetAxisRange(1.e-5, 1.5, "Y") ;
                dummy->DrawCopy();
                //   ppRadiusRho->Draw("same P") ;
                
                RadiusPsi[nbin-ipad] = new TGraphErrors(nfile, rad, psi[nbin-ipad], errR, errpsi[nbin-ipad]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                RadiusPsi[nbin-ipad]->SetMarkerStyle(20);
                RadiusPsi[nbin-ipad]->SetMarkerColor(1);
                RadiusPsi[nbin-ipad]->SetLineColor(1);
                RadiusPsi[nbin-ipad]->Draw("same PE") ;
                
                ak3Psi[nbin-ipad] = new TGraphErrors(nfile, rad, ak3Psi[nbin-ipad], errR, errak3Psi[nbin-ipad]);
                ak3Psi[nbin-ipad]->SetMarkerStyle(22);
                ak3Psi[nbin-ipad]->SetMarkerColor(2);
                ak3Psi[nbin-ipad]->SetLineColor(2);
                ak3Psi[nbin-ipad]->Draw("same PE") ;
                //   ak3Psi[nbin-ipad]->Print();
                
                bkgRadiusPsi[nbin-ipad] = new TGraphErrors(nrbin, rad, bkgpsi[nbin-ipad], errR, errbkgpsi[nbin-ipad]);
                bkgRadiusPsi[nbin-ipad]->SetMarkerStyle(24);
                bkgRadiusPsi[nbin-ipad]->SetMarkerColor(1);
                bkgRadiusPsi[nbin-ipad]->SetLineColor(1);
                bkgRadiusPsi[nbin-ipad]->Draw("same PE") ;
                
                bkgak3Psi[nbin-ipad] = new TGraphErrors(nrbin, rad, ak3bkgpsi[nbin-ipad], errR, errak3bkgpsi[nbin-ipad]);
                //ak3Rho[nbin-ipad]->Print();
                bkgak3Psi[nbin-ipad]->SetMarkerStyle(26);
                bkgak3Psi[nbin-ipad]->SetMarkerColor(2);
                bkgak3Psi[nbin-ipad]->SetLineColor(2);
                bkgak3Psi[nbin-ipad]->Draw("same PE") ;
                
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[ipad]),0.6,0.65, 15);
                if(ipad==1)drawText("Ak PF, R=0.3",0.6,0.68,15);
                if(nbin>1){
                    if(ipad==2)drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.45,0.75,15);
                    if(ipad==4)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.45,0.75,15);
                    
                }
                else {
                    drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.65,0.78,15);
                    drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.65,0.82,15);
                }
                
            }
            
            t1->AddEntry(RadiusPsi[nbin-1],"akPu3PF RAW", "P");
            t1->AddEntry(bkgRadiusPsi[nbin-1],"akPu3PF Bkg", "P");  
            t1->AddEntry(ak3Psi[nbin-1],"ak3PF RAW", "P");  
            t1->AddEntry(bkgak3Psi[nbin-1],"ak3PF Bkg", "P");
            t1->Draw("same");
            if(SavePlot)c1->Print(Form("%s/%sDiffJetAlgoPt%.f_%.fTrkCut%.fIntJSwith%sBkg_%d.eps",plotsdir.Data(),type.Data(), pt[current],pt[current+1],trackcut, bkg.Data(), date));               
            c1->Update();   
            
            break ;
        case kDiffSub:
         //   if(nbin>1)c1->Divide(nbin/2, 2,0, 0);
            //  else c1->Divide(1, 2,0, 0);           
            for(int ipad = 1 ; ipad <=nbin; ipad++){
                c1->cd(ipad);
                c1->cd(ipad)->SetLogy();
                dummy->GetYaxis()->SetTitle("#rho (r)");
                if(trackcut>2.)dummy->SetAxisRange(8.e-3, 50., "Y") ;
                else 
                    dummy->SetAxisRange(5.e-2, 50., "Y") ;
                dummy->DrawCopy();
                //   ppRadiusRho->Draw("same P") ;
                
                
                RhoNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, rhosub[nbin-ipad], errR, errrhosub[nbin-ipad]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                RhoNoBkg[nbin-ipad]->SetMarkerStyle(24);
                RhoNoBkg[nbin-ipad]->SetMarkerColor(1);
                RhoNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                RhoNoBkg[nbin-ipad]->SetLineColor(1);
                RhoNoBkg[nbin-ipad]->Draw("same PE") ;
                
                ak3RhoNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, ak3rhosub[nbin-ipad], errR, errak3rhosub[nbin-ipad]);
                //ak3Rho[nbin-ipad]->Print();
                ak3RhoNoBkg[nbin-ipad]->SetMarkerStyle(29);
                ak3RhoNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                ak3RhoNoBkg[nbin-ipad]->SetMarkerColor(2);
                ak3RhoNoBkg[nbin-ipad]->SetLineColor(2);                   
                ak3RhoNoBkg[nbin-ipad]->Draw("same PE") ;

            //    if(ipad==1)drawCMS(0.45,0.9,pbpbLumi);
                if(ipad==1)drawCMSpp(0.45,0.9,ppLumi);
                if(nbin>1){
                    if(ipad==2){
                  //      t1->AddEntry(ak3RhoNoBkg[nbin-ipad],"2013 pp", "PL");
                        t1->AddEntry(ak3RhoNoBkg[nbin-ipad],"smear:YaxianPara", "PL");
                        t1->AddEntry(RhoNoBkg[nbin-ipad],"smear:PawanPara", "PL");
                 //       t1->AddEntry(RhoNoBkg[nbin-ipad],"2011 pp", "PL");
//                        t1->AddEntry(RhoNoBkg[nbin-ipad],"QM12", "PL");
//                        t1->AddEntry(ak3RhoNoBkg[nbin-ipad],"new smearF", "PL");              
                        t1->Draw("same");
                    }
                             if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.8,17);
                  //  if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.25,0.8,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.7,17);
                    if(ipad==3)drawText("Ak PF, R=0.3",0.55,0.8,17);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.7,17);
                }
                else {
                    t1->AddEntry(RhoNoBkg[nbin-1],"ak3PF", "PL");
                    t1->AddEntry(ak3RhoNoBkg[nbin-1],"akPu3PF", "PL");
//                    t1->AddEntry(RhoNoBkg[nbin-1],"w/o vz shift", "PL");
//                    t1->AddEntry(ak3RhoNoBkg[nbin-1],"w/ vz shift (+0.4875)", "PL");
//                    t1->AddEntry(RhoNoBkg[nbin-1],"2011 pp HI-tracking", "PL");
//                    t1->AddEntry(ak3RhoNoBkg[nbin-1],"2013 pp HI-tracking", "PL");
                    t1->Draw("same");
                //    drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.65,0.75,17);
                    drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.65,0.75,17);
                    drawText("Ak PF, R=0.3",0.65,0.68,17);
                    drawText(Form("|#eta|_{jet} < %.f", etacut),0.65,0.62,17);
                    drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.65,0.55,17);
                }
                
            }
            
            for(int ibin = 0 ; ibin <nbin; ibin++){
                if(nbin>1)c1->cd(ibin+5);
                else c1->cd(2);
          //      dummy->SetAxisRange(0.5, 2., "Y") ;
                dummy->SetAxisRange(0.9, 1.12, "Y") ;
           //     dummy->SetAxisRange(0.8, 1.12, "Y") ;
    //            dummy->GetYaxis()->SetTitle("#rho(r)^{YaxianFit}/#rho(r)^{Pawan}");
       //         if(data==kMea)dummy->GetYaxis()->SetTitle("#rho(r)^{w/}/#rho(r)^{w/o}");
//                if(data==kMea)dummy->GetYaxis()->SetTitle("#rho(r)^{ppTracking}/#rho(r)^{HITracking}");
//                if(data==kMea)dummy->GetYaxis()->SetTitle("#rho(r)^{2013}/#rho(r)^{2011}");
                if(data==kMea)dummy->GetYaxis()->SetTitle("#rho(r)^{akPu3PF}/#rho(r)^{ak3PF}");
                else 
                    dummy->GetYaxis()->SetTitle("#rho(r)^{Yaxian}/#rho(r)^{Pawan}");
           //     dummy->GetYaxis()->SetTitle("#rho(r)^{2013}/#rho(r)^{2011}");
                dummy->DrawCopy();

                RatioRho[nbin-ibin-1] = new TGraphErrors(nrbin, rad, rhoratio[nbin-ibin-1], errR, errrhoratio[nbin-ibin-1]);
                RatioRho[nbin-ibin-1]->SetMarkerStyle(20);
                RatioRho[nbin-ibin-1]->SetMarkerColor(1);
                RatioRho[nbin-ibin-1]->SetLineColor(1);
                RatioRho[nbin-ibin-1]->SetMarkerSize(1.5);
                RatioRho[nbin-ibin-1]->Draw("same PE") ;
                RatioRho[nbin-ibin-1]->Print();
                regSun(0.,1.,0.3,1.,1, 1);
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.35,0.9, 17);
            }            
            if(SavePlot)
                c1->Print(Form("%s/%sDiffSmearParaPtThres%.f_%.fTrkCut%.fDiffJSCenBin%dAfterBkgSub.gif",plotsdir.Data(),type.Data(), pt[current],pt[current+1],trackcut, nbin));
             //   c1->Print(Form("%s/%sRebinJetPtThres%.fTrkCut%.fDiffJSAfterBkgSub.gif",plotsdir.Data(),type.Data(), pt[current],trackcut));               
            c1->Update();   
            break ;
        case kIntSub:
  
            for(int ipad = 1 ; ipad <=nbin; ipad++){
                c1->cd(ipad);
                //   c1->cd(ipad)->SetLogy();
                dummy->GetYaxis()->SetTitle("#psi (r)");
                //if(current==4)dummy->SetAxisRange(5.e-3, 20., "Y") ;
                // else 
                dummy->SetAxisRange(0.3, 1.3, "Y") ;
                dummy->DrawCopy();
                //   ppRadiusRho->Draw("same P") ;
                
                PsiNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, psisub[nbin-ipad], errR, errpsisub[nbin-ipad]);
                PsiNoBkg[nbin-ipad]->SetMarkerStyle(20);
                PsiNoBkg[nbin-ipad]->SetMarkerColor(1);
                PsiNoBkg[nbin-ipad]->SetLineColor(1);
                PsiNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                PsiNoBkg[nbin-ipad]->Draw("same PE") ;
                
                ak3PsiNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, ak3psisub[nbin-ipad], errR, errak3psisub[nbin-ipad]);
                ak3PsiNoBkg[nbin-ipad]->SetMarkerStyle(25);
                ak3PsiNoBkg[nbin-ipad]->SetMarkerColor(1);
                ak3PsiNoBkg[nbin-ipad]->SetLineColor(1);
                ak3PsiNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                ak3PsiNoBkg[nbin-ipad]->Draw("same PE") ;
                if(ipad==1)drawCMSmc(0.5,0.85,type.Data());
                if(nbin>1){
                    if(ipad==2){
                        t1->AddEntry(PsiNoBkg[nbin-1],"akPu3PF", "PL");
                        t1->AddEntry(ak3PsiNoBkg[nbin-1],"ak3PF", "PL");  
                        t1->Draw("same");
                    }
                    //         if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.25,0.85,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.8,17);
                    if(ipad==3)drawText("Ak PF, R=0.3",0.55,0.85,17);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.8,17);
                }
                else {
                    t1->AddEntry(PsiNoBkg[nbin-1],"akPu3PF", "PL");
                    t1->AddEntry(ak3PsiNoBkg[nbin-1],"ak3PF", "PL");  
                    t1->Draw("same");
                    drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.6,0.75,17);
                }
                
            }
            
            for(int ibin = 0 ; ibin <nbin; ibin++){
                if(nbin>1)c1->cd(ibin+5);
                else c1->cd(2);
                dummy->SetAxisRange(0.7, 1.3, "Y") ;
                dummy->GetYaxis()->SetTitle("#psi(r)^{ak3PF}/#psi(r)^{akPu3PF}");
                dummy->DrawCopy();
                
                RatioPsi[nbin-ibin-1] = new TGraphErrors(nrbin, rad, psiratio[nbin-ibin-1], errR, errpsiratio[nbin-ibin-1]);
                RatioPsi[nbin-ibin-1]->SetMarkerStyle(20);
                RatioPsi[nbin-ibin-1]->SetMarkerColor(1);
                RatioPsi[nbin-ibin-1]->SetLineColor(1);
                RatioPsi[nbin-ibin-1]->SetMarkerSize(1.5);
                RatioPsi[nbin-ibin-1]->Draw("same PE") ;
                regSun(0.,1.,0.3,1.,1, 1);
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.25,0.9, 17);
            }            
            if(SavePlot)c1->Print(Form("%s/%sDiffJetAlgoPt%.fTrkCut%.fqRecoIntJSAfter%sBkgSub.pdf",plotsdir.Data(),type.Data(), pt[current],trackcut, bkg.Data()));               
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
void drawCMSmc(float px, float py, TString coll) {
    TLatex *cms = new TLatex(px,py,"CMS Preliminary");
    cms->SetTextFont(63);
    cms->SetTextSize(17);
    cms->SetNDC();
    cms->Draw();
    if(coll=="HI")TLatex *lumi = new TLatex(px,py-0.05,"PYTHIA+HYDJET");
    else TLatex *lumi = new TLatex(px,py-0.05,"PYTHIA");
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
    TLatex *lumi = new TLatex(px,py-0.05,Form("#intL dt = %.1f pb^{-1}",nLumi));
//    TLatex *lumi = new TLatex(px,py-0.10,Form("#intL dt = %.1f nb^{-1}",nLumi));
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
    h->SetTitleSize(16);
    h->SetLabelSize(18,"Y");
    h->SetLabelSize(18,"X");
    h->SetTitleFont(43,"Y");
    h->SetTitleSize(22,"Y");
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
