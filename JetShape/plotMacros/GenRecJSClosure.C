

/*
 *  Compare the JS with the reconsturcted leading corresponding to the generator level jet shape
 * 
 *
 *  Crea
 ted by ymao on 12/28/11.
 
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
//const double pt[]={100., 120.,300., 500.};
//const double pt[]={100., 120., 140., 160., 200., 300., 500.};
const double pt[]={100., 300.};

const int nrbin = 6 ;
double rbin[nrbin+1]={0.0,0.05,0.10, 0.15, 0.20,0.25,0.30};

const Double_t kw = 1000;
const Double_t kh = 560;

const Double_t pbpbLumi = 150 ;
const Double_t ppLumi = 231 ;

enum Option_t {kDiff, kInt, kDiffSub, kIntSub, kBkgRatio} ; 
enum Display_t {kPP, kHI} ;
enum Method_t {kRC, kER} ;
enum Order_t {kRecGen, kGenRec} ; //starting from reconstruction or reverse
//---------------------------------------------------------------------

void GenRecJSClosure(Option_t opt=kDiff,  Method_t met=kER, Display_t data = kPP, Order_t ord = kRecGen, int current=0){
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
    if(data==kPP) type="PP";
    else type = "HI" ;
    TString bkg ;
    if(met==kER) bkg = "ER";
    else bkg = "RC";
    TString effTab = "HistIterTrkCorrtest"; // "HistIterTrkCorrtest" ; // "HistIterTrkCorrv14fXSec";
//    if(data==kPP)
//        const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/GenJet" ;      
//    else 
//        const char * kHomeDir = "/FullStatPbPb/RebinJetPt/TrkCorr4CentBin" ;      
    const char * kHomeDir = "/Users/ymao/group/CMS/anaOutputs/JetShape/CorrectdR/NoMatchGen/NoTrkAlgoCut" ;      

    TString fileName ;
    TFile * f ;
    
    if(data==kPP) {
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
    
    TH2F * GenIntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * GenbkgIntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * GenDiffJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    TH2F * GenbkgDiffJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    
    TH1F * GenIntegrated[nbin][nrbin] ;
    TH1F * GenbkgIntegrated[nbin][nrbin] ;
    
    TH1F * Gendiff[nbin][nrbin] ;
    TH1F * Genbkgdiff[nbin][nrbin] ;
    
    TGraphErrors * RhoGen[nbin];
    TGraphErrors * PsiGen[nbin] ;
    TGraphErrors * bkgRhoGen[nbin];
    TGraphErrors * bkgPsiGen[nbin] ;
    
    TGraphErrors * RhoNoBkg[nbin];
    TGraphErrors * PsiNoBkg[nbin] ;
    TGraphErrors * GenRhoNoBkg[nbin];
    TGraphErrors * GenPsiNoBkg[nbin] ;

    TGraphErrors * RatioRho[nbin];
    TGraphErrors * RatioPsi[nbin] ;
    TGraph * bkgRhoRatio[nbin];
    TGraph * bkgPsiRatio[nbin] ;

    TH1::SetDefaultSumw2();
    
    double rho[nbin][nrbin];
    double errrho[nbin][nrbin];
    double psi[nbin][nrbin];
    double errpsi[nbin][nrbin];
    
    double genrho[nbin][nrbin];
    double errgenrho[nbin][nrbin];
    double genpsi[nbin][nrbin];
    double errgenpsi[nbin][nrbin];
    
    double bkgrho[nbin][nrbin];
    double errbkgrho[nbin][nrbin];
    double bkgpsi[nbin][nrbin];
    double errbkgpsi[nbin][nrbin];
    
    double genbkgrho[nbin][nrbin];
    double errgenbkgrho[nbin][nrbin];
    double genbkgpsi[nbin][nrbin];
    double errgenbkgpsi[nbin][nrbin];
    
    double rhosub[nbin][nrbin];
    double errrhosub[nbin][nrbin];
    double psisub[nbin][nrbin];
    double errpsisub[nbin][nrbin];
    
    double genrhosub[nbin][nrbin];
    double errgenrhosub[nbin][nrbin];
    double genpsisub[nbin][nrbin];
    double errgenpsisub[nbin][nrbin];
    
    double rhoratio[nbin][nrbin];
    double errrhoratio[nbin][nrbin];
    double psiratio[nbin][nrbin];
    double errpsiratio[nbin][nrbin];

    double doublerhoratio[nbin][nrbin];
    double errdoublerhoratio[nbin][nrbin];
    double doublepsiratio[nbin][nrbin];
    double errdoublepsiratio[nbin][nrbin];

    double scalepsi[nbin];
    double scalegenpsi[nbin];
    double scalerho[nbin];
    double scalegenrho[nbin];
    
    TString canv_name = "c1";
    //    canv_name += parti;
    //    canv_name += cent;
    //    makeMultiPanelCanvas(c1,4,1,0.0,0.0,0.2,0.2,0.02);
    if(nbin>1){
        c1 = new TCanvas(canv_name," ",10,10,kw,kh);
        makeMultiPanelCanvas(c1,4,2,0.0,0.0,0.2,0.2,0.02);
    }
    else {
        c1 = new TCanvas(canv_name," ",800, 800);
        c1->Divide(1, 2, 0, 0);
    }
    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);
    
    
    if(nbin>1)TLegend *t1=new TLegend(0.25,0.7,0.55,0.93);
    else TLegend *t1=new TLegend(0.55,0.75,0.75,0.93);
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
        case kPP:
            switch(ord){
                case kRecGen:
                    //   fileName = Form("MCPP_ChargedJetPt80_40Trk4Eta2_BkgEtaBound3_CenBin1_Nrbin6_Pthat80mit_ivan_pp276Dijet80_merged.root");
                    switch(met){
                        case kER:
                   //         fileName = Form("merged_MCPP_Ak3PFIncJetPt100_HistCorrTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_mergedFile.root", trackcut, nbin);
                            fileName = Form("mergedCSdiff_MCPP2013_Ak3PFIncJetPt100_2013%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin1_Nrbin6_LJbin1_trkbin3_pp276Dijet_ppHiIterativeTrack_merged.root", effTab.Data(), trackcut);           
                 //           fileName = Form("mergedCSdiff_MCPP_Ak3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin%d_Nrbin6_LJbin1_trkbin3_pp276Dijet_merged.root", effTab.Data(), trackcut, nbin);
                            break ;
                        case kRC:
                 //   fileName = Form("merged_MCPP_ChargedJetPt80_40Trk%.fEta2_RCbkgEtaBound20_CenBin1_Nrbin6_pp276Dijet.root", trackcut);
                    fileName = Form("MCPP_PFRecAxisJetPt100_50Trk%.fEta2_RCbkgJS_CenBin1_Nrbin6_LJbin6_Pthat80pp276Dijet80_merged.root", trackcut);
                            break ;
                    }
                    break ;
                case kGenRec:
                    switch(met){
                        case kER:
                            fileName = Form("MCPP_AbsGenChargedJetPt80_40Trk%.fEta2_ERbkg_CenBin1_Nrbin6_Pthat80Pythia80_signal01_HiForest2_v19_0.root", trackcut);                 
                            break ;
                        case kRC:
                            fileName = Form("MCPP_AbsGenChargedJetPt80_40Trk%.fEta2_RCbkg_CenBin1_Nrbin6_Pthat80Pythia80_signal01_HiForest2_v19_0.root", trackcut);
                            break ;
                    }
            }
            break ;
        case kHI:
            switch(ord){
                case kRecGen:
                    switch(met){
                case kER:
                 //   fileName = Form("MCHI_PFIncJetPt100_HistCorrTrk%.fEtaLimit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1_Dijet100_HydjetDrum_v27_mergedV1.root", trackcut);
                 //   fileName = Form("MCHI_PFIncJetPt100_HistCorrTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1_pyquenFull100_HYDJET.root", trackcut);
                  //      fileName = Form("merged_MCHI_AkPu3PFIncJetPt100_HistCent4BinCorrTrk%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1_Dijet_HydjetDrum_v27v28.root", trackcut);
                            fileName = Form("mergedCSdiff_MCHI_AkPu3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_Dijet_HydjetDrum_v27v28.root", effTab.Data(), trackcut);
                  //  fileName = Form("mergedCSdiff_MCHI_AkPu3PFIncJetPt100_%s%.fEtaCut20Limit3_ERbkgJSNormJetCone3_CenBin4_Nrbin6_LJbin1trkbin3_Dijet_HydjetDrum_v27v28.root", effTab.Data(), trackcut);
                    break ;
                case kRC:
          //  fileName = Form("MCHI_PFJetPt100_50CorrTrkR%.fEta2_RCbkgJS_CenBin4_Nrbin6_LJbin6_Pthat80Pythia80_HydjetDrum_mix01_HiForest2_v20.root", trackcut);
          fileName=Form("MCHI_PFJetPt100_50Trk%.fEta2_RCbkgJS_CenBin4_Nrbin6_LJbin6_Pthat80Pythia80_HydjetDrum_mix01_HiForest2_v20.root", trackcut);
                            
                    break ;
                    }
                    break ;
                case kGenRec:
                    switch(met){
                        case kER:
                            fileName = Form("MCHI_AbsGenChargedJetPt80_40Trk%.fEta2_ERbkgEtaBound3_CenBin6_Nrbin6_Pthat80Pythia80_HydjetDrum_mix01_HiForest2_v20.root", trackcut);
                            break ;
                        case kRC:
                            fileName = Form("MCHI_AbsGenChargedJetPt80_40Trk%.fEta2_RCbkgEtaBound20_CenBin6_Nrbin6_Pthat80Pythia80_HydjetDrum_mix01_HiForest2_v20.root", trackcut);
                            break ;
                    }
                    break ;
            }

    }
    
    if(data==kPP) f = TFile::Open(Form("%s/MyselfTrkEff/%s", kHomeDir, fileName.Data()), "readonly");
        else f = TFile::Open(Form("%s/%s", kHomeDir, fileName.Data()), "readonly");
    for(int ibin = 0 ; ibin <nbin; ibin++){   
        scalegenpsi[ibin]=0. ;
        scalepsi[ibin]=0. ;
        scalegenrho[ibin]=0. ;
        scalerho[ibin]=0. ;
        double sumrho = 0.0 ;
        for(int ir =0 ; ir <nrbin; ir++){  
            rho[ibin][ir] = 0 ;
            errrho[ibin][ir] = 0 ;
            psi[ibin][ir]=0 ;
            errpsi[ibin][ir]=0 ;  
            genrho[ibin][ir] = 0 ;
            errgenrho[ibin][ir] = 0 ;
            genpsi[ibin][ir]=0 ;
            errgenpsi[ibin][ir]=0 ;  
            bkgrho[ibin][ir] = 0 ;
            errbkgrho[ibin][ir] = 0 ;
            bkgpsi[ibin][ir]=0 ;
            errbkgpsi[ibin][ir]=0 ;  
            genbkgrho[ibin][ir] = 0 ;
            errgenbkgrho[ibin][ir] = 0 ;
            genbkgpsi[ibin][ir]=0 ;
            errgenbkgpsi[ibin][ir]=0 ;  
            rhosub[ibin][ir] = 0 ;
            errrhosub[ibin][ir] = 0 ;
            psisub[ibin][ir]=0 ;
            errpsisub[ibin][ir]=0 ;  
            genrhosub[ibin][ir] = 0 ;
            errgenrhosub[ibin][ir] = 0 ;
            genpsisub[ibin][ir]=0 ;
            errgenpsisub[ibin][ir]=0 ;  
            
            rhoratio[ibin][ir] = 0 ;
            errrhoratio[ibin][ir] = 0 ;
            psiratio[ibin][ir]=0 ;
            errpsiratio[ibin][ir]=0 ;  
            
            GenIntJS[ibin][ir]=(TH2F*)f->Get(Form("genIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            IntJS[ibin][ir]=(TH2F*)f->Get(Form("IntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            GenbkgIntJS[ibin][ir]=(TH2F*)f->Get(Form("genbkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            bkgIntJS[ibin][ir]=(TH2F*)f->Get(Form("bkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            
            GenDiffJS[ibin][ir]=(TH2F*)f->Get(Form("gendifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            DiffJS[ibin][ir]=(TH2F*)f->Get(Form("differentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            GenbkgDiffJS[ibin][ir]=(TH2F*)f->Get(Form("genbkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            bkgDiffJS[ibin][ir]=(TH2F*)f->Get(Form("bkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]));  
            
            GenIntJS[ibin][ir]->Sumw2();
            GenIntegrated[ibin][ir]=(TH1F*)GenIntJS[ibin][ir]->ProjectionY(Form("GenJetPt%.f_%.fPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),GenIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), GenIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            
            IntJS[ibin][ir]->Sumw2();
            Integrated[ibin][ir]=(TH1F*)IntJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),IntJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), IntJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            
            GenbkgIntJS[ibin][ir]->Sumw2();
            GenbkgIntegrated[ibin][ir]=(TH1F*)GenbkgIntJS[ibin][ir]->ProjectionY(Form("GenJetPt%.f_%.fbkgPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),GenbkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), GenbkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            
            bkgIntJS[ibin][ir]->Sumw2();
            bkgIntegrated[ibin][ir]=(TH1F*)bkgIntJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fbkgPsidR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),bkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), bkgIntJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            
            GenDiffJS[ibin][ir]->Sumw2();
            Gendiff[ibin][ir]=(TH1F*)GenDiffJS[ibin][ir]->ProjectionY(Form("GenJetPt%.f_%.fRhodR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),GenDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), GenDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            
            DiffJS[ibin][ir]->Sumw2();
            differential[ibin][ir]=(TH1F*)DiffJS[ibin][ir]->ProjectionY(Form("JetPt%.f_%.fRhodR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),DiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), DiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            
            GenbkgDiffJS[ibin][ir]->Sumw2();
            Genbkgdiff[ibin][ir]=(TH1F*)GenbkgDiffJS[ibin][ir]->ProjectionY(Form("GenJetPt%.f_%.fbkgRhodR%.f_%.f_Cen%.f-%.f%%",pt[current], pt[current+1],rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]),GenbkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current]), GenbkgDiffJS[ibin][ir]->GetXaxis()->FindBin(pt[current+1])-1, "");
            
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
            
            genpsi[ibin][ir]=GenIntegrated[ibin][ir]->GetMean(1);
            errgenpsi[ibin][ir]=GenIntegrated[ibin][ir]->GetMeanError(1);
            genbkgpsi[ibin][ir]=GenbkgIntegrated[ibin][ir]->GetMean(1);
            errgenbkgpsi[ibin][ir]=GenbkgIntegrated[ibin][ir]->GetMeanError(1);
            genrho[ibin][ir]=Gendiff[ibin][ir]->GetMean(1);
            errgenrho[ibin][ir]=Gendiff[ibin][ir]->GetMeanError(1);
            genbkgrho[ibin][ir]=Genbkgdiff[ibin][ir]->GetMean(1);
            errgenbkgrho[ibin][ir]=Genbkgdiff[ibin][ir]->GetMeanError(1);
            
            //            for(int i = 0 ; i<Integrated[ibin][ir]->GetNbinsX(); i++){
            //                psi[ibin][ir]+= Integrated[ibin][ir]->GetBinContent(i)*Integrated[ibin][ir]->GetBinCenter(i);
            //                errpsi[ibin][ir]+= Integrated[ibin][ir]->GetBinError(i)*Integrated[ibin][ir]->GetBinCenter(i); 
            //            }
            //            if(Integrated[ibin][ir]->Integral())psi[ibin][ir]/=Integrated[ibin][ir]->Integral();
            //            if(Integrated[ibin][ir]->Integral())errpsi[ibin][ir]/=Integrated[ibin][ir]->Integral();  
            //            
            //            for(int i = 0 ; i<GenIntegrated[ibin][ir]->GetNbinsX(); i++){
            //                genpsi[ibin][ir]+= GenIntegrated[ibin][ir]->GetBinContent(i)*GenIntegrated[ibin][ir]->GetBinCenter(i);
            //                errgenpsi[ibin][ir]+= GenIntegrated[ibin][ir]->GetBinError(i)*GenIntegrated[ibin][ir]->GetBinCenter(i); 
            //            }
            //            if(GenIntegrated[ibin][ir]->Integral())genpsi[ibin][ir]/=GenIntegrated[ibin][ir]->Integral();
            //            if(GenIntegrated[ibin][ir]->Integral())errgenpsi[ibin][ir]/=GenIntegrated[ibin][ir]->Integral();  
            //            
            //            for(int i = 0 ; i<bkgIntegrated[ibin][ir]->GetNbinsX(); i++){
            //                bkgpsi[ibin][ir]+= bkgIntegrated[ibin][ir]->GetBinContent(i)*bkgIntegrated[ibin][ir]->GetBinCenter(i);
            //                errbkgpsi[ibin][ir]+= bkgIntegrated[ibin][ir]->GetBinError(i)*bkgIntegrated[ibin][ir]->GetBinCenter(i); 
            //            }
            //            if(bkgIntegrated[ibin][ir]->Integral())bkgpsi[ibin][ir]/=bkgIntegrated[ibin][ir]->Integral();
            //            if(bkgIntegrated[ibin][ir]->Integral())errbkgpsi[ibin][ir]/=bkgIntegrated[ibin][ir]->Integral();  
            //            
            //            for(int i = 0 ; i<GenbkgIntegrated[ibin][ir]->GetNbinsX(); i++){
            //                genbkgpsi[ibin][ir]+= GenbkgIntegrated[ibin][ir]->GetBinContent(i)*GenbkgIntegrated[ibin][ir]->GetBinCenter(i);
            //                errgenbkgpsi[ibin][ir]+= GenbkgIntegrated[ibin][ir]->GetBinError(i)*GenbkgIntegrated[ibin][ir]->GetBinCenter(i); 
            //            }
            //            if(GenbkgIntegrated[ibin][ir]->Integral())genbkgpsi[ibin][ir]/=GenbkgIntegrated[ibin][ir]->Integral();
            //            if(GenbkgIntegrated[ibin][ir]->Integral())errgenbkgpsi[ibin][ir]/=GenbkgIntegrated[ibin][ir]->Integral();                 
            //            
            //            for(int i = 0 ; i<differential[ibin][ir]->GetNbinsX(); i++){
            //                rho[ibin][ir]+= differential[ibin][ir]->GetBinContent(i)*differential[ibin][ir]->GetBinCenter(i);
            //                errrho[ibin][ir]+= differential[ibin][ir]->GetBinError(i)*differential[ibin][ir]->GetBinCenter(i); 
            //            }
            //            if(differential[ibin][ir]->Integral())rho[ibin][ir]/=differential[ibin][ir]->Integral();
            //            if(differential[ibin][ir]->Integral())errrho[ibin][ir]/=differential[ibin][ir]->Integral();  
            //            for(int i = 0 ; i<bkgdifferential[ibin][ir]->GetNbinsX(); i++){
            //                bkgrho[ibin][ir]+= bkgdifferential[ibin][ir]->GetBinContent(i)*bkgdifferential[ibin][ir]->GetBinCenter(i);
            //                errbkgrho[ibin][ir]+= bkgdifferential[ibin][ir]->GetBinError(i)*bkgdifferential[ibin][ir]->GetBinCenter(i); 
            //            }
            //            if(bkgdifferential[ibin][ir]->Integral())bkgrho[ibin][ir]/=bkgdifferential[ibin][ir]->Integral();
            //            if(differential[ibin][ir]->Integral())errbkgrho[ibin][ir]/=bkgdifferential[ibin][ir]->Integral();  
            //            
            //            for(int i = 0 ; i<Gendiff[ibin][ir]->GetNbinsX(); i++){
            //                genrho[ibin][ir]+= Gendiff[ibin][ir]->GetBinContent(i)*Gendiff[ibin][ir]->GetBinCenter(i);
            //                errgenrho[ibin][ir]+= Gendiff[ibin][ir]->GetBinError(i)*Gendiff[ibin][ir]->GetBinCenter(i);                         
            //            }
            //            if(Gendiff[ibin][ir]->Integral())genrho[ibin][ir]/=Gendiff[ibin][ir]->Integral();
            //            if(Gendiff[ibin][ir]->Integral())errgenrho[ibin][ir]/=Gendiff[ibin][ir]->Integral();  
            //            for(int i = 0 ; i<Genbkgdiff[ibin][ir]->GetNbinsX(); i++){
            //                genbkgrho[ibin][ir]+= Genbkgdiff[ibin][ir]->GetBinContent(i)*Genbkgdiff[ibin][ir]->GetBinCenter(i);
            //                errgenbkgrho[ibin][ir]+= Genbkgdiff[ibin][ir]->GetBinError(i)*Genbkgdiff[ibin][ir]->GetBinCenter(i);     
            //            }
            //            if(Genbkgdiff[ibin][ir]->Integral())genbkgrho[ibin][ir]/=Genbkgdiff[ibin][ir]->Integral();
            //            if(Genbkgdiff[ibin][ir]->Integral())errgenbkgrho[ibin][ir]/=Genbkgdiff[ibin][ir]->Integral();  
            
              cout <<"ir = "<<ir<<"psi =" << bkgpsi[ibin][ir] <<"rho =" << genrho[ibin][ir]<<"bkgrho =" << genbkgrho[ibin][ir]<<endl ;
            //nromalize to delta cone size to consistent for definition
            rho[ibin][ir]/=deltacone ;
            genrho[ibin][ir]/=deltacone ;
            bkgrho[ibin][ir]/=deltacone ;
            genbkgrho[ibin][ir]/=deltacone ;
            errrho[ibin][ir]/=deltacone ;
            errgenrho[ibin][ir]/=deltacone ;
            errbkgrho[ibin][ir]/=deltacone ;
            errgenbkgrho[ibin][ir]/=deltacone ;
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
            
            genrhosub[ibin][ir]=genrho[ibin][ir]-genbkgrho[ibin][ir];
            errgenrhosub[ibin][ir]=TMath::Sqrt(TMath::Power(errgenrho[ibin][ir],2)+TMath::Power(errgenbkgrho[ibin][ir],2));
            genpsisub[ibin][ir]=genpsi[ibin][ir]-genbkgpsi[ibin][ir];
            errgenpsisub[ibin][ir]=TMath::Sqrt(TMath::Power(errgenpsi[ibin][ir],2)+TMath::Power(errgenbkgpsi[ibin][ir],2));
            
            //            //correct the efficiency 
            //            rho[ibin][ir]/=bkgrhoratio[ibin][ir];
            //            errrho[ibin][ir]/=bkgrhoratio[ibin][ir];
            //            bkgrho[ibin][ir]/=bkgrhoratio[ibin][ir];
            //            errbkgrho[ibin][ir]/=bkgrhoratio[ibin][ir];
            
            //                    scalegenrho[ibin]+=genrhosub[ibin][ir]*deltacone;
            //                    scalerho[ibin]+=rhosub[ibin][ir]*deltacone;                    
            scalegenrho[ibin]+=genrhosub[ibin][ir]*deltacone;
            scalerho[ibin]+=rhosub[ibin][ir]*deltacone;                    
            
            //                    sumrho+=genrhosub[ibin][ir]*deltacone;
                                cout <<"ir = "<<ir<<"psisub =" << genpsisub[ibin][ir] <<"rhosub =" << genrhosub[ibin][ir]<<endl ;
            
        }                
        scalepsi[ibin]=psisub[ibin][nrbin-1];
        scalegenpsi[ibin]=genpsisub[ibin][nrbin-1];
        cout <<"ibin = "<<ibin<<"scale =" << scalerho[ibin] <<"rho =" << scalegenrho[ibin]<<endl ;
        
        for(int ir =0 ; ir <nrbin; ir++){ 
            psisub[ibin][ir]/=scalepsi[ibin] ;           
            rhosub[ibin][ir]/=scalerho[ibin] ;
            errpsisub[ibin][ir]/=scalepsi[ibin] ;           
            errrhosub[ibin][ir]/=scalerho[ibin] ;
            
            genrhosub[ibin][ir]/=scalegenrho[ibin] ;
            genpsisub[ibin][ir]/=scalegenpsi[ibin] ;
            errgenrhosub[ibin][ir]/=scalegenrho[ibin] ;
            errgenpsisub[ibin][ir]/=scalegenpsi[ibin] ;
            sumrho+=rhosub[ibin][ir];
            
            if(opt==kDiff){
                rhoratio[ibin][ir]=rho[ibin][ir]/genrho[ibin][ir];
                errrhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errrho[ibin][ir]/rho[ibin][ir],2)+TMath::Power(errgenrho[ibin][ir]/genrho[ibin][ir],2));
                errrhoratio[ibin][ir]*=rhoratio[ibin][ir];
                psiratio[ibin][ir]=bkgrho[ibin][ir]/genbkgrho[ibin][ir];
                errpsiratio[ibin][ir]=TMath::Sqrt(TMath::Power(errbkgrho[ibin][ir]/bkgrho[ibin][ir],2)+TMath::Power(errgenbkgrho[ibin][ir]/genbkgrho[ibin][ir],2));
                errpsiratio[ibin][ir]*=psiratio[ibin][ir] ; 
                cout <<"ibin = " <<ibin <<"ir = "<<ir<< "ratio =" <<TMath::Abs(rhoratio[ibin][ir]-1)*100<<endl;
                cout <<"ibin = " <<ibin <<"ir = "<<ir<< "ratiobkg =" <<TMath::Abs(psiratio[ibin][ir]-1)*100<<endl;
            }
            else {
            rhoratio[ibin][ir]=rhosub[ibin][ir]/genrhosub[ibin][ir];
            errrhoratio[ibin][ir]=TMath::Sqrt(TMath::Power(errrhosub[ibin][ir]/rhosub[ibin][ir],2)+TMath::Power(errgenrhosub[ibin][ir]/genrhosub[ibin][ir],2));
            errrhoratio[ibin][ir]*=rhoratio[ibin][ir];
            psiratio[ibin][ir]=psisub[ibin][ir]/genpsisub[ibin][ir];
            errpsiratio[ibin][ir]=TMath::Sqrt(TMath::Power(errpsisub[ibin][ir]/psisub[ibin][ir],2)+TMath::Power(errgenpsisub[ibin][ir]/genpsisub[ibin][ir],2));
            errpsiratio[ibin][ir]*=psiratio[ibin][ir] ;
            }
            
         //   cout <<"ir = "<<ir<<"bkg =" << psiratio[ibin][ir] <<"rho =" << rhoratio[ibin][ir]<<endl ;    

        }
        cout <<"sum rho =" <<sumrho<<endl ;
    }
    
//    if(SaveFile){
//        ofstream myfile;
//        myfile.open(Form("%ssysTrkEff.txt", type.Data()));
//        for(int ibin = 0 ; ibin <nbin; ibin++){
//            myfile << 100*TMath::Abs(rhoratio[ibin][0]-1) <<", "<< 100*TMath::Abs(rhoratio[ibin][1]-1)<<", "<< 100*TMath::Abs(rhoratio[ibin][2]-1)<<", "<< 100*TMath::Abs(rhoratio[ibin][3]-1)<<", "<< 100*TMath::Abs(rhoratio[ibin][4]-1)<<", "<< 100*TMath::Abs(rhoratio[ibin][5]-1)<<"\n";
//            myfile << " background closure:" << "\n";
//            myfile << 100*TMath::Abs(psiratio[ibin][0]-1) <<", "<< 100*TMath::Abs(psiratio[ibin][1]-1)<<", "<< 100*TMath::Abs(psiratio[ibin][2]-1)<<", "<< 100*TMath::Abs(psiratio[ibin][3]-1)<<", "<< 100*TMath::Abs(psiratio[ibin][4]-1)<<", "<< 100*TMath::Abs(psiratio[ibin][5]-1)<<"\n";
//        }
//        myfile.close();
//    }
    
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
                if(data==kPP) dummy->SetAxisRange(5.e-4, 50., "Y") ;
                    else dummy->SetAxisRange(5.e-3, 50., "Y") ;
                dummy->DrawCopy();
                //   ppRadiusRho->Draw("same P") ;
                
                RadiusRho[nbin-ipad] = new TGraphErrors(nfile, rad, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                RadiusRho[nbin-ipad]->SetMarkerStyle(24);
                RadiusRho[nbin-ipad]->SetMarkerColor(1);
                RadiusRho[nbin-ipad]->SetLineColor(1);
                RadiusRho[nbin-ipad]->SetMarkerSize(1.5);
                RadiusRho[nbin-ipad]->Draw("same PE") ;
                
                RhoGen[nbin-ipad] = new TGraphErrors(nfile, rad, genrho[nbin-ipad], errR, errgenrho[nbin-ipad]);
                //RhoGen[nbin-ipad]->Print();
                RhoGen[nbin-ipad]->SetMarkerStyle(29);
                RhoGen[nbin-ipad]->SetMarkerColor(1);
                RhoGen[nbin-ipad]->SetLineColor(1);                   
                RhoGen[nbin-ipad]->SetMarkerSize(1.5);
                RhoGen[nbin-ipad]->Draw("same PE") ;
                
                bkgRadiusRho[nbin-ipad] = new TGraphErrors(nrbin, rad, bkgrho[nbin-ipad], errR, errbkgrho[nbin-ipad]);
                bkgRadiusRho[nbin-ipad]->SetMarkerStyle(24);
                bkgRadiusRho[nbin-ipad]->SetMarkerColor(2);
                bkgRadiusRho[nbin-ipad]->SetMarkerSize(1.5);
                bkgRadiusRho[nbin-ipad]->SetLineColor(2);
                bkgRadiusRho[nbin-ipad]->Draw("same PE") ;
                
                bkgRhoGen[nbin-ipad] = new TGraphErrors(nrbin, rad, genbkgrho[nbin-ipad], errR, errgenbkgrho[nbin-ipad]);
                //RhoGen[nbin-ipad]->Print();
                bkgRhoGen[nbin-ipad]->SetMarkerStyle(29);
                bkgRhoGen[nbin-ipad]->SetMarkerColor(2);
                bkgRhoGen[nbin-ipad]->SetLineColor(2);
                bkgRhoGen[nbin-ipad]->SetMarkerSize(1.5);
                bkgRhoGen[nbin-ipad]->Draw("same PE") ;
                
              //  if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[ipad]),0.8,0.9, 17);
                if(nbin>1){
                    if(ipad==1)drawCMSmc(0.5,0.85,type.Data());
                    if(ipad==2)drawText("Ak PF, R=0.3",0.3,0.85,17);
                    if(ipad==2)drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.75,17);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.35,0.85,17);
                    
                }
                else {
                    drawCMSmc(0.75,0.9,type.Data());
                    drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.25,0.45,17);
                    if(ipad==1)drawText("Ak PF, R=0.3",0.25,0.55,17);
                    drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.35,17);
                }
                
            }
            
            t1->AddEntry(RadiusRho[nbin-1],"Reco: (S+B)", "P");
            t1->AddEntry(bkgRadiusRho[nbin-1],"Reco: B", "P");  
            t1->AddEntry(RhoGen[nbin-1],"Gen: (S+B)", "P");  
            t1->AddEntry(bkgRhoGen[nbin-1],"Gen: B", "P");
            t1->Draw("same");
            
            for(int ibin = 0 ; ibin <nbin; ibin++){
                if(nbin>1)c1->cd(ibin+5);
                else c1->cd(2);
                dummy->SetAxisRange(0.8, 1.22, "Y") ;
                dummy->GetYaxis()->SetTitle("#rho(r)^{corrected}/#rho(r)^{gen}");
                dummy->DrawCopy();
                
                RatioRho[nbin-ibin-1] = new TGraphErrors(nrbin, rad, rhoratio[nbin-ibin-1], errR, errrhoratio[nbin-ibin-1]);
                RatioRho[nbin-ibin-1]->SetMarkerStyle(20);
                RatioRho[nbin-ibin-1]->SetMarkerColor(1);
                RatioRho[nbin-ibin-1]->SetLineColor(1);
                RatioRho[nbin-ibin-1]->SetMarkerSize(1.5);
                RatioRho[nbin-ibin-1]->Draw("same PE") ;

                RatioPsi[nbin-ibin-1] = new TGraphErrors(nrbin, rad, psiratio[nbin-ibin-1], errR, errpsiratio[nbin-ibin-1]);
                RatioPsi[nbin-ibin-1]->SetMarkerStyle(29);
                RatioPsi[nbin-ibin-1]->SetMarkerColor(2);
                RatioPsi[nbin-ibin-1]->SetLineColor(2);
                RatioPsi[nbin-ibin-1]->SetMarkerSize(1.5);
                RatioPsi[nbin-ibin-1]->Draw("same PE") ;

                regSun(0.,1.,0.3,1.,1, 1);
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.25,0.9, 17);
            }            
            if(SavePlot)c1->Print(Form("%s/%sJetPt%.f_%.fTrkCut%.f%sDiffJSwith%sBkg.pdf",plotsdir.Data(),type.Data(), pt[current],pt[current+1],trackcut, effTab.Data(), bkg.Data())); 
            c1->Update();   
            if(SaveFile){
                TFile * outf = new TFile(Form("%sIncJet%sDiffJSRatiowith%sbkg.root", type.Data(), effTab.Data(), bkg.Data()), "CREATE");
                for(int ibin = 0 ; ibin <nbin; ibin++){                    
                    RatioRho[ibin]->Write(Form("JetConeRhoERbkg_Cen%d-%d", centr[ibin],centr[ibin+1]));
                    RatioPsi[ibin]->Write(Form("BkgConeRhoERbkg_Cen%d-%d", centr[ibin],centr[ibin+1]));
                }
                outf->Close();  
            }

            break ;
        case kInt:
            if(nbin>1)c1->Divide(nbin/2, 2,0, 0);
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
                
                PsiGen[nbin-ipad] = new TGraphErrors(nfile, rad, genpsi[nbin-ipad], errR, errgenpsi[nbin-ipad]);
                PsiGen[nbin-ipad]->SetMarkerStyle(22);
                PsiGen[nbin-ipad]->SetMarkerColor(2);
                PsiGen[nbin-ipad]->SetLineColor(2);
                PsiGen[nbin-ipad]->Draw("same PE") ;
                //   PsiGen[nbin-ipad]->Print();
                
                bkgRadiusPsi[nbin-ipad] = new TGraphErrors(nrbin, rad, bkgpsi[nbin-ipad], errR, errbkgpsi[nbin-ipad]);
                bkgRadiusPsi[nbin-ipad]->SetMarkerStyle(24);
                bkgRadiusPsi[nbin-ipad]->SetMarkerColor(1);
                bkgRadiusPsi[nbin-ipad]->SetLineColor(1);
                bkgRadiusPsi[nbin-ipad]->Draw("same PE") ;
                
                bkgPsiGen[nbin-ipad] = new TGraphErrors(nrbin, rad, genbkgpsi[nbin-ipad], errR, errgenbkgpsi[nbin-ipad]);
                //RhoGen[nbin-ipad]->Print();
                bkgPsiGen[nbin-ipad]->SetMarkerStyle(26);
                bkgPsiGen[nbin-ipad]->SetMarkerColor(2);
                bkgPsiGen[nbin-ipad]->SetLineColor(2);
                bkgPsiGen[nbin-ipad]->Draw("same PE") ;
                
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
            
            t1->AddEntry(RadiusPsi[nbin-1],"RECO RAW", "P");
            t1->AddEntry(bkgRadiusPsi[nbin-1],"RECO Bkg", "P");  
            t1->AddEntry(PsiGen[nbin-1],"GEN RAW", "P");  
            t1->AddEntry(bkgPsiGen[nbin-1],"GEN Bkg", "P");
            t1->Draw("same");
            //            if(nbin==1){
            //                c1->cd(2);
            //                dummy->SetAxisRange(0.8, 1.2, "Y");
            //                dummy->SetYTitle("Gen_matched/Gen_only");
            //                dummy->DrawCopy();
            //                TGraphErrors * ratio = new TGraphErrors(nfile, rad, psiratio[0], errR, errpsiratio[0]);
            //                //                TH1F * ratio =(TH1F*)RhoGen[0]->Clone("ratio");
            //                //                ratio->Divide(RadiusRho[0]);
            //                ratio->Draw("same PE");
            //                
            //            }
            if(SavePlot)c1->Print(Form("%s/%sJetPt%.f_%.fTrkCut%.fIntJSwith%sBkg_%d.eps",plotsdir.Data(),type.Data(), pt[current],pt[current+1],trackcut, bkg.Data(), date));               
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
                
                GenRhoNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, genrhosub[nbin-ipad], errR, errgenrhosub[nbin-ipad]);
                //RhoGen[nbin-ipad]->Print();
                GenRhoNoBkg[nbin-ipad]->SetMarkerStyle(29);
                GenRhoNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                GenRhoNoBkg[nbin-ipad]->SetMarkerColor(1);
                GenRhoNoBkg[nbin-ipad]->SetLineColor(1);                   
                GenRhoNoBkg[nbin-ipad]->Draw("same PE") ;

                if(nbin>1){
                    if(ipad==1)drawCMSmc(0.5,0.85,type.Data());
                    if(ipad==2){
                        t1->AddEntry(RhoNoBkg[nbin-1],"Reco corrected", "PL");
                        t1->AddEntry(GenRhoNoBkg[nbin-1],"Gen level", "PL");  
                        t1->Draw("same");
                    }
                    //         if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.25,0.8,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.7,17);
                    if(ipad==3)drawText("AkPu3PF, R=0.3",0.45,0.8,17);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.7,17);
                }
                else {
                    drawCMSmc(0.8,0.9,type.Data());
                    t1->AddEntry(RhoNoBkg[nbin-1],"Corrected Reco", "PL");
                    t1->AddEntry(GenRhoNoBkg[nbin-1],"Gen level", "PL");  
                    t1->Draw("same");
                    drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.65,0.7,17);
                    drawText("Ak3PF, R=0.3",0.65,0.6,17);
                    drawText(Form("|#eta|_{jet} < %.f", etacut),0.65,0.5,17);
                    drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.65,0.4,17);
                }
                
            }
            
            for(int ibin = 0 ; ibin <nbin; ibin++){
                if(nbin>1)c1->cd(ibin+5);
                else c1->cd(2);
                dummy->SetAxisRange(0.8, 1.22, "Y") ;
                dummy->GetYaxis()->SetTitle("#rho(r)^{corrected}/#rho(r)^{gen}");
                dummy->DrawCopy();

                RatioRho[nbin-ibin-1] = new TGraphErrors(nrbin, rad, rhoratio[nbin-ibin-1], errR, errrhoratio[nbin-ibin-1]);
                RatioRho[nbin-ibin-1]->SetMarkerStyle(20);
                RatioRho[nbin-ibin-1]->SetMarkerColor(1);
                RatioRho[nbin-ibin-1]->SetLineColor(1);
                RatioRho[nbin-ibin-1]->SetMarkerSize(1.5);
                RatioRho[nbin-ibin-1]->Draw("same PE") ;
                regSun(0.,1.,0.3,1.,1, 1);
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.25,0.9, 17);
            }            
            if(SavePlot)c1->Print(Form("%s/%sJetPtThres%.f%s%.fGenRecDiffJSAfter%sBkgSub.pdf",plotsdir.Data(),type.Data(), pt[current],effTab.Data(), trackcut, bkg.Data(), date));               
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
                
                GenPsiNoBkg[nbin-ipad] = new TGraphErrors(nrbin, rad, genpsisub[nbin-ipad], errR, errgenpsisub[nbin-ipad]);
                GenPsiNoBkg[nbin-ipad]->SetMarkerStyle(25);
                GenPsiNoBkg[nbin-ipad]->SetMarkerColor(1);
                GenPsiNoBkg[nbin-ipad]->SetLineColor(1);
                GenPsiNoBkg[nbin-ipad]->SetMarkerSize(1.5);
                GenPsiNoBkg[nbin-ipad]->Draw("same PE") ;
                if(ipad==1)drawCMSmc(0.5,0.85,type.Data());
                if(nbin>1){
                    if(ipad==2){
                        t1->AddEntry(PsiNoBkg[nbin-1],"Reco corrcted", "PL");
                        t1->AddEntry(GenPsiNoBkg[nbin-1],"Gen level", "PL");  
                        t1->Draw("same");
                    }
                    //         if(ipad==4)drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.8,17);
                    if(ipad==4)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.25,0.85,17);
                    if(ipad==4)drawText(Form("|#eta|_{jet} < %.f", etacut),0.45,0.8,17);
                    if(ipad==3)drawText("Ak PF, R=0.3",0.55,0.85,17);
                    if(ipad==3)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.8,17);
                }
                else {
                    t1->AddEntry(PsiNoBkg[nbin-1],"Corrected Reco", "PL");
                    t1->AddEntry(GenPsiNoBkg[nbin-1],"Gen level", "PL");  
                    t1->Draw("same");
                    drawText(Form("%.f<p_{T}^{jet1}<%.f (GeV/c)", pt[current],pt[current+1]),0.6,0.75,17);
                }
                
            }
            
            for(int ibin = 0 ; ibin <nbin; ibin++){
                if(nbin>1)c1->cd(ibin+5);
                else c1->cd(2);
                dummy->SetAxisRange(0.7, 1.3, "Y") ;
                dummy->GetYaxis()->SetTitle("#psi(r)^{corrected}/#psi(r)^{gen}");
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
            if(SavePlot)c1->Print(Form("%s/%sJetPt%.fTrkCut%.fGenRecoIntJSAfter%sBkgSub.pdf",plotsdir.Data(),type.Data(), pt[current],trackcut, bkg.Data()));               
            c1->Update();               
            break ;
        case kBkgRatio:
            if(nbin>1)c1->Divide(nbin/2, 2,0, 0);
            //  else c1->Divide(1, 2,0, 0);           
            for(int ipad = 1 ; ipad <=nbin; ipad++){
                c1->cd(ipad);
             //   c1->cd(ipad)->SetLogy();
                dummy->GetYaxis()->SetTitle("ratio");
                //if(current==4)dummy->SetAxisRange(5.e-3, 20., "Y") ;
                // else 
                dummy->SetAxisRange(0., 2., "Y") ;
                dummy->DrawCopy();
                //   ppRadiusRho->Draw("same P") ;
                
                bkgRhoRatio[nbin-ipad] = new TGraph(nfile, rad, bkgrhoratio[nbin-ipad]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                bkgRhoRatio[nbin-ipad]->SetMarkerStyle(20);
                bkgRhoRatio[nbin-ipad]->SetMarkerColor(1);
                bkgRhoRatio[nbin-ipad]->SetLineColor(1);
                bkgRhoRatio[nbin-ipad]->Draw("same PE") ;

                RhoRatio[nbin-ipad] = new TGraph(nfile, rad, rhoratio[nbin-ipad]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                RhoRatio[nbin-ipad]->SetMarkerStyle(22);
                RhoRatio[nbin-ipad]->SetMarkerColor(2);
                RhoRatio[nbin-ipad]->SetLineColor(2);
                RhoRatio[nbin-ipad]->Draw("same PE") ;
                
                bkgPsiRatio[nbin-ipad] = new TGraph(nfile, rad, bkgpsiratio[nbin-ipad]);
                //               RadiusRho[nbin-ipad] = new TGraphErrors(nfile, ratio, rho[nbin-ipad], errR, errrho[nbin-ipad]);
                bkgPsiRatio[nbin-ipad]->SetMarkerStyle(24);
                bkgPsiRatio[nbin-ipad]->SetMarkerColor(4);
                bkgPsiRatio[nbin-ipad]->SetLineColor(4);
             //   bkgPsiRatio[nbin-ipad]->Draw("same PE") ;
                
                if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[ipad]),0.6,0.65, 15);
                if(ipad==1)drawText("Ak PF, R=0.3",0.4,0.68,15);
                if(nbin>1){
                    if(ipad==2)drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.45,0.75,15);
                    if(ipad==4)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.45,0.75,15);
                    
                }
                else {
                    drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.4,0.65,15);
                    drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.4,0.62,15);
                }
                
            }
            
            t1->AddEntry(RhoRatio[nbin-1],"#rho (r)", "P");
            t1->AddEntry(bkgRhoRatio[nbin-1]," bkg #rho (r)", "P");
        //    t1->AddEntry(bkgPsiRatio[nbin-1],"bkg #psi (r)", "P");
            t1->Draw("same");
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
            if(SavePlot)c1->Print(Form("%s/%sJetPt%.f_%.fTrkCut%.fDiffJSwith%sBkg_%d.pdf",plotsdir.Data(),type.Data(), pt[current],pt[current+1],trackcut, bkg.Data(), date));               
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
