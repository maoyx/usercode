/*
 *  Draw the background subtracted JS ratio from event mixing and eta reflection
 *  outputs from Pelin
 *
 *  
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
const double dphicut = 7*TMath::TwoPi()/8. ;
const double leadingjetcut = 100. ;
const double subleadjetcut = 40. ;
const double trackcut = 1.0;
bool SavePlot = kTRUE ;

const bool IsMC=kTRUE ;
const bool SaveFile=kFALSE ;

const int nfile = 6 ; //8 until 0.3, 13 until 0.5
//const int nfile = (conesize-deltacone/2.)/deltacone +1 ;
//const double rad[] = {0.05, 0.15, 0.25};
//const double ratio[] = {0.05/conesize, 0.15/conesize, 0.25/conesize};
//const double errR[] = {0.05, 0.05, 0.05};
double rad[nfile];
double ratio[nfile];
double errR[nfile]; 
const int nptbin = 6 ;
//const double pt[]={100., 110., 120.,130., 140., 150.,160.,200., 500.};
//const double pt[]={100., 120., 140., 160., 200., 300., 500.};
const double pt[]={100., 300.};

//const double pt[]={100., 120., 140.,500.};
const double subpt[]={40., 50., 60., 70., 80., 100., 120};
const double radRebin[] = {0.05, 0.15, 0.25};

const int nrbin = 6 ;
double rbin[nrbin+1]={0.0,0.05,0.10, 0.15, 0.20,0.25,0.30};

const double kAj[]={0.0,0.1,0.4,1.0};

const Double_t pbpbLumi = 150 ;
const Double_t ppLumi = 231 ;

//---------------------------------------------------------------------

void DrawTrackClosureDoubleRatio(){
    for(int i = 0 ; i <nrbin; i++){
        //        if(opt==kInt) rad[i]=rbin[i+1] ;
        //        else         
        rad[i]=deltacone/2.+i*deltacone ;
        ratio[i]=rad[i]/conesize;
        errR[i]=deltacone/2. ;
    }
    TString plotsdir = "/Users/ymao/group/CMS/plots/finalJS/CorrectDR";
//    TString coll ;
//    if(data==kPP) coll="PP";
//    else coll = "HI" ;
//    if(data==kPP) {
//        ////if it is pp, no centrality bins, only one
//        const int nbin = 1 ;
//        const int centr[] ={0,100};  
//    }
//    else {
        // for HI centrality bin
        const int nbin = 4 ;
        const int centr[] ={0,10,30,50,100}; 
//    }
    const char * kHomeDir = "/Users/ymao/group/CMS/macros" ;      
 //   const char * kHomeDir2 = "/Users/ymao/group/CMS/macros" ;      
    TString fileName ;
    TFile * f ;

    TString JSfileName ;
    TFile * JSf ;

    TString ppfileName ;
    TFile * ppf ;
    
    TString ppJSfileName ;
    TFile * ppJSf ;

    TH1F * pptrkjetcone ;
    TH1F * pptrkbkgcone ;
    
    TGraphErrors * pprhojetcone ;
    TGraphErrors * pprhobkgcone ;

    TH1F * trkjetcone[nbin] ;
    TH1F * trkbkgcone[nbin] ;
    
    TGraphErrors * rhojetcone[nbin] ;
    TGraphErrors * rhobkgcone[nbin] ;
    
    TH1F * trkjetconeratio[nbin] ;
    TH1F * trkbkgconeratio[nbin] ;
    
    TGraph * rhojetconeratio[nbin] ;
    TGraph * rhobkgconeratio[nbin] ;
    

    double ppjetrho[nrbin] ;
    double ppbkgrho[nrbin] ;
    double jetrho[nbin][nrbin];
    double bkgrho[nbin][nrbin];
    double tmp[nbin][nrbin];
    
    double jetrhoratio[nbin][nrbin];
    double bkgrhoratio[nbin][nrbin];
    for(int ir =0 ; ir <6; ir++){ 
        ppjetrho[ir]=0.;
        ppbkgrho[ir]=0.;
    }
    for(int ibin = 0 ; ibin <nbin; ibin++){  
        for(int ir =0 ; ir <6; ir++){ 
            tmp[ibin][ir]=0.;
            jetrho[ibin][ir]= 0.;
            bkgrho[ibin][ir] =0. ;
            jetrhoratio[ibin][ir]= 0.;
            bkgrhoratio[ibin][ir] =0. ;            
        }
    }
    
    TH1::SetDefaultSumw2();
    
    TString canv_name = "c1";
//    if(data==kHI){
        const Double_t kw = 1000;
        const Double_t kh = 560;
        c1 = new TCanvas(canv_name," ",10,10,kw,kh);
        //    makeMultiPanelCanvas(c1,4,1,0.0,0.0,0.2,0.2,0.02);
        makeMultiPanelCanvas(c1,nbin,1,0.0,0.0,0.2,0.2,0.02);
//    }
//    else {
//        const Double_t kw = 800;
//        const Double_t kh = 800;
//        c1 = new TCanvas(canv_name," ",10,10,kw,kh);
//        //    makeMultiPanelCanvas(c1,4,1,0.0,0.0,0.2,0.2,0.02);
//        makeMultiPanelCanvas(c1,nbin,1,0.0,0.0,0.1,0.2,0.05);  
//    }
    
    gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);
    
    
    if(nbin>1)TLegend *t1=new TLegend(0.25,0.7,0.45,0.93);
    else TLegend *t1=new TLegend(0.35,0.75,0.55,0.93);

    //   TLegend *t1=new TLegend(0.25,0.6,0.8,0.92);
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
    
    TString effTab = "HistIterTrkCorrtest"; // "HistIterTrkCorrtest" ; // "HistIterTrkCorrv14fXSec";

    ppfileName =  Form("MCPP%swtTrackRadius.root", effTab.Data()); 
    ppJSfileName =  Form("PPIncJet%sDiffJSRatiowithERbkg.root", effTab.Data()); 

    fileName =  Form("MCHI%swtTrackRadius.root", effTab.Data()); 
    JSfileName =  Form("HIIncJet%sDiffJSRatiowithERbkg.root",  effTab.Data()); 
    
    ppf = TFile::Open(Form("%s/%s", kHomeDir, ppfileName.Data()), "readonly");
    ppJSf = TFile::Open(Form("%s/%s", kHomeDir, ppJSfileName.Data()), "readonly");

    f = TFile::Open(Form("%s/%s", kHomeDir, fileName.Data()), "readonly");
    JSf = TFile::Open(Form("%s/%s", kHomeDir, JSfileName.Data()), "readonly");

    pptrkjetcone=(TH1F*)ppf->Get(Form("wtTrkCut1_300ClosureRadiusJet_Cen%d-%d%%", centr[0], centr[nbin]));  
    pptrkjetcone->Sumw2();
    
    pptrkbkgcone=(TH1F*)ppf->Get(Form("wtTrkCut1_300ClosureRadiusBkg_Cen%d-%d%%", centr[0], centr[nbin]));  
    pptrkbkgcone->Sumw2();
    
    pprhojetcone= (TGraphErrors*)ppJSf->Get(Form("JetConeRhoERbkg_Cen%d-%d", centr[0], centr[nbin])); 
    pprhobkgcone= (TGraphErrors*)ppJSf->Get(Form("BkgConeRhoERbkg_Cen%d-%d", centr[0], centr[nbin]));    
    for(int ir =0 ; ir <6; ir++){ 
        pprhojetcone->GetPoint(ir, tmp[0][ir], ppjetrho[ir]); 
        pprhobkgcone->GetPoint(ir, tmp[0][ir], ppbkgrho[ir]); 
    }
    
    for(int ibin = 0 ; ibin <nbin; ibin++){  
        trkjetcone[ibin]=(TH1F*)f->Get(Form("wtTrkCut1_300ClosureRadiusJet_Cen%d-%d%%", centr[ibin], centr[ibin+1]));  
        trkjetcone[ibin]->Sumw2();

        trkbkgcone[ibin]=(TH1F*)f->Get(Form("wtTrkCut1_300ClosureRadiusBkg_Cen%d-%d%%", centr[ibin], centr[ibin+1]));  
        trkbkgcone[ibin]->Sumw2();

        rhojetcone[ibin]= (TGraphErrors*)JSf->Get(Form("JetConeRhoERbkg_Cen%d-%d", centr[ibin], centr[ibin+1])); 
        rhobkgcone[ibin]= (TGraphErrors*)JSf->Get(Form("BkgConeRhoERbkg_Cen%d-%d", centr[ibin], centr[ibin+1]));    
        
        trkjetconeratio[ibin] =(TH1F*)trkjetcone[ibin]->Clone(Form("wtTrkCut1_300ClosureJetRadiusRatio_Cen%d-%d%%", centr[ibin], centr[ibin+1]));
        trkjetconeratio[ibin]->Divide(pptrkjetcone);

        trkbkgconeratio[ibin] =(TH1F*)trkbkgcone[ibin]->Clone(Form("wtTrkCut1_300ClosureBkgRadiusRatio_Cen%d-%d%%", centr[ibin], centr[ibin+1]));
        trkbkgconeratio[ibin]->Divide(pptrkbkgcone);
        
        for(int ir =0 ; ir <6; ir++){ 
            rhojetcone[ibin]->GetPoint(ir, tmp[ibin][ir], jetrho[ibin][ir]); 
            rhobkgcone[ibin]->GetPoint(ir, tmp[ibin][ir], bkgrho[ibin][ir]); 
            
        }
        for(int ir =0 ; ir <6; ir++){ 
            jetrhoratio[ibin][ir]=jetrho[ibin][ir]/ppjetrho[ir];
            bkgrhoratio[ibin][ir]=bkgrho[ibin][ir]/ppbkgrho[ir];
        }
        
        rhojetconeratio[ibin] =new TGraph(nrbin, rad, jetrhoratio[ibin]);
                
        rhobkgconeratio[ibin] =new TGraph(nrbin, rad, bkgrhoratio[ibin]);     
    }
    
    TH1F * dummy = new TH1F("dummy", "dummy", 100, 0., 1.);
    dummy->SetTitle("") ;
    dummy->SetStats(kFALSE) ;
    dummy->SetAxisRange(0., conesize, "X") ;
    dummy->GetXaxis()->SetTitle("radius (r)");
    fixedFontHist(dummy,1.1, 1.3);
    
    for(int ibin = 0 ; ibin <nbin; ibin++){
        c1->cd(ibin+1);
        dummy->SetAxisRange(0.9, 1.12, "Y") ;
//        dummy->SetAxisRange(0.8, 1.22, "Y") ;
        dummy->GetYaxis()->SetTitle("PbPb/pp");
        dummy->DrawCopy();
        trkjetconeratio[nbin-ibin-1]->SetMarkerStyle(24);
        trkjetconeratio[nbin-ibin-1]->SetMarkerColor(1);
        trkjetconeratio[nbin-ibin-1]->SetLineColor(1);
        trkjetconeratio[nbin-ibin-1]->SetMarkerSize(1.5);
        
        trkjetconeratio[nbin-ibin-1]->Draw("same PE") ;

        rhojetconeratio[nbin-ibin-1]->SetMarkerStyle(29);
        rhojetconeratio[nbin-ibin-1]->SetMarkerColor(1);
        rhojetconeratio[nbin-ibin-1]->SetLineColor(1);
        rhojetconeratio[nbin-ibin-1]->SetMarkerSize(1.5);
        rhojetconeratio[nbin-ibin-1]->Draw("same PE") ;
        
        rhojetconeratio[nbin-ibin-1]->Print();

        trkbkgconeratio[nbin-ibin-1]->SetMarkerStyle(24);
        trkbkgconeratio[nbin-ibin-1]->SetMarkerColor(2);
        trkbkgconeratio[nbin-ibin-1]->SetLineColor(2);
        trkbkgconeratio[nbin-ibin-1]->SetMarkerSize(1.5);
        
    //    trkbkgconeratio[nbin-ibin-1]->Draw("same PE") ;
        
        rhobkgconeratio[nbin-ibin-1]->SetMarkerStyle(29);
        rhobkgconeratio[nbin-ibin-1]->SetMarkerColor(2);
        rhobkgconeratio[nbin-ibin-1]->SetLineColor(2);
        rhobkgconeratio[nbin-ibin-1]->SetMarkerSize(1.5);
        
     //   rhobkgconeratio[nbin-ibin-1]->Draw("same PE") ;
        rhobkgconeratio[nbin-ibin-1]->Print();
        regSun(0.,1.,0.3,1.,1, 1);
        if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.25,0.25, 17);
    //    if(ibin==0)drawCMS(0.5,0.9,pbpbLumi);
        //    if(ipad==1)drawCMS(0.5,0.85,ppLumi);
        //     if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[nbin-ipad+1]),0.65,0.8, 17);
      //  if(ibin==1){

 //       }
        if(nbin>1){
            if(ibin==0)drawCMSmc(0.5,0.9,"HI");
            if(ibin==1)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.25,0.9,17);
            if(ibin==1)drawText(Form("|#eta|_{jet} < %.f", etacut),0.25,0.83,17);
            if(ibin==1)drawText("AkPu3PF, R=0.3",0.25,0.75,17);
            if(ibin==2)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.25,0.9,17);
        }
        else {
            drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.55,0.28,17);
            drawText(Form("|#eta|_{jet} < %.f", etacut),0.77,0.28,17);
             drawText("Ak3PF, R=0.3",0.55,0.35,17);
             drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.22,17);
        }
        
    }            
    t1->AddEntry(trkjetconeratio[nbin-1],"dN/dr: jet cone", "PL");
    t1->AddEntry(rhojetconeratio[nbin-1],"#rho(r): jet cone", "PL");
//    t1->AddEntry(trkbkgconeratio[nbin-1],"dN/dr: bkg cone", "PL");
//    t1->AddEntry(rhobkgconeratio[nbin-1],"#rho(r): bkg cone", "PL");
    t1->Draw("same");   
    if(SavePlot)
        c1->Print(Form("%s/JetTrkCut%.f%sEffUncertaintyDoubleRatio.pdf",plotsdir.Data(),trackcut, effTab.Data()));  
    
    if(SaveFile){
        TFile * outf = new TFile(Form("%sIncJet%sTrkDoubleRatio.root", coll.Data(), effTab.Data()), "CREATE");
        for(int ibin = 0 ; ibin <nbin; ibin++){                    
            trkjetconeratio[ibin]->Write();
            rhojetconeratio[ibin]->Write();
            trkbkgconeratio[ibin]->Write();
            rhobkgconeratio[ibin]->Write();
        }
        outf->Close();  
    }
    c1->Update();   
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
