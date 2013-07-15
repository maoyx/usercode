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

const bool IsMC=kFALSE ;
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
//const double pt[]={100., 120., 150., 200., 300.};

//const double pt[]={100., 120., 140.,500.};
const double subpt[]={40., 50., 60., 70., 80., 100., 120};
const double radRebin[] = {0.05, 0.15, 0.25};

const int nrbin = 6 ;
double rbin[nrbin+1]={0.0,0.05,0.10, 0.15, 0.20,0.25,0.30};

const double kAj[]={0.0,0.1,0.4,1.0};

const Double_t pbpbLumi = 150 ;
const Double_t ppLumi = 231 ;

//---------------------------------------------------------------------

void DrawFinalJSDataMC(int current = 0){
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

    TString mcfileName ;
    TFile * mcf ;

    TGraphErrors * rhodata[nbin] ;
    TGraphErrors * rhomc[nbin] ;
        

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
        const Double_t kh = 400;
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
    
    
    if(nbin>1)TLegend *t1=new TLegend(0.15,0.85,0.35,0.95);
    else TLegend *t1=new TLegend(0.15,0.85,0.35,0.93);

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


    fileName =  Form("DATARECOIncJetTrkEffCorrectedTrkHistIterTrkCorrtestDiffJSbkgSub.root"); 
 //   mcfileName =  Form("DATARECOTrig80IncJetTrkEffCorrectedTrkHistIterTrkCorrtestDiffJSbkgSub.root");
 //   mcfileName =  Form("MCRECOIncJetTrkEffCorrectedTrkHistIterTrkCorrtestDiffJSbkgSub.root");
//    mcfileName =  Form("DATARECOIncJetTrkEffCorrectedTrkHistIterTrkCorrtestDiffJSbkgSub2011Ref.root");
    mcfileName =  Form("DATARECOIncJetNoTrkEffCorrTrkTrkDiffJSbkgSub.root"); 

    f = TFile::Open(Form("%s/%s", kHomeDir, fileName.Data()), "readonly");
    mcf = TFile::Open(Form("%s/%s", kHomeDir, mcfileName.Data()), "readonly");
    
    for(int ibin = 0 ; ibin <nbin; ibin++){  

        rhodata[ibin]= (TGraphErrors*)f->Get(Form("RatioRhoERbkgsubJetPt%.f_%.f_Cen%d-%d", pt[current],pt[current+1], centr[ibin],centr[ibin+1])); 
        rhomc[ibin]= (TGraphErrors*)mcf->Get(Form("RatioRhoERbkgsubJetPt%.f_%.f_Cen%d-%d", pt[current],pt[current+1], centr[ibin],centr[ibin+1]));    
        
    }
    
//    double syspercent[nbin][nrbin]={{1.26314, 2.7714, 5.15745, 7.18162, 8.69646, 11.1731},
//        {1.24532, 2.90773, 5.23464, 7.38403, 8.96785, 10.9962},
//        {1.22754, 1.58031, 3.86354, 3.26574, 5.39212, 4.74478},
//        {1.18007, 1.40265, 2.58602, 2.74553, 2.79803, 3.08416}};   
    //! when using 2013 pp reference
    double syspercent[nbin][nrbin]={{1.72241, 3.00862, 5.28871, 7.27646, 8.77494, 11.2343},
        {1.70939, 3.13465, 5.36401, 7.4763, 9.04397, 11.0584},
        {1.69647, 1.96686, 4.0371, 3.46932, 5.5178, 4.88714},
        {1.66245, 1.82719, 2.83878, 2.98481, 3.03317, 3.29897}};

    
    TH1F * dummy = new TH1F("dummy", "dummy", 100, 0., 1.);
    dummy->SetTitle("") ;
    dummy->SetStats(kFALSE) ;
    dummy->SetAxisRange(0., conesize, "X") ;
    dummy->GetXaxis()->SetTitle("radius (r)");
    fixedFontHist(dummy,1.1, 1.3);
    
    for(int ibin = 0 ; ibin <nbin; ibin++){
        c1->cd(ibin+1);
        dummy->SetAxisRange(0.3, 1.7, "Y") ;
        dummy->GetYaxis()->SetTitle("PbPb/pp");
        dummy->DrawCopy();
        rhodata[nbin-ibin-1]->SetMarkerStyle(20);
        rhodata[nbin-ibin-1]->SetMarkerColor(1);
        rhodata[nbin-ibin-1]->SetLineColor(1);
        rhodata[nbin-ibin-1]->SetMarkerSize(1.5);
        drawSys(rhodata[nbin-ibin-1], syspercent[nbin-ibin-1],deltacone/2.,7, 1001, 1);  
        
        rhodata[nbin-ibin-1]->Draw("same PE") ;

       rhomc[nbin-ibin-1]->SetMarkerStyle(24);
        rhomc[nbin-ibin-1]->SetMarkerColor(2);
        rhomc[nbin-ibin-1]->SetLineColor(2);
        rhomc[nbin-ibin-1]->SetMarkerSize(1.5);
 //       drawSys(rhomc[nbin-ibin-1], syspercent[nbin-ibin-1],deltacone/2.,7, 1001, 1);  
        rhomc[nbin-ibin-1]->Draw("same PE") ;
        
        rhomc[nbin-ibin-1]->Print();

        regSun(0.,1.,0.3,1.,1, 1);
        if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ibin-1],centr[nbin-ibin]),0.25,0.25, 17);
    //    if(ibin==0)drawCMS(0.5,0.9,pbpbLumi);
        //    if(ipad==1)drawCMS(0.5,0.85,ppLumi);
        //     if(nbin>1)drawText(Form("%d-%d%%",centr[nbin-ipad],centr[nbin-ipad+1]),0.65,0.8, 17);
      //  if(ibin==1){

 //       }
        if(nbin>1){
            if(ibin==0)drawCMS(0.5,0.9,pbpbLumi);
            if(ibin==1)drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.25,0.9,17);
//            if(ibin==1)drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.25,0.9,17);
            if(ibin==1)drawText(Form("|#eta|_{jet} < %.f", etacut),0.25,0.83,17);
            if(ibin==1)drawText("Anti-k_{T} PF, R=0.3",0.25,0.75,17);
            if(ibin==2)drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.25,0.9,17);
        }
        else {
            drawText(Form("%.f<p_{T}^{jet}<%.f (GeV/c)", pt[current],pt[current+1]),0.55,0.28,17);
//            drawText(Form("p_{T}^{jet}>%.f GeV/c", leadingjetcut),0.55,0.28,17);
            drawText(Form("|#eta|_{jet} < %.f", etacut),0.77,0.28,17);
             drawText("Ak3PF, R=0.3",0.55,0.35,17);
             drawText(Form("p_{T}^{trk} >%.f GeV/c",trackcut),0.55,0.22,17);
        }
        
    }            
//    t1->AddEntry(rhodata[nbin-1],"Data", "PL");
// //   t1->AddEntry(rhomc[nbin-1],"MC: Gen", "PL");
//    t1->AddEntry(rhomc[nbin-1],"MC: Reco", "PL");
    t1->AddEntry(rhodata[nbin-1],"W/ Trk eff. Corr", "PL");
    t1->AddEntry(rhomc[nbin-1],"w/o Trk eff. Corr", "PL");
//    t1->AddEntry(rhodata[nbin-1],"HLT jet >90 GeV/c", "PL");
//    t1->AddEntry(rhomc[nbin-1],"HLT jet > 80 GeV/c", "PL");
    t1->Draw("same");
    if(SavePlot)
        c1->Print(Form("%s/JetTrkCut%.fJSRatioJetPt%.f_%.fWithWithoutEff.gif",plotsdir.Data(),trackcut, pt[current],pt[current+1]));
    
    if(SaveFile){
        TFile * outf = new TFile(Form("%sIncJet%sTrkDoubleRatio.root", coll.Data(), effTab.Data()), "CREATE");
        for(int ibin = 0 ; ibin <nbin; ibin++){                    
            rhodata[ibin]->Write();
            rhodataratio[ibin]->Write();
            trkmcratio[ibin]->Write();
            rhomcratio[ibin]->Write();
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
