//
//  DrawSmearingFactorPawan.C
//  
//
//  Created by Yaxian Mao on 12/3/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include "SmearingFactors.h"
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom.h>
#include <cstdlib>
#include <cmath>
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
#include "./rootlogon.h"    //// A good drawing style is defined here
#include<cstdlib>
#include<cstdio>
#include<cmath>
#include<iostream>
#include<fstream>
#include<string>
#include <iostream>
using namespace std;

void DrawSmearingFactorPawan()
{
    const int nbin = 6 ;
    const int centr[] ={0, 5, 10,30,50,70, 90}; 
    const Double_t pt[]={30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,400};
    const int nptbin = 33 ;
    
    //! Load Smearing factors
    LoadParameters();
    
    TH1F * smf[nbin];
    TH1F * ms[nbin] ;
    for(int ibin = 0 ; ibin <nbin; ibin++){
        smf[ibin] = new TH1F(Form("SmearingFactor_%d-%d%%",centr[ibin],centr[ibin+1]), Form("SmearingFactor_%d-%d%%",centr[ibin],centr[ibin+1]), nptbin, pt);
        smf[ibin]->Sumw2();
        ms[ibin] = new TH1F(Form("jetMeanShift_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetMeanShift_%d-%d%%",centr[ibin],centr[ibin+1]), nptbin, pt);
        ms[ibin]->Sumw2();
    }
    
    TString canv_name = "c1";
    const Double_t kw = 1000;
    const Double_t kh = 780;
    c1 = new TCanvas(canv_name," ",10,10,kw,kh);
    makeMultiPanelCanvas(c1,3,2,0.0,0.0,0.2,0.2,0.02); 
    //gStyle->SetOptStat(0);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadTopMargin   (0.025);
    gStyle->SetPadLeftMargin  (0.15);
    gStyle->SetPadRightMargin (0.025);
    gStyle->SetPadTickX       (1);
    gStyle->SetPadTickY       (1);
    
    TRandom ran ; 
    Double_t mpt ;
    Double_t jersmf ;
    Double_t jesmsf ;
    for(int ibin = 0 ; ibin <nbin; ibin++){ 
        for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
            jersmf = 0. ;
            jesmsf = 0. ;
            Int_t kCounts = 500; 
            for (Int_t index = 1; index < =kCounts; index++) {
            mpt = ran.Uniform(pt[ipt],pt[ipt+1]); 
            jersmf = GetSmFactor(2, ibin, mpt);
            jesmsf = GetMeanShift(2, ibin, mpt);
                cout<<"jer =" <<jesmsf*mpt<<endl;
            smf[ibin]->Fill(mpt, jersmf);
            ms[ibin]->Fill(mpt, jesmsf*mpt);
            }
        }
//        smf[ibin]->Scale(1/kCounts);
//        ms[ibin]->Scale(1/kCounts);
//        for(Int_t i = 1 ; i < =nptbin ; i++){
//            smf[ibin]->SetBinContent(i,smf[ibin]->GetBinContent(i)/smf[ibin]->GetBinWidth(i)); 
//            smf[ibin]->SetBinError(i,smf[ibin]->GetBinError(i)/smf[ibin]->GetBinWidth(i)); 
//            ms[ibin]->SetBinContent(i,ms[ibin]->GetBinContent(i)/ms[ibin]->GetBinWidth(i)); 
//            ms[ibin]->SetBinError(i,ms[ibin]->GetBinError(i)/ms[ibin]->GetBinWidth(i)); 
//        }
    }
    
    TH1F * dummy = new TH1F("dummy", "dummy", nptbin, pt);
    dummy->SetTitle("") ;
    dummy->SetStats(kFALSE) ;
    // dummy->GetYaxis()->SetTitle("counts");
    dummy->GetXaxis()->SetTitle("p_{T} GeV/c");
    dummy->SetAxisRange(30., 400., "X") ;
    //   dummy->SetAxisRange(1.e-5, 1., "Y") ;
    fixedFontHist(dummy,1.5, 2.0);
    for(int ipad = 1 ; ipad <=nbin; ipad++){        
        c1->cd(ipad)->SetLogy();
        dummy->SetAxisRange(1.e-6., 1.e2, "Y") ;
        dummy->DrawCopy();
//        smf[nbin-ipad]->SetMarkerStyle(20);
//        smf[nbin-ipad]->SetMarkerSize(1.5);
//        smf[nbin-ipad]->SetMarkerColor(1);
//        smf[nbin-ipad]->SetLineColor(1);
//        smf[nbin-ipad]->DrawCopy("same");
        ms[nbin-ipad]->SetMarkerStyle(20);
        ms[nbin-ipad]->SetMarkerSize(1.5);
        ms[nbin-ipad]->SetMarkerColor(1);
        ms[nbin-ipad]->SetLineColor(1);
        ms[nbin-ipad]->Draw("same");
        drawText(Form("%d-%d%%",centr[nbin-ipad],centr[nbin-ipad+1]),0.2,0.92, 15);

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
