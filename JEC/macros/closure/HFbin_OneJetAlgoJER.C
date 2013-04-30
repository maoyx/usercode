#include <cstdlib>
#include <cmath>
#include <iostream>
#include <TROOT.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TLatex.h>
#include "TCanvas.h"
#include <TMarker.h>
#include <TString.h>
#include "MultiCanvas.h"
using namespace std;

const double pi=acos(-1.);
const double pi2=2*pi -1;

int bins =500;
float minval = 0;
float maxval = 1000;
float binw   = 2;

const int knj = 16;
enum Algo_t {kPF, kCalo};
enum Display_t {kMultJER, kMultJES, kHFJER, kHFJES} ;
enum Generator_t {kHIJING, kPYTHIA};

const int nmult=6;
const char *cmult[nmult] = {"N_{trk}^{offline}<60",
			    "60<N_{trk}^{offline}<90",
			    "90<N_{trk}^{offline}<110",
			    "110<N_{trk}^{offline}<150",
			    "150<N_{trk}^{offline}<180",
			    "180<N_{trk}^{offline}"
};
const Int_t multBin[nmult+1] = {0, 60, 90, 110, 150, 180, 500};
const int nhf=6;
const Double_t hfBin[nhf+1] = {0., 5., 10., 15., 20., 30., 100};
//const char *chf[nhf] = {"E_{T}^{HF[#eta>4]}<5",
//    "5<E_{T}^{HF[#eta>4]}<10",
//    "10<E_{T}^{HF[#eta>4]}<15",
//    "15<E_{T}^{HF[#eta>4]}<20",
//    "20<E_{T}^{HF[#eta>4]}<30",
//    "30<E_{T}^{HF[#eta>4]}"
//};
const char *chf[nhf] = {"E_{T}^{HF[|#eta|>4]}<20",
    "20<E_{T}^{HF[|#eta|>4]}<25",
    "25<E_{T}^{HF[|#eta|>4]}<30",
    "30<E_{T}^{HF[|#eta|>4]}<40",
    "E_{T}^{HF[|#eta|>4]}>40",
    "inclusive"
};
const char *ksp  [nmult] = {"pPb"    ,"pPb"       ,"pPb"       ,"pPb"        ,"pPb"         ,"pPb"};

//const char *fopt="MLRQ+"; 
int iFit=0; 
const char *fopt="RQ+";
const int knpx=2000;
float fitmin=0.00;
float fitmax=5.00;

double xmin=30.;
double xmax=340.;
double xfitmin=0.3;
double xfitmax=2.0;
int maxEntry=5;

int GetPtBin(float /*pt*/);
void MakeHist(TH1 */*hist*/,int /*istat*/);
void MakeHistRMS(TH1 */*hRMS*/,float /*max*/,float /*min*/);
void MakeHistMean(TH1 */*Mean*/,float /*max*/,float /*min*/);
void CleanHist(TH1 */*h1D*/,float /*lxrange*/,float /*hxrange*/);
void FillMeanSigma(int /*ip*/,TH1 */*h1*/,TH1 */*ArM*/,TH1 */*RMS*/,TH1 */*Mean*/,TH1 */*Sigma*/);

int HFbin_OneJetAlgoJER(Algo_t algo = kPF, Display_t dis=kHFJES, Generator_t cmode=kHIJING, int rfit=0)
{
  int statop=0;
  int kRebin=1;
  const char *reta = "eta2.0";
  //bool iSave=false;
  bool iMult=false;

    TString plotsdir = "/Users/ymao/group/CMS/plots/JEC";
    bool SavePlot = kTRUE ;
    const bool SaveFile = kFALSE ;
    
    TString gen ;
    if(cmode==kPYTHIA) gen="PYTHIA";
    else gen="PYTHIAHIJING";    

    TString algn ;
    Int_t nj ;
    switch(algo){
        case kPF:
            algn="PF";
            nj = 1 ;
            break ;
        case kCalo:
            algn="Calo";
            nj = 9 ;
            break ;
    }
    
  float ketacut=2.0;
  if(strcmp(reta,"eta2.0")==0)ketacut=2.5;
  bool iSigma=false;

  if(rfit==0){fitmin=0.01;fitmax=2.00;}   
  else if(rfit==1){fitmin=0.00;fitmax=1.25;}
  else if(rfit==2){fitmin=0.00;fitmax=1.50;}
  else if(rfit==3){fitmin=0.00;fitmax=1.75;}
  else if(rfit==4){fitmin=0.00;fitmax=2.50;}
  else if(rfit==5){fitmin=0.00;fitmax=3.50;} 
  else if(rfit==6){fitmin=0.50;fitmax=1.25;}
  else if(rfit==7){fitmin=0.50;fitmax=1.50;}
  else if(rfit==8){fitmin=0.50;fitmax=1.75;}
  else if(rfit==9){fitmin=0.50;fitmax=2.50;}
  else if(rfit==10){fitmin=0.50;fitmax=3.50;}
  else if(rfit==11){fitmin=0.75;fitmax=1.25;}
  else if(rfit==12){fitmin=0.75;fitmax=1.50;}
  else if(rfit==13){fitmin=0.75;fitmax=1.75;}
  else if(rfit==14){fitmin=0.75;fitmax=2.50;}
  else if(rfit==15){fitmin=0.75;fitmax=3.50;}
  else if(rfit==16){fitmin=1.00;fitmax=1.25;}
  else if(rfit==17){fitmin=1.00;fitmax=1.50;}
  else if(rfit==18){fitmin=1.00;fitmax=1.75;}
  else if(rfit==19){fitmin=1.00;fitmax=2.50;}
  else if(rfit==20){fitmin=1.00;fitmax=3.50;}

  if(kRebin){
    if(bins%kRebin!=0){
      cout<<"Cannot be divided in these bins chose another combination : "<<endl;
      return 0;
    }
    bins /= kRebin;
    binw = (maxval  - minval)/(1.0*bins);
    cout<<"kRebin : "<<kRebin<<"\t bins : "<<bins<<"\t binw  : "<<binw<<endl;
  }
 
    TFile *fin_pPb ;
    //! Input files 
    if(cmode==kPYTHIA)
        fin_pPb     = TFile::Open("../Output/pPbResponse_PYTHIA_merged.root","r");
    else 
        fin_pPb     = TFile::Open("../Output/pPbResponse_PYTHIAHIJING_merged.root","r");
    
  cout<<"\t"<<endl
      <<"rfit : "<<rfit<<"\t fitmin : "<<fitmin<<"\t fitmax : "<<fitmax<<endl
      <<"Input file name  pPb   : "<<fin_pPb->GetName()<<endl
    <<"# of pt bins : "<<bins<<endl
      <<"\t"<<endl;
  //return 0;

  double ptbins[] = {30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,340};
  const int histb  = sizeof(ptbins)/sizeof(Double_t) - 1;
  const int nbins = histb;
  cout<<"nbins = "<<nbins<<endl;
    cout <<"nj ="<<nj<<endl ;
  //return 0;
  int maxr=2;

  //xfitmin = (ptbins[1]+ptbins[2])/2.;
  //xfitmax = (ptbins[nbins] + ptbins[nbins-1])/2.;


  int ipad=0;
  gStyle->SetPadTickY(1);
  gStyle->SetPadTickX(1);
  gStyle->SetErrorX(0);
  gStyle->SetCanvasColor(10);
  gStyle->SetFillColor(10);
  gStyle->SetFrameFillColor(10);
  gStyle->SetPadColor(10);
  gStyle->SetStatColor(10);
  gStyle->SetTitleFillColor(10);
  gStyle->SetPadBorderSize(0);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetPalette(1);
  gStyle->SetPadBottomMargin(0.1);
  gStyle->SetPadLeftMargin(0.1);
  gStyle->SetOptStat(1);
  gStyle->SetOptFit(1);
  gStyle->SetPadBorderMode(0);
  
  //! pp Response
  TH2F *hratiocorrrefpt_pp[nmult];
  TH1F *hratiocorrrefpt1D_pp[nmult][nbins];  
    TH2F *hPuratiocorrrefpt_pp[nmult];
    TH1F *hPuratiocorrrefpt1D_pp[nmult][nbins];  
    TH1F *hMean_pp [nmult], *hArM_pp[nmult], *hSigma_pp [nmult], *hRMS_pp[nmult];
  TH1F *hPuMean_pp [nmult], *hPuArM_pp[nmult], *hPuSigma_pp [nmult], *hPuRMS_pp[nmult];

    TH2F *hratiocorrrefpt_pphf[nhf];
    TH1F *hratiocorrrefpt1D_pphf[nhf][nbins];  
    TH1F *hMean_pphf[nhf], *hArM_pphf[nhf], *hSigma_pphf[nhf], *hRMS_pphf[nhf];
    
    TH2F *hPuratiocorrrefpt_pphf[nhf];
    TH1F *hPuratiocorrrefpt1D_pphf[nhf][nbins];  
    TH1F *hPuMean_pphf[nhf], *hPuArM_pphf[nhf], *hPuSigma_pphf[nhf], *hPuRMS_pphf[nhf];
 

//  for(int nj=0;nj<knj;nj++){
    //cout<<"nj : "<<nj<<Form("\t %s",algn.Data())<<endl;
    
    
    for(int imult=0;imult<nmult;imult++){
      //! pp /////////////////////////////
      hMean_pp [imult] = new TH1F(Form("hMean_pp_%d",imult),Form("Mean <reco p_{T}/gen p_{T}> %s pp %d",algn.Data(),imult),nbins,ptbins);
      MakeHist(hMean_pp[imult],statop);

      hArM_pp [imult] = new TH1F(Form("hArM_pp_%d",imult),Form("<reco p_{T}/gen p_{T}> %s pp %d",algn.Data(),imult),nbins,ptbins);
      MakeHist(hArM_pp[imult],statop);

      hSigma_pp [imult] = new TH1F(Form("hSigma_pp_%d",imult),Form("#sigma(reco p_{T}/gen p_{T}) %s pp %d",algn.Data(),imult),nbins,ptbins);
      MakeHist(hSigma_pp[imult],statop);    

      hRMS_pp [imult] = new TH1F(Form("hRMS_pp_%d",imult),Form("RMS(reco p_{T}/gen p_{T}) %s pp %d",algn.Data(),imult),nbins,ptbins);
      MakeHist(hRMS_pp[imult],statop);        

        hPuMean_pp [imult] = new TH1F(Form("hPuMean_pp_%d",imult),Form("Pu Mean <reco p_{T}/gen p_{T}> %s pp %d",algn.Data(),imult),nbins,ptbins);
        MakeHist(hPuMean_pp[imult],statop);
        
        hPuArM_pp [imult] = new TH1F(Form("hPuArM_pp_%d",imult),Form("PU <reco p_{T}/gen p_{T}> %s pp %d",algn.Data(),imult),nbins,ptbins);
        MakeHist(hPuArM_pp[imult],statop);
        
        hPuSigma_pp [imult] = new TH1F(Form("hPuSigma_pp_%d",imult),Form("PU #sigma(reco p_{T}/gen p_{T}) %s pp %d",algn.Data(),imult),nbins,ptbins);
        MakeHist(hPuSigma_pp[imult],statop);    
        
        hPuRMS_pp [imult] = new TH1F(Form("hPuRMS_pp_%d",imult),Form("PU RMS(reco p_{T}/gen p_{T}) %s pp %d",algn.Data(),imult),nbins,ptbins);
        MakeHist(hPuRMS_pp[imult],statop);     
      
      hratiocorrrefpt_pp[imult]  = (TH2F*)fin_pPb->Get(Form("hrescrpt_genm%d_%d",nj,imult));
      hratiocorrrefpt_pp[imult]->SetName(Form("hratiocorrrefpt_pPb%d_%d",nj,imult));
      //hratiocorrrefpt_pp[imult]->Rebin2D(kRebin,2);

        hPuratiocorrrefpt_pp[imult]  = (TH2F*)fin_pPb->Get(Form("hrescrpt_genm%d_%d",(nj+4),imult));
        hPuratiocorrrefpt_pp[imult]->SetName(Form("hratiocorrrefpt_pPb%d_%d",(nj+4),imult));
        
  
      for(int ip=0;ip<nbins;ip++){
	int lbin = (int)(ptbins[ip]   - minval)/binw +1;
	int hbin = (int)(ptbins[ip+1] - minval)/binw +1;
	hratiocorrrefpt1D_pp[imult][ip]  = (TH1F*)hratiocorrrefpt_pp [imult]->ProjectionY(Form("hratiocorrrefpt1D_pp%d_%d_%d",nj,imult,ip),lbin,hbin);
	if(hratiocorrrefpt1D_pp[imult][ip]->GetEntries()<maxEntry)continue;
	FillMeanSigma(ip,hratiocorrrefpt1D_pp[imult][ip],hArM_pp[imult],hRMS_pp[imult],hMean_pp[imult],hSigma_pp[imult]);

          hPuratiocorrrefpt1D_pp[imult][ip]  = (TH1F*)hPuratiocorrrefpt_pp [imult]->ProjectionY(Form("hPuratiocorrrefpt1D_pp%d_%d_%d",(nj+4),imult,ip),lbin,hbin);
          if(hPuratiocorrrefpt1D_pp[imult][ip]->GetEntries()<maxEntry)continue;
          FillMeanSigma(ip,hPuratiocorrrefpt1D_pp[imult][ip],hPuArM_pp[imult],hPuRMS_pp[imult],hPuMean_pp[imult],hPuSigma_pp[imult]);
      }//! ip bin

    }//! imult loop ends
      //for HFplusEta4 bin loop
      for(int ihf=0;ihf<nhf;ihf++){
          //! pp /////////////////////////////

          hMean_pphf[ihf] = new TH1F(Form("hMean_pphf_%d",ihf),Form("Mean <reco p_{T}/gen p_{T}> %s pp hf %d",algn.Data(),ihf),nbins,ptbins);
          MakeHist(hMean_pphf[ihf],statop);
        
          hArM_pphf[ihf] = new TH1F(Form("hArM_pphf_%d",ihf),Form("<reco p_{T}/gen p_{T}> %s pp hf %d",algn.Data(),ihf),nbins,ptbins);
          MakeHist(hArM_pphf[ihf],statop);
          
          hSigma_pphf[ihf] = new TH1F(Form("hSigma_pphf_%d",ihf),Form("#sigma(reco p_{T}/gen p_{T}) %s pp hf%d",algn.Data(),ihf),nbins,ptbins);
          MakeHist(hSigma_pphf[ihf],statop);    
      
          hRMS_pphf[ihf] = new TH1F(Form("hRMS_pphf_%d",ihf),Form("RMS(reco p_{T}/gen p_{T}) %s pp hf %d",algn.Data(),ihf),nbins,ptbins);
          MakeHist(hRMS_pphf[ihf],statop); 
          
          hPuMean_pphf[ihf] = new TH1F(Form("hPuMean_pphf_%d",ihf),Form("PU Mean <reco p_{T}/gen p_{T}> %s pp hf %d",algn.Data(),ihf),nbins,ptbins);
          MakeHist(hPuMean_pphf[ihf],statop);
          
          hPuArM_pphf[ihf] = new TH1F(Form("hPuArM_pphf_%d",ihf),Form("PU <reco p_{T}/gen p_{T}> %s pp hf %d",algn.Data(),ihf),nbins,ptbins);
          MakeHist(hPuArM_pphf[ihf],statop);
          
          hPuSigma_pphf[ihf] = new TH1F(Form("hPuSigma_pphf_%d",ihf),Form("PU #sigma(reco p_{T}/gen p_{T}) %s pp hf%d",algn.Data(),ihf),nbins,ptbins);
          MakeHist(hPuSigma_pphf[ihf],statop);    
          
          hPuRMS_pphf[ihf] = new TH1F(Form("hPuRMS_pphf_%d",ihf),Form("PU RMS(reco p_{T}/gen p_{T}) %s pp hf %d",algn.Data(),ihf),nbins,ptbins);
          MakeHist(hPuRMS_pphf[ihf],statop); 
          
          hratiocorrrefpt_pphf[ihf]  = (TH2F*)fin_pPb->Get(Form("hratiocorrrefpt_genhfbin%d_%d",nj,ihf));
          hratiocorrrefpt_pphf[ihf]->SetName(Form("hratiocorrrefpt_pPbhf%d_%d",nj,ihf));
          hPuratiocorrrefpt_pphf[ihf]  = (TH2F*)fin_pPb->Get(Form("hratiocorrrefpt_genhfbin%d_%d",(nj+4),ihf));
          hPuratiocorrrefpt_pphf[ihf]->SetName(Form("hPuratiocorrrefpt_pPbhf%d_%d",(nj+4),ihf));
          for(int ip=0;ip<nbins;ip++){
              int lbin = (int)(ptbins[ip]   - minval)/binw +1;
              int hbin = (int)(ptbins[ip+1] - minval)/binw +1;
         //     cout <<"lbin=" << lbin <<"  hbin="<<hbin <<endl;
              hratiocorrrefpt1D_pphf[ihf][ip]  = (TH1F*)hratiocorrrefpt_pphf [ihf]->ProjectionY(Form("hratiocorrrefpt1D_pphf%d_%d_%d",nj,ihf,ip),lbin,hbin);
              if(hratiocorrrefpt1D_pphf[ihf][ip]->GetEntries()<maxEntry)continue;
              FillMeanSigma(ip,hratiocorrrefpt1D_pphf[ihf][ip],hArM_pphf[ihf],hRMS_pphf[ihf],hMean_pphf[ihf],hSigma_pphf[ihf]);
              hPuratiocorrrefpt1D_pphf[ihf][ip]  = (TH1F*)hPuratiocorrrefpt_pphf [ihf]->ProjectionY(Form("hPuratiocorrrefpt1D_pphf%d_%d_%d",(nj+4),ihf,ip),lbin,hbin);
              if(hPuratiocorrrefpt1D_pphf[ihf][ip]->GetEntries()<maxEntry)continue;
              FillMeanSigma(ip,hPuratiocorrrefpt1D_pphf[ihf][ip],hPuArM_pphf[ihf],hPuRMS_pphf[ihf],hPuMean_pphf[ihf],hPuSigma_pphf[ihf]);

          }
      } //! hf loop ends
      
//  }//! nj loop ends
  //return 0;

//  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<nmult;ic++){
      MakeHistMean(hArM_pp  [ic],1.108,0.904);
      MakeHistRMS (hRMS_pp  [ic],0.33,0.001);
      MakeHistMean(hMean_pp [ic],1.108,0.904);
      MakeHistRMS (hSigma_pp[ic],0.33,0.001);

      MakeHistMean(hPuArM_pp  [ic],1.108,0.904);
      MakeHistRMS (hPuRMS_pp  [ic],0.33,0.001);
      MakeHistMean(hPuMean_pp [ic],1.108,0.904);
      MakeHistRMS (hPuSigma_pp[ic],0.33,0.001);
    } 
      for(int ihf=0;ihf<nhf;ihf++){
        MakeHistMean(hArM_pphf[ihf],1.108,0.904);
      MakeHistRMS (hRMS_pphf[ihf],0.33,0.001);
      MakeHistMean(hMean_pphf[ihf],1.108,0.904);
      MakeHistRMS (hSigma_pphf[ihf],0.33,0.001);
          MakeHistMean(hPuArM_pphf[ihf],1.108,0.9044);
          MakeHistRMS (hPuRMS_pphf[ihf],0.33,0.001);
          MakeHistMean(hPuMean_pphf[ihf],1.108,0.904);
          MakeHistRMS (hPuSigma_pphf[ihf],0.33,0.001);
      }
//  }


    TLegend *t1=new TLegend(0.35,0.65,0.7,0.82);
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

    int maxc=3;
    maxr=2;
    TCanvas *c1 = new TCanvas(Form("c1"),Form("%s",algn.Data()),3,230,1465,605);
    makeMultiPanelCanvas(c1,maxc,maxr,0.0,0.0,0.22,0.22,0.02);

      TPaveText *pt3 = new TPaveText(0.3273642,0.777537,0.6165358,0.9545994,"brNDC");
      pt3->SetBorderSize(0);
      pt3->SetFillColor(10);
      pt3->SetTextFont(42);
      TText *text3 = pt3->AddText("CMS Preliminary");
      text3->SetTextSize(0.09);
      if(cmode==kPYTHIA)
          TText *text4 = pt3->AddText("PYTHIA");
      else 
          TText *text4 = pt3->AddText("PYTHIA+HIJING");
      text4->SetTextSize(0.09);

    TPaveText *pt = new TPaveText(0.2626385,0.02689651,0.5303767,0.1458786,"brNDC");
    pt->SetBorderSize(0);
    pt->SetFillColor(10);
    pt->SetTextFont(42);
      
    
    TPaveText *pt1 = new TPaveText(0.709289,0.627715,0.941172,0.9182789,"brNDC");
    pt1->SetBorderSize(0);
    pt1->SetFillColor(10);
    pt1->SetTextFont(42);
    TText *text1 ;
    TText *text2;    
    if(algo==kPF){
       text1 = pt1->AddText(Form("ak3PF")); 
        text2 = pt1->AddText(Form("akPu3PF"));
    }
    else {
        text1 = pt1->AddText(Form("ak3Calo")); 
        text2 = pt1->AddText(Form("akPu3Calo"));        
    }
    text1->SetTextSize(0.09);
    text2->SetTextSize(0.09);
    text2->SetTextColor(2);
    TText *text3 = pt1->AddText("|#eta|<3.0");
    text3->SetTextSize(0.09);

    
    TLine *line = new TLine(xmin,1,xmax,1);
    line->SetLineWidth(1);
    line->SetLineStyle(2);
   // line->Draw();
    
    switch(dis){
        case kMultJER:
            for(int ic=0;ic<nmult;ic++){ 
                c1->cd(ic+1);
                hRMS_pp[ic]->SetMarkerStyle(20);
                hRMS_pp[ic]->SetMarkerColor(1);
                hRMS_pp[ic]->SetLineColor(1);
                hRMS_pp[ic]->SetMarkerSize(1.3);
                
                hPuRMS_pp[ic]->SetMarkerStyle(24);
                hPuRMS_pp[ic]->SetMarkerColor(2);
                hPuRMS_pp[ic]->SetLineColor(2);
                hPuRMS_pp[ic]->SetMarkerSize(1.3);                
                //gPad->SetLogx();
                hRMS_pp[ic]->Draw("p");
                hPuRMS_pp[ic]->Draw("psame");
                if(ic<maxc) drawText(Form("%s",cmult[ic]),0.3,0.10,17);
                else drawText(Form("%s",cmult[ic]),0.3,0.3,17);     
                //
//                TText *text = pt->AddText(Form("%s",cmult[ic]));
//                text->SetTextSize(0.09);
//                pt->Draw();  
                if(ic==0){
                    pt3->Draw();
                    pt1->Draw();
                }
                
                
            }
            if(SavePlot){
                c1->Print(Form("%s/%sMultBin%sJetJER.pdf", plotsdir.Data(),gen.Data(), algn.Data()));  
                
            }
            c1->Update();
            break;
        case kMultJES:
            for(int ic=0;ic<nmult;ic++){ 
                c1->cd(ic+1);
                hArM_pp[ic]->SetMarkerStyle(20);
                hArM_pp[ic]->SetMarkerColor(1);
                hArM_pp[ic]->SetLineColor(1);
                hArM_pp[ic]->SetMarkerSize(1.3);
                
                hPuArM_pp[ic]->SetMarkerStyle(24);
                hPuArM_pp[ic]->SetMarkerColor(2);
                hPuArM_pp[ic]->SetLineColor(2);
                hPuArM_pp[ic]->SetMarkerSize(1.3);
                //gPad->SetLogx();
                hArM_pp[ic]->Draw("p");
                hPuArM_pp[ic]->Draw("psame");
                if(ic<maxc) drawText(Form("%s",cmult[ic]),0.3,0.10,17);
                else drawText(Form("%s",cmult[ic]),0.3,0.3,17);                
//
//                TText *text = pt->AddText(Form("%s",cmult[ic]));
//                text->SetTextSize(0.09);
//                pt->Draw();  
                if(ic==0){
                    pt3->Draw();
                    pt1->Draw();
                }
                line->Draw();  
                
            }
            if(SavePlot){
                c1->Print(Form("%s/%sMultBin%sJetJES.pdf", plotsdir.Data(),gen.Data(), algn.Data()));  
                
            }
            c1->Update();
            break ;
        case kHFJER:
            for(int ic=0;ic<nhf;ic++){ 
                c1->cd(ic+1);
                hRMS_pphf[ic]->SetMarkerStyle(20);
                hRMS_pphf[ic]->SetMarkerColor(1);
                hRMS_pphf[ic]->SetLineColor(1);
                hRMS_pphf[ic]->SetMarkerSize(1.3);
                
                hPuRMS_pphf[ic]->SetMarkerStyle(24);
                hPuRMS_pphf[ic]->SetMarkerColor(2);
                hPuRMS_pphf[ic]->SetLineColor(2);
                hPuRMS_pphf[ic]->SetMarkerSize(1.3);
                
                //gPad->SetLogx();
                hRMS_pphf[ic]->Draw("p");
                hPuRMS_pphf[ic]->Draw("psame");
                if(ic<maxc) drawText(Form("%s",chf[ic]),0.3,0.10,17);
                else drawText(Form("%s",chf[ic]),0.3,0.3,17);
                //                TText *text = pt->AddText(Form("%s",chf[ic]));
//                text->SetTextSize(0.09);
//                pt->Draw();  
                if(ic==0){
                    pt3->Draw();
                    pt1->Draw();
                }
                
  
            }
            if(SavePlot){
                c1->Print(Form("%s/%sHFplusEta%sJetJER.pdf", plotsdir.Data(),gen.Data(), algn.Data()));  
                
            }
            c1->Update();
            
            break ;
        case kHFJES:
            for(int ic=0;ic<nhf;ic++){ 
                c1->cd(ic+1);
                hArM_pphf[ic]->SetMarkerStyle(20);
                hArM_pphf[ic]->SetMarkerColor(1);
                hArM_pphf[ic]->SetLineColor(1);
                hArM_pphf[ic]->SetMarkerSize(1.3);
                
                hPuArM_pphf[ic]->SetMarkerStyle(24);
                hPuArM_pphf[ic]->SetMarkerColor(2);
                hPuArM_pphf[ic]->SetLineColor(2);
                hPuArM_pphf[ic]->SetMarkerSize(1.3);
                //gPad->SetLogx();
                hArM_pphf[ic]->Draw("p");
                hPuArM_pphf[ic]->Draw("psame");
                if(ic<maxc) drawText(Form("%s",chf[ic]),0.3,0.10,17);
                else drawText(Form("%s",chf[ic]),0.3,0.3,17);
//                TText *text = pt->AddText(Form("%s",chf[ic]));
//                text->SetTextSize(0.09);
//                pt->Draw();  
                if(ic==0){
                    pt3->Draw();
                    pt1->Draw();
                }
                line->Draw();  
            }
            if(SavePlot){
                c1->Print(Form("%sHFplusEta%sJetJES.C", gen.Data(), algn.Data()));
//                c1->Print(Form("%s/%sHFplusEta%sJetJES.pdf", plotsdir.Data(),gen.Data(), algn.Data()));
                
            }
            if(SaveFile){
                ofstream myfile;
                myfile.open(Form("akPu3PFJESNonClosure.txt"));
                for(int ic=0;ic<nhf;ic++){ 
                myfile <<"HF bins: " <<chf[ic]<<"\n" ;
                for(int ip=0;ip<nbins;ip++){                
                    cout <<"ip=" <<ip << "pt ="<<hPuArM_pphf[ic]->GetBinCenter(ip)<<"JES =" <<hPuArM_pphf[ic]->GetBinContent(ip) <<endl ; 
                    myfile <<"iHFbin = "<< ic <<", ipt ="<< hPuArM_pphf[ic]->GetBinCenter(ip) <<", ratio ="<< hPuArM_pphf[ic]->GetBinContent(ip) <<"\n" ;
                    
                }
                myfile <<"\n";
            }
            }

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
void MakeHist(TH1 *histo,int istat)
{
  histo->SetStats(istat);
  histo->SetMarkerStyle(24);
  histo->SetMarkerColor(1);
  histo->SetLineColor(1);
  histo->SetLineStyle(1);
  histo->GetXaxis()->SetTitle("p_{T}^{GenJet} (GeV/c)");
  histo->GetXaxis()->CenterTitle(true);
  histo->GetYaxis()->SetTitle("<p_{T}^{RecoJet}/p_{T}^{GenJet}>");
  histo->GetYaxis()->CenterTitle(true);
}
void FillMeanSigma(int ip,TH1 *h1F,TH1 *hArM,TH1 *hRMS,TH1 *hMean,TH1 *hSigma)
{

  TF1 *f1 = new TF1("f1","(([0]/(2*pi*[1]*[1]))*exp(-1.*((x-[2])*(x-[2])/(2*[1]*[1]))))",0,5);
  f1->SetParameters(1,0.1,1);
  f1->SetParNames("A","#sigma","mean");
  f1->SetParLimits(0,0,20);     //! A
  f1->SetParLimits(1,0.,10.0);  //! Sigma
  f1->SetParLimits(2,fitmin,fitmax); //! Mean
  f1->SetLineWidth(1);
  f1->SetNpx(knpx);
  
  float mm=0,ss=0,p0=0;  
  if(h1F->GetEntries()<maxEntry){
    h1F->Scale(0.);
    hArM  ->SetBinContent(ip,-9999);
    hArM  ->SetBinError  (ip,0);
    hRMS  ->SetBinContent(ip,-9999);
    hRMS  ->SetBinError  (ip,0);
  }
  if(h1F->Integral()>0){

    h1F->Scale(1./h1F->Integral());
    if(iFit==0){
      h1F->Fit("gaus",fopt,"",fitmin,fitmax);
      TF1* f2 = (TF1*)h1F->GetFunction("gaus");
      f2->SetLineWidth(1);
      f2->SetLineStyle(2);
      f2->SetNpx(knpx);
      hMean ->SetBinContent(ip+1,f2->GetParameter(1));
      hSigma->SetBinContent(ip+1,f2->GetParameter(2));
      
      if(strcmp(fopt,"MLRQ+")==0){
	hMean ->SetBinError  (ip+1,h1F->GetMeanError());
	hSigma->SetBinError  (ip+1,h1F->GetRMSError());
      }else{
	hMean ->SetBinError(ip+1,f2->GetParError(1));
	hSigma->SetBinError(ip+1,f2->GetParError(2));
      }
    }else{
      mm = h1F->GetMean();
      ss = h1F->GetRMS();
      p0 = h1F->GetMaximum();
      f1->SetParameters(p0,ss,mm);
      f1->SetParLimits(0,0,2*p0);
      f1->SetParLimits(1,0,2*ss);
      f1->SetParLimits(2,fitmin,fitmax);
      //f1->SetParLimits(2,mm-2.5*ss,mm+2.5*ss);
      
      h1F->Fit("f1",fopt,"",fitmin,fitmax);
      hMean ->SetBinContent(ip+1,f1->GetParameter(2));
      hSigma->SetBinContent(ip+1,f1->GetParameter(1));
      
      if(strcmp(fopt,"MLRQ+")==0){
	hMean ->SetBinError  (ip+1,h1F->GetMeanError());
	hSigma->SetBinError  (ip+1,h1F->GetRMSError());
      }else{
	hMean ->SetBinError  (ip+1,f1->GetParError(2));
	hSigma->SetBinError  (ip+1,f1->GetParError(1));
      }
    }
    hArM  ->SetBinContent(ip+1,h1F->GetMean());
    hArM  ->SetBinError  (ip+1,h1F->GetMeanError());
    hRMS  ->SetBinContent(ip+1,h1F->GetRMS());
    hRMS  ->SetBinError  (ip+1,h1F->GetRMSError());
  }
}
void MakeHistRMS(TH1 *h1,float ymax,float ymin)
{

  h1->SetTitle("");
  h1->SetMaximum(ymax);
  h1->SetMinimum(ymin);
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetNdivisions(507);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetLabelOffset(0.01);
  h1->GetXaxis()->SetLabelSize(0.07);
  h1->GetXaxis()->SetTitleSize(0.07);
  h1->GetXaxis()->SetTitleOffset(1.15);
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetTitle("#sigma (RecoJet p_{T} / GenJet p_{T})");
  h1->GetYaxis()->CenterTitle(true);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetLabelFont(42);
  h1->GetYaxis()->SetLabelOffset(0.01);
  h1->GetYaxis()->SetLabelSize(0.09);
  h1->GetYaxis()->SetTitleSize(0.09);
  h1->GetYaxis()->SetTitleOffset(1.12);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetDecimals(true);

}
void MakeHistMean(TH1 *h1,float ymax,float ymin)
{
  h1->SetMaximum(ymax);
  h1->SetMinimum(ymin);
  h1->SetTitle("");
  h1->GetXaxis()->SetRangeUser(xmin,xmax);
  h1->GetXaxis()->SetTitle("GenJet p_{T} (GeV/c)");
  h1->GetXaxis()->CenterTitle(true);
  h1->GetXaxis()->SetMoreLogLabels();
  h1->GetXaxis()->SetNoExponent();
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetXaxis()->SetLabelFont(42);
  h1->GetXaxis()->SetTitleSize(0.07);
  h1->GetXaxis()->SetTitleOffset(1.15);
  h1->GetXaxis()->SetLabelSize(0.07);
  h1->GetXaxis()->SetLabelOffset(0.005);
  h1->GetXaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetTitle("<RecoJet p_{T} / GenJet p_{T}>");
  h1->GetYaxis()->SetTitleSize(0.07);
  h1->GetYaxis()->SetTitleOffset(1.50);
  h1->GetYaxis()->SetLabelSize(0.07);
  h1->GetYaxis()->SetNdivisions(507);
  h1->GetYaxis()->SetDecimals(true);
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetLabelFont(42);
}
int GetPtBin(float pt)
{
  int ibin=-1;
  ibin = (int)(pt + minval)/binw;
  return ibin;
}
void CleanHist(TH1 *h1F,float lxval, float hxval)
{
  for(int ix=1;ix<=h1F->GetNbinsX();ix++){
    double val = h1F->GetBinCenter(ix);
    if(val<lxval){
      h1F->SetBinContent(ix,0);
      h1F->SetBinError(ix,0);
    }
    else  if(val>hxval){
      h1F->SetBinContent(ix,0);
      h1F->SetBinError(ix,0);
    }
  }
}
