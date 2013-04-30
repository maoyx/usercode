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
const char *calgo[knj] = {"ak2PF","ak3PF","ak4PF","ak5PF",
    "akPu2PF","akPu3PF","akPu4PF","akPu5PF",
    "ak2Calo","ak3Calo","ak4Calo","ak5Calo",
    "akPu2Calo","akPu3Calo","akPu4Calo","akPu5Calo"};
const char *algn[knj] = {"ak2PF","ak3PF","ak4PF","ak5PF",
    "akPu2PF","akPu3PF","akPu4PF","akPu5PF",
    "ak2Calo","ak3Calo","ak4Calo","ak5Calo",
    "akPu2Calo","akPu3Calo","akPu4Calo","akPu5Calo"};

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
const char *chf[nhf] = {"E_{T}^{HF(#eta=4)}<5",
    "5<E_{T}^{HF(#eta=4)}<10",
    "10<E_{T}^{HF(#eta=4)}<15",
    "15<E_{T}^{HF(#eta=4)}<20",
    "20<E_{T}^{HF(#eta=4)}<30",
    "30<E_{T}^{HF(#eta=4)}"
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

int Mult_JetEnergyResolution(const char *cmode="hijing",int rfit=0)
{
  int statop=0;
  int kRebin=1;
  const char *reta = "eta2.0";
  //bool iSave=false;
  bool iMult=false;

    TString plotsdir = "/Users/ymao/group/CMS/plots/pA";
    bool SavePlot = kFALSE ;
    TString gen ;
    if(cmode=="pythia") gen="PYTHIA";
    else gen="PYTHIAHIJING";    

  float ketacut=1.6;
  if(strcmp(reta,"eta2.0")==0)ketacut=2;
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
    if(cmode=="pythia")
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
  TH2F *hratiocorrrefpt_pp[knj][nmult];
  TH1F *hratiocorrrefpt1D_pp[knj][nmult][nbins];  
  TH1F *hMean_pp [knj][nmult], *hArM_pp[knj][nmult], *hSigma_pp [knj][nmult], *hRMS_pp[knj][nmult];

    TH2F *hratiocorrrefpt_pphf[knj][nhf];
    TH1F *hratiocorrrefpt1D_pphf[knj][nhf][nbins];  
    TH1F *hMean_pphf[knj][nhf], *hArM_pphf[knj][nhf], *hSigma_pphf[knj][nhf], *hRMS_pphf[knj][nhf];
 
  TH2F *hratiorawrefpt_pp[knj][nmult];
  TH1F *hratiorawrefpt1D_pp[knj][nmult][nbins];  
  TH1F *hMean_r_pp [knj][nmult], *hArM_r_pp[knj][nmult], *hSigma_r_pp [knj][nmult], *hRMS_r_pp[knj][nmult];

  TProfile *hrecoraw[knj][nmult];


  for(int nj=0;nj<knj;nj++){
    //cout<<"nj : "<<nj<<Form("\t %s",calgo[nj])<<endl;
    
    
    for(int imult=0;imult<nmult;imult++){
      //! pp /////////////////////////////
      hMean_pp [nj][imult] = new TH1F(Form("hMean_pp%d_%d",nj,imult),Form("<reco p_{T}/gen p_{T}> %s pp %d",calgo[nj],imult),nbins,ptbins);
      MakeHist(hMean_pp[nj][imult],statop);

      hArM_pp [nj][imult] = new TH1F(Form("hArM_pp%d_%d",nj,imult),Form("<reco p_{T}/gen p_{T}> %s pp %d",calgo[nj],imult),nbins,ptbins);
      MakeHist(hArM_pp[nj][imult],statop);

      hSigma_pp [nj][imult] = new TH1F(Form("hSigma_pp%d_%d",nj,imult),Form("#sigma(reco p_{T}/gen p_{T}) %s pp %d",calgo[nj],imult),nbins,ptbins);
      MakeHist(hSigma_pp[nj][imult],statop);    

      hRMS_pp [nj][imult] = new TH1F(Form("hRMS_pp%d_%d",nj,imult),Form("RMS(reco p_{T}/gen p_{T}) %s pp %d",calgo[nj],imult),nbins,ptbins);
      MakeHist(hRMS_pp[nj][imult],statop);        

      hMean_r_pp [nj][imult] = new TH1F(Form("hMean_r_pp%d_%d",nj,imult),Form("<raw p_{T}/gen p_{T}> %s pp %d",calgo[nj],imult),nbins,ptbins);
      MakeHist(hMean_r_pp[nj][imult],statop);

      hArM_r_pp [nj][imult] = new TH1F(Form("hArM_r_pp%d_%d",nj,imult),Form("<raw p_{T}/gen p_{T}> %s pp %d",calgo[nj],imult),nbins,ptbins);
      MakeHist(hArM_r_pp[nj][imult],statop);

      hSigma_r_pp [nj][imult] = new TH1F(Form("hSigma_r_pp%d_%d",nj,imult),Form("#sigma(raw p_{T}/gen p_{T}) %s pp %d",calgo[nj],imult),nbins,ptbins);
      MakeHist(hSigma_r_pp[nj][imult],statop);    

      hRMS_r_pp [nj][imult] = new TH1F(Form("hRMS_r_pp%d_%d",nj,imult),Form("RMS(raw p_{T}/gen p_{T}) %s pp %d",calgo[nj],imult),nbins,ptbins);
      MakeHist(hRMS_r_pp[nj][imult],statop);        
      
      hratiocorrrefpt_pp[nj][imult]  = (TH2F*)fin_pPb->Get(Form("hrescrpt_genm%d_%d",nj,imult));
      hratiocorrrefpt_pp[nj][imult]->SetName(Form("hratiocorrrefpt_pPb%d_%d",nj,imult));
      //hratiocorrrefpt_pp[nj][imult]->Rebin2D(kRebin,2);

      hratiorawrefpt_pp[nj][imult]  = (TH2F*)fin_pPb->Get(Form("hresrrpt_genm%d_%d",nj,imult));
      hratiorawrefpt_pp[nj][imult]->SetName(Form("hratiorawrefpt_pPb%d_%d",nj,imult));
      //hratiorawrefpt_pp[nj][imult]->Rebin2D(kRebin,2);

      hrecoraw[nj][imult] = (TProfile*)fin_pPb->Get(Form("hrecoraw%d_%d",nj,imult));
  
      for(int ip=0;ip<nbins;ip++){
	int lbin = (int)(ptbins[ip]   - minval)/binw +1;
	int hbin = (int)(ptbins[ip+1] - minval)/binw +1;
	hratiocorrrefpt1D_pp[nj][imult][ip]  = (TH1F*)hratiocorrrefpt_pp [nj][imult]->ProjectionY(Form("hratiocorrrefpt1D_pp%d_%d_%d",nj,imult,ip),lbin,hbin);
	if(hratiocorrrefpt1D_pp[nj][imult][ip]->GetEntries()<maxEntry)continue;
	FillMeanSigma(ip,hratiocorrrefpt1D_pp[nj][imult][ip],hArM_pp[nj][imult],hRMS_pp[nj][imult],hMean_pp[nj][imult],hSigma_pp[nj][imult]);

	//! Raw/Gen
	hratiorawrefpt1D_pp[nj][imult][ip]  = (TH1F*)hratiorawrefpt_pp [nj][imult]->ProjectionY(Form("hratiorawrefpt1D_pp%d_%d_%d",nj,imult,ip),lbin,hbin);
	if(hratiorawrefpt1D_pp[nj][imult][ip]->GetEntries()<maxEntry)continue;
	FillMeanSigma(ip,hratiorawrefpt1D_pp[nj][imult][ip],hArM_r_pp[nj][imult],hRMS_r_pp[nj][imult],hMean_r_pp[nj][imult],hSigma_r_pp[nj][imult]);
      }//! ip bin

    }//! imult loop ends
      //for HFplusEta4 bin loop
      for(int ihf=0;ihf<nhf;ihf++){
          //! pp /////////////////////////////

          hMean_pphf[nj][ihf] = new TH1F(Form("hMean_pphf%d_%d",nj,ihf),Form("<reco p_{T}/gen p_{T}> %s pp hf %d",calgo[nj],ihf),nbins,ptbins);
          MakeHist(hMean_pphf[nj][ihf],statop);
        
          hArM_pphf[nj][ihf] = new TH1F(Form("hArM_pphf%d_%d",nj,ihf),Form("<reco p_{T}/gen p_{T}> %s pp hf %d",calgo[nj],ihf),nbins,ptbins);
          MakeHist(hArM_pphf[nj][ihf],statop);
          
          hSigma_pphf[nj][ihf] = new TH1F(Form("hSigma_pphf%d_%d",nj,ihf),Form("#sigma(reco p_{T}/gen p_{T}) %s pp hf%d",calgo[nj],ihf),nbins,ptbins);
          MakeHist(hSigma_pphf[nj][ihf],statop);    
      
          hRMS_pphf[nj][ihf] = new TH1F(Form("hRMS_pphf%d_%d",nj,ihf),Form("RMS(reco p_{T}/gen p_{T}) %s pp hf %d",calgo[nj],ihf),nbins,ptbins);
          MakeHist(hRMS_pphf[nj][ihf],statop); 
        //  cout<<"nj : "<<nj<<Form("\t %s",calgo[nj])<<endl;

       //   cout <<"ihf =" <<ihf<<endl;
          hratiocorrrefpt_pphf[nj][ihf]  = (TH2F*)fin_pPb->Get(Form("hratiocorrrefpt_genhfbin%d_%d",nj,ihf));
          hratiocorrrefpt_pphf[nj][ihf]->SetName(Form("hratiocorrrefpt_pPbhf%d_%d",nj,ihf));
          for(int ip=0;ip<nbins;ip++){
              int lbin = (int)(ptbins[ip]   - minval)/binw +1;
              int hbin = (int)(ptbins[ip+1] - minval)/binw +1;
         //     cout <<"lbin=" << lbin <<"  hbin="<<hbin <<endl;
              hratiocorrrefpt1D_pphf[nj][ihf][ip]  = (TH1F*)hratiocorrrefpt_pphf [nj][ihf]->ProjectionY(Form("hratiocorrrefpt1D_pphf%d_%d_%d",nj,ihf,ip),lbin,hbin);
              if(hratiocorrrefpt1D_pphf[nj][ihf][ip]->GetEntries()<maxEntry)continue;
              FillMeanSigma(ip,hratiocorrrefpt1D_pphf[nj][ihf][ip],hArM_pphf[nj][ihf],hRMS_pphf[nj][ihf],hMean_pphf[nj][ihf],hSigma_pphf[nj][ihf]);
          }
      } //! hf loop ends
      
  }//! nj loop ends
  //return 0;

  for(int nj=0;nj<knj;nj++){
    for(int ic=0;ic<nmult;ic++){
      MakeHistMean(hArM_pp  [nj][ic],1.088,0.84);
      MakeHistRMS (hRMS_pp  [nj][ic],0.43,0.0001);
      MakeHistMean(hMean_pp [nj][ic],1.088,0.84);
      MakeHistRMS (hSigma_pp[nj][ic],0.43,0.0001);

      MakeHistMean(hArM_r_pp  [nj][ic],1.088,0.84);
      MakeHistRMS (hRMS_r_pp  [nj][ic],0.43,0.0001);
      MakeHistMean(hMean_r_pp [nj][ic],1.088,0.84);
      MakeHistRMS (hSigma_r_pp[nj][ic],0.43,0.0001);
    } 
      for(int ihf=0;ihf<nhf;ihf++){
        MakeHistMean(hArM_pphf[nj][ihf],1.088,0.84);
      MakeHistRMS (hRMS_pphf[nj][ihf],0.43,0.0001);
      MakeHistMean(hMean_pphf[nj][ihf],1.088,0.84);
      MakeHistRMS (hSigma_pphf[nj][ihf],0.43,0.0001);
      }
  }


  TLegend *l3=0;
  l3 = new TLegend(0.2453033,0.6247655,0.4838159,0.9399625,NULL,"BRNDC");
  l3->SetHeader("");
  l3->SetBorderSize(0);
  l3->SetTextFont(42);
  l3->SetTextSize(0.09);
  l3->SetLineColor(1);
  l3->SetLineStyle(1);
  l3->SetLineWidth(1);
  l3->SetFillColor(10);
  l3->SetFillStyle(1001);
  l3->SetHeader("");			  
  l3->AddEntry(hRMS_pp[0][0],"ak5Calo","p"); 
  l3->AddEntry(hRMS_pp[1][0],"ak5PF","p"); 

      TPaveText *pt3 = new TPaveText(0.3273642,0.777537,0.6165358,0.9545994,"brNDC");
      pt3->SetBorderSize(0);
      pt3->SetFillColor(10);
      pt3->SetTextFont(42);
      TText *text3 = pt3->AddText("CMS Preliminary");
      text3->SetTextSize(0.09);
      if(cmode=="pythia")
          TText *text4 = pt3->AddText("PYTHIA");
      else 
          TText *text4 = pt3->AddText("PYTHIA+HIJING");
      text4->SetTextSize(0.09);

      int maxc=4;
      maxr=2;
      
      TCanvas *c1[nhf];
      
      for(int ic=0;ic<nhf;ic++){
          c1[ic] = new TCanvas(Form("c1_%d",ic),Form("%s PF JES JER",chf[ic]),3,230,1465,605);
          makeMultiPanelCanvas(c1[ic],maxc,maxr,0.0,0.0,0.22,0.22,0.02);
          ipad=0;
          for(int nj=0;nj<4;nj++){
              
              hRMS_pphf[nj][ic]->SetMarkerStyle(20);
              hRMS_pphf[nj][ic]->SetMarkerColor(1);
              hRMS_pphf[nj][ic]->SetLineColor(1);
              hRMS_pphf[nj][ic]->SetMarkerSize(1.3);
              
              hRMS_pphf[nj+4][ic]->SetMarkerStyle(24);
              hRMS_pphf[nj+4][ic]->SetMarkerColor(2);
              hRMS_pphf[nj+4][ic]->SetLineColor(2);
              hRMS_pphf[nj+4][ic]->SetMarkerSize(1.3);
              
              hArM_pphf[nj][ic]->SetMarkerStyle(20);
              hArM_pphf[nj][ic]->SetMarkerColor(1);
              hArM_pphf[nj][ic]->SetLineColor(1);
              hArM_pphf[nj][ic]->SetMarkerSize(1.3);
              
              hArM_pphf[nj+4][ic]->SetMarkerStyle(24);
              hArM_pphf[nj+4][ic]->SetMarkerColor(2);
              hArM_pphf[nj+4][ic]->SetLineColor(2);
              hArM_pphf[nj+4][ic]->SetMarkerSize(1.3);
              c1[ic]->cd(++ipad);
              //gPad->SetLogx();
              hRMS_pphf[nj][ic]->Draw("p");
              hRMS_pphf[nj+4][ic]->Draw("psame");
              
              if(ipad==1){
                  pt3->Draw();
              }
              
              if(ipad==1){
                  TPaveText *pt = new TPaveText(0.2626385,0.02689651,0.5303767,0.1458786,"brNDC");
                  pt->SetBorderSize(0);
                  pt->SetFillColor(10);
                  pt->SetTextFont(42);
                  TText *text = pt->AddText(Form("%s",chf[ic]));
                  text->SetTextSize(0.09);
                  pt->Draw();
              }
              
              TPaveText *pt1 = new TPaveText(0.709289,0.627715,0.941172,0.9182789,"brNDC");
              pt1->SetBorderSize(0);
              pt1->SetFillColor(10);
              pt1->SetTextFont(42);
              TText *text1 = pt1->AddText(Form("%s",algn[nj]));
              text1->SetTextSize(0.09);
              TText *text2 = pt1->AddText(Form("%s",algn[nj+4]));
              text2->SetTextSize(0.09);
              text2->SetTextColor(2);
              TText *text3 = pt1->AddText("|#eta|<2.0");
              text3->SetTextSize(0.09);
              pt1->Draw();
              
              
              c1[ic]->cd(ipad+maxc);
              
              hArM_pphf[nj][ic]->Draw("p");
              hArM_pphf[nj+4][ic]->Draw("psame");
              //hArM_r_pp[nj][ic]->Draw("psame");      
              
              if(ipad==5){
                  TLatex *tex1 = new TLatex(175,0.68,"Raw");
                  tex1->SetTextFont(42);
                  tex1->SetTextSize(0.08);
                  tex1->SetLineWidth(2);
                  tex1->Draw();
                  TLatex *tex2 = new TLatex(174,0.63,"Corrected");
                  tex2->SetTextFont(42);
                  tex2->SetTextSize(0.08);
                  tex2->SetLineWidth(2);
                  tex2->Draw();
                  TMarker *marker1 = new TMarker(162,0.65,25);
                  marker1->SetMarkerStyle(25);
                  marker1->Draw();
                  TMarker *marker2 = new TMarker(139,0.65,21);
                  marker2->SetMarkerStyle(21);
                  marker2->Draw();
                  TMarker *marker3 = new TMarker(162,0.70,24);
                  marker3->SetMarkerStyle(24);
                  marker3->Draw();
                  TMarker *marker4 = new TMarker(139,0.70,20);
                  marker4->SetMarkerStyle(20);
                  marker4->Draw();
              }
              TLine *line = new TLine(xmin,1,xmax,1);
              line->SetLineWidth(1);
              line->SetLineStyle(2);
              line->Draw();
              
          }//! algo
          if(SavePlot){
              c1[ic]->Print(Form("%s/%sHFplusEta4%.f_%.fPFJetJESJER.pdf", plotsdir.Data(),gen.Data(), hfBin[ic], hfBin[ic+1]));  
              
          }
          c1[ic]->Update();
      }
          TCanvas *c2[nhf];
          for(int ic=0;ic<nhf;ic++){
              c2[ic] = new TCanvas(Form("c2_%d",ic),Form("%s Calo JES JER",chf[ic]),3,230,1465,605);
              makeMultiPanelCanvas(c2[ic],maxc,maxr,0.0,0.0,0.22,0.22,0.02);
              ipad=0;
              for(int nj=8;nj<12;nj++){
                  
                  hRMS_pphf[nj][ic]->SetMarkerStyle(20);
                  hRMS_pphf[nj][ic]->SetMarkerColor(1);
                  hRMS_pphf[nj][ic]->SetLineColor(1);
                  hRMS_pphf[nj][ic]->SetMarkerSize(1.3);
                  
                  hRMS_pphf[nj+4][ic]->SetMarkerStyle(24);
                  hRMS_pphf[nj+4][ic]->SetMarkerColor(2);
                  hRMS_pphf[nj+4][ic]->SetLineColor(2);
                  hRMS_pphf[nj+4][ic]->SetMarkerSize(1.3);
                  
                  hArM_pphf[nj][ic]->SetMarkerStyle(20);
                  hArM_pphf[nj][ic]->SetMarkerColor(1);
                  hArM_pphf[nj][ic]->SetLineColor(1);
                  hArM_pphf[nj][ic]->SetMarkerSize(1.3);
                  
                  hArM_pphf[nj+4][ic]->SetMarkerStyle(24);
                  hArM_pphf[nj+4][ic]->SetMarkerColor(2);
                  hArM_pphf[nj+4][ic]->SetLineColor(2);
                  hArM_pphf[nj+4][ic]->SetMarkerSize(1.3);
                  
                  c2[ic]->cd(++ipad);
                  //gPad->SetLogx();
                  hRMS_pphf[nj][ic]->Draw("p");
                  hRMS_pphf[nj+4][ic]->Draw("psame");
                  
                  if(ipad==1){
//                      TPaveText *pt3 = new TPaveText(0.3273642,0.777537,0.6165358,0.9545994,"brNDC");
//                      pt3->SetBorderSize(0);
//                      pt3->SetFillColor(10);
//                      pt3->SetTextFont(42);
//                      TText *text3 = pt3->AddText("CMS Preliminary");
//                      text3->SetTextSize(0.09);
//                      if(cmode=="pythia") TText *text4 = pt3->AddText("PYTHIA");
//                      else
//                          TText *text4 = pt3->AddText("PYTHIA+HIJING");
//                      text4->SetTextSize(0.09);
                      pt3->Draw();
                  }
                  
                  if(ipad==1){
                      TPaveText *pt = new TPaveText(0.2626385,0.02689651,0.5303767,0.1458786,"brNDC");
                      pt->SetBorderSize(0);
                      pt->SetFillColor(10);
                      pt->SetTextFont(42);
                      TText *text = pt->AddText(Form("%s",chf[ic]));
                      text->SetTextSize(0.09);
                      pt->Draw();
                  }
                  TPaveText *pt1 = new TPaveText(0.709289,0.627715,0.941172,0.9182789,"brNDC");
                  pt1->SetBorderSize(0);
                  pt1->SetFillColor(10);
                  pt1->SetTextFont(42);
                  TText *text1 = pt1->AddText(Form("%s",algn[nj]));
                  text1->SetTextSize(0.09);
                  TText *text2 = pt1->AddText(Form("%s",algn[nj+4]));
                  text2->SetTextSize(0.09);
                  text2->SetTextColor(2);
                  TText *text3 = pt1->AddText("|#eta|<2.0");
                  text3->SetTextSize(0.09);
                  pt1->Draw();
                  
                  c2[ic]->cd(ipad+maxc);
                  
                  hArM_pphf[nj][ic]->Draw("p");
                  hArM_pphf[nj+4][ic]->Draw("psame");
                  //hArM_r_pp[nj][ic]->Draw("psame");      
                  
                  if(ipad==5){
                      TLatex *tex1 = new TLatex(175,0.68,"Raw");
                      tex1->SetTextFont(42);
                      tex1->SetTextSize(0.08);
                      tex1->SetLineWidth(2);
                      tex1->Draw();
                      TLatex *tex2 = new TLatex(174,0.63,"Corrected");
                      tex2->SetTextFont(42);
                      tex2->SetTextSize(0.08);
                      tex2->SetLineWidth(2);
                      tex2->Draw();
                      TMarker *marker1 = new TMarker(162,0.65,25);
                      marker1->SetMarkerStyle(25);
                      marker1->Draw();
                      TMarker *marker2 = new TMarker(139,0.65,21);
                      marker2->SetMarkerStyle(21);
                      marker2->Draw();
                      TMarker *marker3 = new TMarker(162,0.70,24);
                      marker3->SetMarkerStyle(24);
                      marker3->Draw();
                      TMarker *marker4 = new TMarker(139,0.70,20);
                      marker4->SetMarkerStyle(20);
                      marker4->Draw();
                  }
                  TLine *line = new TLine(xmin,1,xmax,1);
                  line->SetLineWidth(1);
                  line->SetLineStyle(2);
                  line->Draw();
                  
              }//! algo
              if(SavePlot){
                  c2[ic]->Print(Form("%s/%sHFplusEta4%.f_%.fCaloJetJESJER.pdf", plotsdir.Data(),gen.Data(), hfBin[ic], hfBin[ic+1]));  
                  
              }
              c2[ic]->Update();
          }
  TCanvas *c3[nmult];

  for(int ic=0;ic<nmult;ic++){
    c3[ic] = new TCanvas(Form("c3_%d",ic),Form("%s JES JER",cmult[ic]),3,230,1465,605);
    makeMultiPanelCanvas(c3[ic],maxc,maxr,0.0,0.0,0.22,0.22,0.02);
    ipad=0;
    for(int nj=0;nj<4;nj++){
	
	hRMS_pp[nj][ic]->SetMarkerStyle(20);
	hRMS_pp[nj][ic]->SetMarkerColor(1);
	hRMS_pp[nj][ic]->SetLineColor(1);
	hRMS_pp[nj][ic]->SetMarkerSize(1.3);
      
	hRMS_pp[nj+4][ic]->SetMarkerStyle(24);
	hRMS_pp[nj+4][ic]->SetMarkerColor(2);
	hRMS_pp[nj+4][ic]->SetLineColor(2);
	hRMS_pp[nj+4][ic]->SetMarkerSize(1.3);
	
	hArM_pp[nj][ic]->SetMarkerStyle(20);
	hArM_pp[nj][ic]->SetMarkerColor(1);
	hArM_pp[nj][ic]->SetLineColor(1);
	hArM_pp[nj][ic]->SetMarkerSize(1.3);
	
	hArM_pp[nj+4][ic]->SetMarkerStyle(24);
	hArM_pp[nj+4][ic]->SetMarkerColor(2);
	hArM_pp[nj+4][ic]->SetLineColor(2);
	hArM_pp[nj+4][ic]->SetMarkerSize(1.3);
	
	hArM_r_pp[nj][ic]->SetMarkerStyle(21);
	hArM_r_pp[nj][ic]->SetMarkerColor(1);
	hArM_r_pp[nj][ic]->SetLineColor(1);
	hArM_r_pp[nj][ic]->SetMarkerSize(1.3);

	c3[ic]->cd(++ipad);
	//gPad->SetLogx();
	hRMS_pp[nj][ic]->Draw("p");
	hRMS_pp[nj+4][ic]->Draw("psame");

	if(ipad==1){
//	  TPaveText *pt3 = new TPaveText(0.3273642,0.777537,0.6165358,0.9545994,"brNDC");
//	  pt3->SetBorderSize(0);
//	  pt3->SetFillColor(10);
//	  pt3->SetTextFont(42);
//	  TText *text3 = pt3->AddText("CMS Preliminary");
//	  text3->SetTextSize(0.09);
//        if(cmode=="pythia")
//            TText *text4 = pt3->AddText("PYTHIA");
//        else 
//            TText *text4 = pt3->AddText("PYTHIA+HIJING");
//	  text4->SetTextSize(0.09);
	  pt3->Draw();
	}

	if(ipad==1){
	  TPaveText *pt = new TPaveText(0.2626385,0.02689651,0.5303767,0.1458786,"brNDC");
	  pt->SetBorderSize(0);
	  pt->SetFillColor(10);
	  pt->SetTextFont(42);
	  TText *text = pt->AddText(Form("%s",cmult[ic]));
	  text->SetTextSize(0.09);
	  pt->Draw();
	}
	
	TPaveText *pt1 = new TPaveText(0.709289,0.627715,0.941172,0.9182789,"brNDC");
	pt1->SetBorderSize(0);
	pt1->SetFillColor(10);
	pt1->SetTextFont(42);
	TText *text1 = pt1->AddText(Form("%s",algn[nj]));
	text1->SetTextSize(0.09);
	TText *text2 = pt1->AddText(Form("%s",algn[nj+4]));
	text2->SetTextSize(0.09);
	text2->SetTextColor(2);
	TText *text3 = pt1->AddText("|#eta|<2.0");
	text3->SetTextSize(0.09);
	pt1->Draw();
	
	
	c3[ic]->cd(ipad+maxc);

	hArM_pp[nj][ic]->Draw("p");
	hArM_pp[nj+4][ic]->Draw("psame");
	//hArM_r_pp[nj][ic]->Draw("psame");      

	if(ipad==5){
	  TLatex *tex1 = new TLatex(175,0.68,"Raw");
	  tex1->SetTextFont(42);
	  tex1->SetTextSize(0.08);
	  tex1->SetLineWidth(2);
	  tex1->Draw();
	  TLatex *tex2 = new TLatex(174,0.63,"Corrected");
	  tex2->SetTextFont(42);
	  tex2->SetTextSize(0.08);
	  tex2->SetLineWidth(2);
	  tex2->Draw();
	  TMarker *marker1 = new TMarker(162,0.65,25);
	  marker1->SetMarkerStyle(25);
	  marker1->Draw();
	  TMarker *marker2 = new TMarker(139,0.65,21);
	  marker2->SetMarkerStyle(21);
	  marker2->Draw();
	  TMarker *marker3 = new TMarker(162,0.70,24);
	  marker3->SetMarkerStyle(24);
	  marker3->Draw();
	  TMarker *marker4 = new TMarker(139,0.70,20);
	  marker4->SetMarkerStyle(20);
	  marker4->Draw();
	}
      TLine *line = new TLine(xmin,1,xmax,1);
      line->SetLineWidth(1);
      line->SetLineStyle(2);
      line->Draw();
      
    }//! algo
      if(SavePlot){
          c3[ic]->Print(Form("%s/%sMult%d_%dPFJetJESJER.pdf", plotsdir.Data(),gen.Data(), multBin[ic], multBin[ic+1]));  

      }
      c3[ic]->Update();
    /*
      if(iSave){
      c3[ic]->SaveAs(Form("AN/JERJES/JERJES_%s.png",cmult[ic]));
      c3[ic]->SaveAs(Form("AN/JERJES/JERJES_%s.pdf",cmult[ic]));
      c3[ic]->SaveAs(Form("AN/JERJES/JERJES_%s.C",cmult[ic));
      c3[ic]->SaveAs(Form("AN/JERJES/JERJES_%s.eps",cmult[ic]));
      }
    */
  }//! imult
  //return 0;

  TCanvas *c4[nmult];
  for(int ic=0;ic<nmult;ic++){
    c4[ic] = new TCanvas(Form("c4_%d",ic),Form("%s Calo JES JER",cmult[ic]),3,230,1465,605);
    makeMultiPanelCanvas(c4[ic],maxc,maxr,0.0,0.0,0.22,0.22,0.02);
    ipad=0;
    for(int nj=8;nj<12;nj++){

      hRMS_pp[nj][ic]->SetMarkerStyle(20);
      hRMS_pp[nj][ic]->SetMarkerColor(1);
      hRMS_pp[nj][ic]->SetLineColor(1);
      hRMS_pp[nj][ic]->SetMarkerSize(1.3);

      hRMS_pp[nj+4][ic]->SetMarkerStyle(24);
      hRMS_pp[nj+4][ic]->SetMarkerColor(2);
      hRMS_pp[nj+4][ic]->SetLineColor(2);
      hRMS_pp[nj+4][ic]->SetMarkerSize(1.3);

      hArM_pp[nj][ic]->SetMarkerStyle(20);
      hArM_pp[nj][ic]->SetMarkerColor(1);
      hArM_pp[nj][ic]->SetLineColor(1);
      hArM_pp[nj][ic]->SetMarkerSize(1.3);

      hArM_pp[nj+4][ic]->SetMarkerStyle(24);
      hArM_pp[nj+4][ic]->SetMarkerColor(2);
      hArM_pp[nj+4][ic]->SetLineColor(2);
      hArM_pp[nj+4][ic]->SetMarkerSize(1.3);

      hArM_r_pp[nj][ic]->SetMarkerStyle(21);
      hArM_r_pp[nj][ic]->SetMarkerColor(1);
      hArM_r_pp[nj][ic]->SetLineColor(1);
      hArM_r_pp[nj][ic]->SetMarkerSize(1.3);

      c4[ic]->cd(++ipad);
      //gPad->SetLogx();
      hRMS_pp[nj][ic]->Draw("p");
      hRMS_pp[nj+4][ic]->Draw("psame");

      if(ipad==1){
	TPaveText *pt3 = new TPaveText(0.3273642,0.777537,0.6165358,0.9545994,"brNDC");
	pt3->SetBorderSize(0);
	pt3->SetFillColor(10);
	pt3->SetTextFont(42);
	TText *text3 = pt3->AddText("CMS Preliminary");
	text3->SetTextSize(0.09);
          if(cmode=="pythia") TText *text4 = pt3->AddText("PYTHIA");
          else
              TText *text4 = pt3->AddText("PYTHIA+HIJING");
	text4->SetTextSize(0.09);
	pt3->Draw();
      }

      if(ipad==1){
	TPaveText *pt = new TPaveText(0.2626385,0.02689651,0.5303767,0.1458786,"brNDC");
	pt->SetBorderSize(0);
	pt->SetFillColor(10);
	pt->SetTextFont(42);
	TText *text = pt->AddText(Form("%s",cmult[ic]));
	text->SetTextSize(0.09);
	pt->Draw();
      }
      TPaveText *pt1 = new TPaveText(0.709289,0.627715,0.941172,0.9182789,"brNDC");
      pt1->SetBorderSize(0);
      pt1->SetFillColor(10);
      pt1->SetTextFont(42);
      TText *text1 = pt1->AddText(Form("%s",algn[nj]));
      text1->SetTextSize(0.09);
      TText *text2 = pt1->AddText(Form("%s",algn[nj+4]));
      text2->SetTextSize(0.09);
      text2->SetTextColor(2);
      TText *text3 = pt1->AddText("|#eta|<2.0");
      text3->SetTextSize(0.09);
      pt1->Draw();

      c4[ic]->cd(ipad+maxc);
      
      hArM_pp[nj][ic]->Draw("p");
      hArM_pp[nj+4][ic]->Draw("psame");
      //hArM_r_pp[nj][ic]->Draw("psame");      

      if(ipad==5){
	TLatex *tex1 = new TLatex(175,0.68,"Raw");
	tex1->SetTextFont(42);
	tex1->SetTextSize(0.08);
	tex1->SetLineWidth(2);
	tex1->Draw();
	TLatex *tex2 = new TLatex(174,0.63,"Corrected");
	tex2->SetTextFont(42);
	tex2->SetTextSize(0.08);
	tex2->SetLineWidth(2);
	tex2->Draw();
	TMarker *marker1 = new TMarker(162,0.65,25);
	marker1->SetMarkerStyle(25);
	marker1->Draw();
	TMarker *marker2 = new TMarker(139,0.65,21);
	marker2->SetMarkerStyle(21);
	marker2->Draw();
	TMarker *marker3 = new TMarker(162,0.70,24);
	marker3->SetMarkerStyle(24);
	marker3->Draw();
	TMarker *marker4 = new TMarker(139,0.70,20);
	marker4->SetMarkerStyle(20);
	marker4->Draw();
      }
      TLine *line = new TLine(xmin,1,xmax,1);
      line->SetLineWidth(1);
      line->SetLineStyle(2);
      line->Draw();

    }//! algo
      if(SavePlot){
          c4[ic]->Print(Form("%s/%sMult%d_%dCaloJetJESJER.pdf", plotsdir.Data(),gen.Data(), multBin[ic], multBin[ic+1]));  
          
      }
      c4[ic]->Update();
    /*
    if(iSave){
      c3[ic]->SaveAs(Form("AN/JERJES/JERJES_%s.png",cmult[ic]));
      c3[ic]->SaveAs(Form("AN/JERJES/JERJES_%s.pdf",cmult[ic]));
      c3[ic]->SaveAs(Form("AN/JERJES/JERJES_%s.C"  ,cmult[ic]));
      c3[ic]->SaveAs(Form("AN/JERJES/JERJES_%s.eps",cmult[ic]));
    }
    */
  }//! imult
  return 0;



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
