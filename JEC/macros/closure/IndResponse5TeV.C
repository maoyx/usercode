#include "../HiForest_V3/hiForest.h"
#include <TF1.h>
#include <TProfile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TFile.h>
#include <TStyle.h>
#include <TStopwatch.h>
#include <TRandom3.h>
#include <cstdlib>
#include <cmath>
#include <iostream>
using namespace std;


#ifdef _MAKECINT_
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ class HiForest;
#pragma link C++ class Evts;
#pragma link C++ class Jets;
#endif


int GetEtaBin(float /*eta*/);
int GetDetEtaBin(float /*eta*/);
int GetPhiBin(float /*phi*/);
int GetMultBin (int /*ntracks*/);
int GetHFplusEta4Bin (float /*HFplusEta4*/);
int GetPtBin(float /*pt*/);
bool selecJet(Jets */*jetc*/, int /*ij*/);
void FindLeadSubLeadJets(Jets */*jetc*/, int */*ljet*/);

//! pt binning

double ptbins[] ={10,20,30,50,70,90,120,160,200,240,280,340,400};
const int bins = sizeof(ptbins)/sizeof(double) - 1;

//! data pt binning
//double ptbins_data[] ={100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240, 300};
//const int dbins = sizeof(ptbins_data)/sizeof(double) - 1;


//! constants
const int iYear=2013;
const double pi=acos(-1.);
const double pi2=2*pi -1;

const int nmult=7; //! 6th is minbias, multiplicity bin
const int nhfbin=7; //! 6th is minbias, hiHFplusEta4 bin
const int ketar=4;
const int maxe=5;
const int maxph=5;
const float ketacut=3.0;
const double kptrecocut=10.;
const double kptgencut =0.;
const double kdRcut=0.3;
const double kVzcut=15.;

const int rbins=50;
double rbinl=0.0;
double rbinh=5.0;

TStopwatch timer;
int IndResponse5TeV(double kPt=30,const char *ksp="ppb",int ifile=0)
{

  timer.Start();

  //! Load Lib
  gSystem->Load("../HiForest_V3/hiForest_h.so");
  
  //! Load the hiforest input root file
  HiForest *c = 0;
  //! Latest
//  if(strcmp(ksp,"ppb")==0)c = new HiForest(Form("/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod13/Hijing_Pythia_pt%.f/HiForest_v72_v01_merged01/pt%.f_prod13_v72_merged_forest_%d.root",kPt,kPt,ifile),"PythiaHijing",cPPb); 
  if(strcmp(ksp,"ppb")==0)c = new HiForest(Form("/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt%.f/HiForest_v77_merged01/pt%.f_HP04_prod16_v77_merged_forest_%d.root",kPt,kPt,ifile),"PythiaHijing",cPPb); 
  else {
    c = new HiForest(Form("/mnt/hadoop/cms/store/user/kjung/pPbForest/Signal_Pythia_pt%.f/merged/pt%.f_HP04_prod16_v77_merged_forest_Sum.root",kPt,kPt),"Pythia",cPPb);
//    if(kPt==280)c = new HiForest(Form("/mnt/hadoop/cms/store/user/dgulhan/pPb/HP03/prod10/Signal_Pythia_pt%0.0f/HiForest_v63_v02_merged03/pt%0.0f_HP03_prod09_merged_forest_%d.root",kPt,kPt,ifile),"Pythia",cPPb);
//    else c = new HiForest(Form("/mnt/hadoop/cms/store/user/dgulhan/pPb/HP03/prod10/Signal_Pythia_pt%0.0f/HiForest_v63_merged03/pt%0.0f_HP03_prod09_merged_forest_%d.root",kPt,kPt,ifile),"Pythia",cPPb);
  }
  double xsection=0;
  double xup=0;
  double xsub=0;
  double maxpthat=9999;

  if(kPt==15){
    maxpthat=30;
    xup =5.335e-01;
    xsub=3.378e-02;
  }
  else if(kPt==30){
    maxpthat=50;
    xup =3.378e-02;
    xsub=3.778e-03;
  }
  else if(kPt==50){
    maxpthat=80;
    xup =3.778e-03;
    xsub=4.412e-04;
  }
  else if(kPt==80){
    maxpthat=120;
    xup =4.412e-04;
    xsub=6.147e-05;
  }
  else if(kPt==120){
    maxpthat=170;
    xup=6.147e-05;
    xsub=1.018e-05;
  }else if(kPt==170){
    maxpthat=220;
    xup=1.018e-05;
    xsub=2.477e-06;
  }else if(kPt==220){
    maxpthat=280;
    xup =2.477e-06;
    xsub=6.160e-07;
  }else if(kPt==280){
//    maxpthat=9999;
//    xup =6.160e-07;
//    xsub=0;
//  }
    maxpthat=370;
    xup =6.160e-07;
    xsub=1.088e-07;
    }else if(kPt==370){
    maxpthat=9999;
    xup=1.088e-07; 
    xsub=0;
    }

  xsection = xup-xsub;

  std::cout<<std::endl;
  std::cout<<std::endl;

  
  /*
  const int knj = 25;
  const char *cjets[knj] = {"icPu5",
			    "ak1PF","ak2PF","ak3PF","ak4PF","ak5PF","ak6PF",
			    "akPu1PF","akPu2PF","akPu3PF","akPu4PF","akPu5PF","akPu6PF",			    
			    "ak1Calo","ak2Calo","ak3Calo","ak4Calo","ak5Calo","ak6Calo",
			    "akPu1Calo","akPu2Calo","akPu3Calo","akPu4Calo","akPu5Calo","akPu6Calo"};

  const char *calgo[knj] = {"icPu5",
			    "ak1PF","ak2PF","ak3PF","ak4PF","ak5PF","ak6PF",
			    "akPu1PF","akPu2PF","akPu3PF","akPu4PF","akPu5PF","akPu6PF",			    
			    "ak1Calo","ak2Calo","ak3Calo","ak4Calo","ak5Calo","ak6Calo",
			    "akPu1Calo","akPu2Calo","akPu3Calo","akPu4Calo","akPu5Calo","akPu6Calo"};
  */
  const int knj = 16;
  const char *cjets[knj] = {"ak2PF","ak3PF","ak4PF","ak5PF",
			    "akPu2PF","akPu3PF","akPu4PF","akPu5PF",
			    "ak2Calo","ak3Calo","ak4Calo","ak5Calo",
			    "akPu2Calo","akPu3Calo","akPu4Calo","akPu5Calo"};
  const char *calgo[knj] = {"ak2PF","ak3PF","ak4PF","ak5PF",
			    "akPu2PF","akPu3PF","akPu4PF","akPu5PF",
			    "ak2Calo","ak3Calo","ak4Calo","ak5Calo",
			    "akPu2Calo","akPu3Calo","akPu4Calo","akPu5Calo"};

  c->hasIcPu5JetTree=0;
  c->hasAk1PFJetTree=0;
  c->hasAk6PFJetTree=0;
  c->hasAkPu1PFJetTree=0;
  c->hasAkPu6PFJetTree=0;
  c->hasAk1CaloJetTree=0;
  c->hasAk6CaloJetTree=0;
  c->hasAkPu1CaloJetTree=0;
  c->hasAkPu6CaloJetTree=0;


  /*
  c->hasMetTree=0;
  c->hasPFTree=0;
  c->hasIcPu5JetTree=0;
  c->hasAkPu2JetTree=0;
  c->hasAkPu3JetTree=0;
  c->hasAkPu4JetTree=0;
  c->hasAkPu2CaloJetTree=0;
  c->hasAkPu3CaloJetTree=0;
  c->hasAkPu4CaloJetTree=0;
  c->hasHltTree=0;
  c->hasTrackTree=0;
  c->hasPixTrackTree=0;
  c->hasSkimTree=0;
  c->hasTowerTree=0;
  c->hasHbheTree=0;
  c->hasEbTree=0;
  c->hasGenpTree=0;
  c->hasGenParticleTree=0;
  c->hasPhotonTree=0;
  */


  //! To get the jet object from hiforest
  Jets *iJet=0;

  std::cout<<"Loaded all tree variables and # of jet algorithms : "<<knj<<std::endl;
  std::cout<<"\t"<<std::endl;


  //! Open a output file for histos
  TFile *fout = new TFile(Form("Output/Response_%s_pT%0.0fGeV_%d.root",ksp,kPt,ifile),"RECREATE");


  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<Form("Running for %s ",ksp)<<std::endl;
  std::cout<<Form("pT  cut for %0.3f ",kptrecocut)<<std::endl;
  std::cout<<Form("eta cut for %0.3f ",ketacut)<<std::endl;
  std::cout<<"My hiForest TTree : " <<c->GetName()<<std::endl;
  std::cout<<"Output file  "<<fout->GetName()<<std::endl;
  std::cout<<"**************************************************** "<<std::endl;
  std::cout<<"\t"<<std::endl;
  std::cout<<"\t"<<std::endl;


  //! 
  //! Define histograms here
  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  TProfile::SetDefaultSumw2();

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  TH1F *hEvt    = new TH1F("hEvt","# of events ",4,0,4);
  TH1F *hVz     = new TH1F("hVz","# of events ",80,-20,20);
  TH1F *hBin    = new TH1F("hBin","Centrality bin",100,-0.5,100-0.5);
  TH1F *hTotEve = new TH1F("hTotEve","# of events in the skimmed files",4,0,4);

  
  TH1F *hHF       = new TH1F("hHF","HF distribution",500,0,500);
  TH1F *hNTracks  = new TH1F("hNTracks","hiNTracks",500,0,500);
  TH2F *hNTracksHF = new TH2F("hNTracksHF","Ntrack vs HF ",500,0,500,500,0,500);
  TH2F *hBinHF      = new TH2F("hBinHF","HF vs hiBin",100,-0.5,100-0.5,500,0,500);
  TH2F *hBinNTracks = new TH2F("hBinNTracks","HF vs NTracks",100,-0.5,100-0.5,500,0,500);
    TH2F *hNTracksHFplusEta4 = new TH2F("hNTracksHFplusEta4","NTracks vs HFplusEta4",500,0,500, 500, 0., 500);
  
  TH1F *hgenpt_genm [knj][nmult], *hrecopt_genm[knj][nmult], *hrawpt_genm[knj][nmult];
  TH1F *hjeteta  [knj][nmult], *hjetphi[knj][nmult];
  TH2F *hjetpteta[knj][nmult], *hjetptphi[knj][nmult], *hjetetaphi[knj][nmult];

  //! Ratios of the pt distributions
  TProfile *hrecogen[knj][nmult], *hrecoraw[knj][nmult], *hrawgen[knj][nmult];

  //! Resposnse
  TH2F *hrescrpt_genm[knj][nmult], *hresrrpt_genm[knj][nmult], *hresrcrpt_genm[knj][nmult];
  TH2F *hrescreta_genm[knj][nmult], *hresrreta_genm[knj][nmult], *hresrcreta_genm[knj][nmult];

  TH2F *hratiorawrefpt_eta[knj][nmult][2], *hratiocorrrefpt_eta[knj][nmult][2];

  TH2F *hratiocorrrefpt_genm[knj][nmult];
  TH2F *hratiocorrrefpt_lead[knj][nmult], *hratiocorrrefpt_slead[knj][nmult];

    TH2F *hratiocorrrefpt_genhfb[knj][nhfbin];
    TH2F *hratiocorrrefpt_etaptbin[knj][bins];

  //! For comparison with data
  TH2F *hJetEnergyScale[knj][nmult];

  TH2F *hpteta[knj][nmult][maxe] ;
  TH2F *hptphi[knj][nmult][maxph] ;

  TH2F *hgenjrecoj[knj][nmult];
  TH1F *hNjets_genm[knj][nmult];
  TH1F *hNevt [knj][nmult];

  //! Background for jets 
  //! (photonSum+neutralSum+chargedSum-rawpt)
  TH2F *hjetptpu_genm[knj][nmult];             //! centrality                                                                       
  TH2F *hjetptpu_etab_genm[knj][nmult][ketar]; //! eta dependence             
  TH1F *hjetbkgd_genm[knj][nmult];
  TH2F *hjetptbkgd_genm[knj][nmult];
  TH2F *hjetptbkgd_etab_genm[knj][nmult][ketar];
  TH2F *hPFFraction_genm[knj][nmult][3];

  //! Efficency histos
  TH1F *hPtAll [knj][nmult], *hPtSel[knj][nmult];
  TH1F *hEtaAll[knj][nmult], *hEtaSel[knj][nmult];
  TH1F *hPhiAll[knj][nmult], *hPhiSel[knj][nmult];

  //! Response vs deltar
  TH1F *hRspVsDeltaR[knj][nmult][25];

  //! DeltaR efficiency                                                                                                                                                              
  TH1F *hDeltaR[knj][nmult];
  TH1F *hDeltaRAll[knj][nmult];
  TH1F *hDeltaRSel[knj][nmult];

    for(int nj=0;nj<knj;nj++){
        
        for(int icen=0;icen<nmult;icen++){
            hNevt [nj][icen] = new TH1F(Form("hNevt%d_%d",nj,icen),Form("# of events cent %d %s",icen,cjets[nj]),40,-40,40);
            
            hNjets_genm [nj][icen] = new TH1F(Form("hNjets_genm%d_%d",nj,icen),Form("# of jets cent %d jets %s",icen,cjets[nj]),200,30,630);
            hgenpt_genm [nj][icen] = new TH1F(Form("hgenpt_genm%d_%d",nj,icen),Form("Gen matched gen p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);
            hrecopt_genm[nj][icen] = new TH1F(Form("hrecopt_genm%d_%d",nj,icen),Form("Gen matched reco p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);
            hrawpt_genm [nj][icen] = new TH1F(Form("hrawpt_genm%d_%d",nj,icen),Form("Gen matched raw p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);
            
            //! Ratios
            hrecogen[nj][icen] = new TProfile(Form("hrecogen%d_%d",nj,icen),Form("reco/gen p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);
            hrecoraw[nj][icen] = new TProfile(Form("hrecoraw%d_%d",nj,icen),Form("reco/raw p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);
            hrawgen[nj][icen]  = new TProfile(Form("hrawgen%d_%d",nj,icen),Form("raw/gen p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000);
            
            //! Gen matched
            hrescrpt_genm[nj][icen]= new TH2F(Form("hrescrpt_genm%d_%d",nj,icen),Form("Gen jet:(Reco/Gen) jet p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000,150,rbinl,rbinh);
            hresrrpt_genm[nj][icen]= new TH2F(Form("hresrrpt_genm%d_%d",nj,icen),Form("Gen jet:(Raw/Gen) jet p_{T}  distribution jet centb %d %s",icen,cjets[nj]),500,0,1000,150,rbinl,rbinh);
            hresrcrpt_genm[nj][icen]= new TH2F(Form("hresrcrpt_genm%d_%d",nj,icen),Form("Reco jet:(Reco/Raw) jet p_{T} distribution jet centb %d %s",icen,cjets[nj]),500,0,1000,150,rbinl,rbinh);
            
            hrescreta_genm[nj][icen]= new TH2F(Form("hrescreta_genm%d_%d",nj,icen),Form("Gen jet:(Reco/Gen) jet p_{T} distribution jet centb %d %s",icen,cjets[nj]),100,-5.0,5.0,150,rbinl,rbinh);
            hresrreta_genm[nj][icen]= new TH2F(Form("hresrreta_genm%d_%d",nj,icen),Form("Gen jet:(Raw/Gen) jet p_{T}  distribution jet centb %d %s",icen,cjets[nj]),100,-5.0,5.0,150,rbinl,rbinh);
            hresrcreta_genm[nj][icen]= new TH2F(Form("hresrcreta_genm%d_%d",nj,icen),Form("Reco jet:(Reco/Raw) jet p_{T} distribution jet centb %d %s",icen,cjets[nj]),100,-5.0,5.0,150,rbinl,rbinh);
            
            hjeteta[nj][icen] = new TH1F(Form("hjeteta%d_%d",nj,icen),Form("jet eta distribution jet centb %d %s",icen,calgo[nj]),100,-5.0,5.0);
            hjetphi[nj][icen] = new TH1F(Form("hjetphi%d_%d",nj,icen),Form("jet phi distribution jet centb %d %s",icen,calgo[nj]),72,-pi,pi);
            
            hjetetaphi[nj][icen] = new TH2F(Form("hjetetaphi%d_%d",nj,icen),Form("jet eta-phi distribution jet centb %d %s",icen,calgo[nj]),100,-5.0,5.0,72,-pi,pi);
            hjetpteta[nj][icen] = new TH2F(Form("hjetpteta%d_%d",nj,icen),Form("jet pt-eta distribution jet centb %d %s",icen,calgo[nj]),500,0,1000,100,-5.0,5.0);
            hjetptphi[nj][icen] = new TH2F(Form("hjetptphi%d_%d",nj,icen),Form("jet pt-phi distribution jet centb %d %s",icen,calgo[nj]),500,0,1000,72,-pi,pi);
            
            
            hratiocorrrefpt_lead[nj][icen]= new TH2F(Form("hratiocorrrefpt_lead%d_%d",nj,icen),Form("Leading jet Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",icen,cjets[nj]),
                                                     500,0,1000,rbins,rbinl,rbinh);
            hratiocorrrefpt_slead[nj][icen]= new TH2F(Form("hratiocorrrefpt_slead%d_%d",nj,icen),Form("sub-Leading jet Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",icen,cjets[nj]),
                                                      500,0,1000,rbins,rbinl,rbinh);
            hratiocorrrefpt_genm[nj][icen]= new TH2F(Form("hratiocorrrefpt_genm%d_%d",nj,icen),Form("Gen matched jet Reco jet / Gen jet p_{T} (corr.) distribution jet centb %d %s",icen,cjets[nj]),
                                                     500,0,1000,rbins,rbinl,rbinh);
            hJetEnergyScale[nj][icen] = new TH2F(Form("hJetEnergyScale%d_%d",nj,icen),Form("hJetEnergyScale%d_%d",nj,icen),500,0,1000,50,-1.00,1.00);
            
            for(int ie=0;ie<2;ie++){
                hratiorawrefpt_eta[nj][icen][ie]= new TH2F(Form("hratiorawrefpt_eta%d_%d_%d",nj,icen,ie),
                                                           Form("Raw jet / Gen jet p_{T} (raw) distribution jet centb %d %s etabin%d",icen,cjets[nj],ie),
                                                           500,0,1000,rbins,rbinl,rbinh);
                hratiocorrrefpt_eta[nj][icen][ie]= new TH2F(Form("hratiocorrrefpt_eta%d_%d_%d",nj,icen,ie),
                                                            Form("Reco jet / Gen jet p_{T} (raw) distribution jet centb %d %s etabin%d",icen,cjets[nj],ie),
                                                            500,0,1000,rbins,rbinl,rbinh);
            }
            
            //hjetptpu_genm[nj][icen] = new TH2F(Form("hjetptpu_genm%d_%d",nj,icen),Form("jet(pt:pu) distribution jet centb %d %s",icen,cjets[nj]),dbins,ptbins_data,100,0,300);
            //for(int ie=0;ie<ketar;ie++){
            // hjetptpu_etab_genm[nj][icen][ie] = new TH2F(Form("hjetptpu_etab_genm%d_%d_%d",nj,icen,ie),Form("jet(pt:pu) distribution jet %s etabin %d cen %d",cjets[nj],icen,ie),
            ///					    dbins,ptbins_data,100,0,300);
            //}
            
            hjetptpu_genm[nj][icen] = new TH2F(Form("hjetptpu_genm%d_%d",nj,icen),Form("jet(pt:pu) distribution jet centb %d %s",icen,cjets[nj]),500,0,1000,100,0,300);
            for(int ie=0;ie<ketar;ie++){
                hjetptpu_etab_genm[nj][icen][ie] = new TH2F(Form("hjetptpu_etab_genm%d_%d_%d",nj,icen,ie),Form("jet(pt:pu) distribution jet %s etabin %d cen %d",cjets[nj],icen,ie),
                                                            500,0,1000,100,0,300);
                hjetptbkgd_etab_genm[nj][icen][ie] = new TH2F(Form("hjetptbkgd_etab_genm%d_%d_%d",nj,icen,ie),Form("jet(pt:bkgd) distribution jet etabin %d centb %d %s",ie,icen,calgo[nj]),500,0,1000,100,0,300);
            }
            hjetbkgd_genm  [nj][icen] = new TH1F(Form("hjetbkgd_genm%d_%d",nj,icen),Form("jet(pu) p_{T} distribution jet centb %d %s",icen,calgo[nj]),100,0,300);
            //hjetptbkgd_genm[nj][icen] = new TH2F(Form("hjetptbkgd_genm%d_%d",nj,icen),Form("jet(pt:bkgd) distribution jet centb %d %s",icen,calgo[nj]),dbins,ptbins_data,100,0,300);
            hjetptbkgd_genm[nj][icen] = new TH2F(Form("hjetptbkgd_genm%d_%d",nj,icen),Form("jet(pt:bkgd) distribution jet centb %d %s",icen,calgo[nj]),500,0,1000,100,0,300);
            //! PF fractions
            for(int ipf=0;ipf<3;ipf++){
                hPFFraction_genm[nj][icen][ipf] = new TH2F(Form("hPFFraction_genm%d_%d_%d",nj,icen,ipf),Form("PF fraction distribution jet centb %d %s %d",icen,calgo[nj],ipf),500,0,1000,15,0,1.5);
            }
            
            
            for(int m=0;m<maxe;m++){
                hpteta[nj][icen][m] = new TH2F(Form("hpteta%d_%d_%d",nj,icen,m),Form("resolution  pt(eta) distribution cent %d jet %s etabin%d",icen,cjets[nj],m),
                                               500,0,1000,rbins,rbinl,rbinh);
            }
            
            for(int m=0;m<maxph;m++){
                hptphi[nj][icen][m] = new TH2F(Form("hptphi%d_%d_%d",nj,icen,m),Form("resolution pt(phi) distribution cent %d jet %s phibin%d",icen,cjets[nj],m),
                                               500,0,1000,rbins,rbinl,rbinh);
            }
            
            hgenjrecoj[nj][icen] =new TH2F(Form("hgenjrecoj%d_%d",nj,icen),Form("gen jet2 : reco jet2 %s cent %d",cjets[nj],icen),500,0.,1000.,500,0.,1000.);
            
            //! efficency histograms
            hPtAll [nj][icen] = new TH1F(Form("hPtAll%d_%d",nj,icen),Form("Denominator pT for algorithm %s cent %d",cjets[nj],icen),40,10,110);
            hEtaAll[nj][icen] = new TH1F(Form("hEtaAll%d_%d",nj,icen),Form("Denominator eta  for algorithm %s cent %d",cjets[nj],icen),100,-5.0,5.0);
            hPhiAll[nj][icen] = new TH1F(Form("hPhiAll%d_%d",nj,icen),Form("Denominator  phi  for algorithm %s cent %d",cjets[nj],icen),20,-pi,pi);
            
            hPtSel [nj][icen] = new TH1F(Form("hPtSel%d_%d",nj,icen),Form("Numerator pT for algorithm %s cent %d",cjets[nj],icen),40,10,110);
            hEtaSel[nj][icen] = new TH1F(Form("hEtaSel%d_%d",nj,icen),Form("Numerator eta  for algorithm %s cent %d",cjets[nj],icen),100,-5.0,5.0);
            hPhiSel[nj][icen] = new TH1F(Form("hPhiSel%d_%d",nj,icen),Form("Numerator  phi  for algorithm %s cent %d",cjets[nj],icen),20,-pi,pi);
            
            hDeltaR[nj][icen]    = new TH1F(Form("hDeltaR%d_%d",nj,icen),Form("#DeltaR for algorithm %s cent %d",cjets[nj],icen),100,0,1);
            hDeltaRAll[nj][icen] = new TH1F(Form("hDeltaRAll%d_%d",nj,icen),Form("#DeltaR (all) for algorithm %s cent %d",cjets[nj],icen),100,0,1);
            hDeltaRSel[nj][icen] = new TH1F(Form("hDeltaRSel%d_%d",nj,icen),Form("#DeltaR (sel) for algorithm %s cent %d",cjets[nj],icen),100,0,1);
            
            for(int ir=0;ir<25;ir++){
                //! Response vs DeltaR
                hRspVsDeltaR[nj][icen][ir] = new TH1F(Form("hRspVsDeltaR%d_%d_%d",nj,icen,ir),Form(" <recopt/refpt> vs. #DeltaR (%d) algorithm %s cent %d",ir,cjets[nj],icen),rbins,rbinl,rbinh);
            }
        }//! icen
        //for eta clsoure in different jet pt bins
        for(int ipt=0;ipt<bins;ipt++){
            hratiocorrrefpt_etaptbin[nj][ipt]= new TH2F(Form("hratiocorrrefpt_etaptbin%d_%d",nj,ipt),Form("Gen matched jet Reco jet / Gen jet p_{T} (corr.) distribution jet vs eta in pt bin %d %s",ipt,cjets[nj]),
                                                  100,-5.0,5.0,rbins,rbinl,rbinh);
            
        }
        //for different HFplusEta4 bins
        for(int ihf=0;ihf<nhfbin;ihf++){
            hratiocorrrefpt_genhfb[nj][ihf]= new TH2F(Form("hratiocorrrefpt_genhfbin%d_%d",nj,ihf),Form("Gen matched jet Reco jet / Gen jet p_{T} (corr.) distribution jet hiHF bin %d %s",ihf,cjets[nj]),
                                                      500,0,1000,rbins,rbinl,rbinh);
        }
        
    }//! nj
  std::cout<<"Initialized the histograms " <<std::endl;
  ///////////////////////////////////////////////////////////////////////////////////////// 


  //! Centrality reweighting function

  //! vertex z reweighting




  Long64_t nentries = c->GetEntries();
  std::cout<<Form("# of entries in TTree for %s : ",ksp)<<nentries<<std::endl;
  std::cout<<std::endl;
  hEvt->Fill(2,nentries);

  //! weight  for the merging of the samples for different pT hat bins
  Float_t wxs = xsection/(nentries/100000.);
  
  Int_t iEvent=0; 
  for (Long64_t ievt=0; ievt<nentries;ievt++) {//! event loop
    //for (Long64_t ievt=0; ievt<5000;ievt++) {//! event loop
    //! load the hiForest event
    c->GetEntry(ievt);

    int hiBin       = c->evt.hiBin;
    float vz        = c->evt.vz;
    float hiHF      = c->evt.hiHF;
    int ntracks     = c->evt.hiNtracks;
    float HFplusEta = c->evt.hiHFplusEta4;
    float HFminusEta = c->evt.hiHFminusEta4 ; 
    float HFplusEta4 = HFplusEta+HFminusEta;
    //! testing
    //if(hiBin>4 && strcmp(ksp,"pbpb")==0)continue;

    //! apply vertex cut
    if(fabs(vz)>kVzcut)continue;


    //! Centrality bin
    if(hiBin<0 || hiBin>100)continue;
    int multb=GetMultBin(ntracks);
    int hiHFb=GetHFplusEta4Bin(HFplusEta4);  
      
    double wcen=1;
    double wvz=1;
    //wxs=1;
    if(strcmp(ksp,"pp")==0)multb=0;
  
    if(ievt%10000==0)std::cout<<" ********** Event # " <<ievt<<std::endl;
    //std::cout<<" ********** Event # " <<ievt<<"\t vz : "<<vz<<"\t hiBin : "<<hiBin<<"\t wxs : "<<wxs<<std::endl;


    //! Centrality from 0-100% 
    if(multb==-1 || multb==nmult)continue;
    if(hiHFb==-1 || hiHFb==nhfbin)continue;

    int istat=0;
    for(int nj=0;nj<knj;nj++){ //! loop over different jet algorithms
      

      if(nj==0)iJet = &(c->ak2PF);
      else if(nj==1)iJet = &(c->ak3PF);
      else if(nj==2)iJet = &(c->ak4PF);
      else if(nj==3)iJet = &(c->ak5PF);

      else if(nj==4)iJet = &(c->akPu2PF);
      else if(nj==5)iJet = &(c->akPu3PF);
      else if(nj==6)iJet = &(c->akPu4PF);
      else if(nj==7)iJet = &(c->akPu5PF);

      else if(nj==8)iJet = &(c->ak2Calo);
      else if(nj==9)iJet = &(c->ak3Calo);
      else if(nj==10)iJet = &(c->ak4Calo);
      else if(nj==11)iJet = &(c->ak5Calo);

      else if(nj==12)iJet = &(c->akPu2Calo);
      else if(nj==13)iJet = &(c->akPu3Calo);
      else if(nj==14)iJet = &(c->akPu4Calo);
      else if(nj==15)iJet = &(c->akPu5Calo);


      //! xsec-weight
      double pthat = iJet->pthat;
      if(pthat > maxpthat)continue;
      istat=1;
      
      //std::cout<<"\t Jet Algorithm : "<<cjets[nj]<<"\t # of Jets  : "<<iJet->nref<<"\t pthat : "<<pthat<<std::endl;
      if(nj==0)hTotEve->Fill(1); //! akPu3PF      


      int *ljet = new int[3];
      FindLeadSubLeadJets(iJet,ljet);
      if(ljet[0]>=0){
	hratiocorrrefpt_lead[nj][multb]->Fill(iJet->refpt[ljet[0]],iJet->jtpt[ljet[0]]/iJet->refpt[ljet[0]],wxs*wcen*wvz);
	hratiocorrrefpt_lead[nj][nmult-1]->Fill(iJet->refpt[ljet[0]],iJet->jtpt[ljet[0]]/iJet->refpt[ljet[0]],wxs*wcen*wvz);
      }
      if(ljet[1]>=0){
	hratiocorrrefpt_slead[nj][multb]->Fill(iJet->refpt[ljet[1]],iJet->jtpt[ljet[1]]/iJet->refpt[ljet[1]],wxs*wcen*wvz);
	hratiocorrrefpt_slead[nj][nmult-1]->Fill(iJet->refpt[ljet[1]],iJet->jtpt[ljet[1]]/iJet->refpt[ljet[1]],wxs*wcen*wvz);
      }

      //! Jet energy scale comparison with data
      if(ljet[0]>=0 && ljet[1]>=0 && iJet->jtpt[ljet[0]]>50. && iJet->jtpt[ljet[1]]>50.){//! atleas a dijet
        int mstat=1;
        double ptdij = (iJet->jtpt[ljet[0]] + iJet->jtpt[ljet[1]])/2.;
        if(ljet[2]>=0){
          //if(iJet->jtpt[ljet[2]]/ptdij > 0.2)mstat=0;
	  mstat=0;
        }
        if(mstat){
          double B=-9999;
          double rn1 = gRandom->Rndm();
          double rn2 = gRandom->Rndm();
          if(rn1 > rn2){
            B = (iJet->jtpt[ljet[0]] - iJet->jtpt[ljet[1]])/(iJet->jtpt[ljet[0]] + iJet->jtpt[ljet[1]]);
          }else{
            B = (iJet->jtpt[ljet[1]] - iJet->jtpt[ljet[0]])/(iJet->jtpt[ljet[1]] + iJet->jtpt[ljet[0]]);
          }
          hJetEnergyScale[nj][multb]->Fill(ptdij,B);
          hJetEnergyScale[nj][nmult-1]->Fill(ptdij,B);
	}
      }

      //! Gen matched jets
      for(int igen=0; igen<iJet->nref; igen++){
	if( iJet->subid[igen] != 0) continue;
	int gj = igen;
	
	float rawpt   = iJet->rawpt[gj];
	float refpt   = iJet->refpt[gj];
	float refeta  = iJet->refeta[gj];
	float refphi  = iJet->refphi[gj];
	float recopt  = iJet->jtpt[gj];
	float recoeta = iJet->jteta[gj];
	float delr    = iJet->refdrjt[gj];
          
       if(rawpt <15.) continue ;
     //  if(recopt<30.) continue ;
	if(fabs(refeta)<ketacut && refpt>10){
          //! Denominator for matching efficiency
	  hPtAll [nj][multb]->Fill(refpt,wxs*wcen*wvz);
          hEtaAll[nj][multb]->Fill(refeta,wxs*wcen*wvz);
          hPhiAll[nj][multb]->Fill(refphi,wxs*wcen*wvz);
	  
	  hPtAll [nj][nmult-1]->Fill(refpt,wxs*wcen*wvz);
          hEtaAll[nj][nmult-1]->Fill(refeta,wxs*wcen*wvz);
          hPhiAll[nj][nmult-1]->Fill(refphi,wxs*wcen*wvz);
	  
	  //! DeltaR efficiency
	  hDeltaR[nj][multb]->Fill(delr,wxs*wcen*wvz);
	  hDeltaR[nj][nmult-1]->Fill(delr,wxs*wcen*wvz);
	  for (int idrmax=0;idrmax<100;idrmax++) {
	    float drmax = idrmax*0.01+0.005;
	    hDeltaRAll[nj][multb]->Fill(drmax,wxs*wcen*wvz);
	    hDeltaRAll[nj][nmult-1]->Fill(drmax,wxs*wcen*wvz);
	    if (delr<drmax){
	      hDeltaRSel[nj][multb]->Fill(drmax,wxs*wcen*wvz);
	      hDeltaRSel[nj][nmult-1]->Fill(drmax,wxs*wcen*wvz);
	    }
	  }
        }

	
	if(recopt<kptrecocut || refpt<kptgencut || refpt==0 || fabs(recoeta)>ketacut || fabs(delr)>kdRcut)continue;

        if(fabs(refeta)<ketacut && refpt>10){
          //! Numerator for matching efficiency
          hPtSel [nj][multb]->Fill(refpt,wxs*wcen*wvz);
          hEtaSel[nj][multb]->Fill(refeta,wxs*wcen*wvz);
          hPhiSel[nj][multb]->Fill(refphi,wxs*wcen*wvz);

          hPtSel [nj][nmult-1]->Fill(refpt,wxs*wcen*wvz);
          hEtaSel[nj][nmult-1]->Fill(refeta,wxs*wcen*wvz);
          hPhiSel[nj][nmult-1]->Fill(refphi,wxs*wcen*wvz);
        }
	

	//! Response
	for (int idr=0;idr<25;idr++) {
	  double drcut = 0.0+idr*(0.25-0.00)/(25-1);
	  if (delr>drcut) continue;
	  hRspVsDeltaR[nj][multb][idr]->Fill(recopt/refpt,wxs*wcen*wvz);
	  hRspVsDeltaR[nj][nmult-1][idr]->Fill(recopt/refpt,wxs*wcen*wvz);
	}
	

    hratiocorrrefpt_genm[nj][multb]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	hratiocorrrefpt_genhfb[nj][hiHFb]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	hrecogen[nj][multb]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	hrecoraw[nj][multb]->Fill(recopt,recopt/rawpt,wxs*wcen*wvz);
        hrawgen [nj][multb]->Fill(refpt,rawpt/refpt,wxs*wcen*wvz);
	hgenjrecoj [nj][multb]->Fill(refpt,recopt,wxs*wcen*wvz);


	hratiocorrrefpt_genm[nj][nmult-1]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
        hratiocorrrefpt_genhfb[nj][nhfbin-1]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
        hratiocorrrefpt_genhfb[nj][nhfbin-2]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	hrecogen[nj][nmult-1]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	hrecoraw[nj][nmult-1]->Fill(recopt,recopt/rawpt,wxs*wcen*wvz);
        hrawgen [nj][nmult-1]->Fill(refpt,rawpt/refpt,wxs*wcen*wvz);
	hgenjrecoj [nj][nmult-1]->Fill(refpt,recopt,wxs*wcen*wvz);

	int ieta=-1;
	if(fabs(recoeta)<1.3)ieta=0; //! barrel region
	else ieta=1; //! HCAL region

   int ptb=GetPtBin(refpt);
   hratiocorrrefpt_etaptbin[nj][ptb]->Fill(refeta,recopt/refpt,wxs*wcen*wvz);
          
          
	//! Response in eta
	hratiocorrrefpt_eta[nj][multb][ieta]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	hratiorawrefpt_eta [nj][multb][ieta]->Fill(refpt,rawpt/refpt,wxs*wcen*wvz);
	hratiocorrrefpt_eta[nj][nmult-1][ieta]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	hratiorawrefpt_eta [nj][nmult-1][ieta]->Fill(refpt,rawpt/refpt,wxs*wcen*wvz);

	hjeteta   [nj][multb]->Fill(iJet->jteta[gj]);
	hjetphi   [nj][multb]->Fill(iJet->jtphi[gj]);
	hjetpteta [nj][multb]->Fill(iJet->jtpt[gj],iJet->jteta[gj]);
	hjetptphi [nj][multb]->Fill(iJet->jtpt[gj],iJet->jtphi[gj]);
	hjetetaphi[nj][multb]->Fill(iJet->jteta[gj],iJet->jtphi[gj]);

	hjeteta   [nj][nmult-1]->Fill(iJet->jteta[gj]);
	hjetphi   [nj][nmult-1]->Fill(iJet->jtphi[gj]);
	hjetpteta [nj][nmult-1]->Fill(iJet->jtpt[gj],iJet->jteta[gj]);
	hjetptphi [nj][nmult-1]->Fill(iJet->jtpt[gj],iJet->jtphi[gj]);
	hjetetaphi[nj][nmult-1]->Fill(iJet->jteta[gj],iJet->jtphi[gj]);



	//! Fill the background pT info for the signal matched jets only
	//! pileup study
	hjetptpu_genm[nj][multb]->Fill(recopt,iJet->jtpu[gj],wxs*wcen*wvz);
	hjetptpu_etab_genm[nj][multb][ieta]->Fill(recopt,iJet->jtpu[gj],wxs*wcen*wvz);

	hjetptpu_genm[nj][nmult-1]->Fill(recopt,iJet->jtpu[gj],wxs*wcen*wvz);
	hjetptpu_etab_genm[nj][nmult-1][ieta]->Fill(recopt,iJet->jtpu[gj],wxs*wcen*wvz);
	
	//! Jet bkgd estimation
	//! (photonSum+neutralSum+chargedSum-rawpt)
	double jbkgd  = (iJet->photonSum[gj]+iJet->neutralSum[gj]+iJet->chargedSum[gj]) - rawpt;
	hjetbkgd_genm  [nj][multb]->Fill(jbkgd,wxs*wcen*wvz);
	hjetptbkgd_genm[nj][multb]->Fill(recopt,jbkgd,wxs*wcen*wvz);
	hjetptbkgd_etab_genm[nj][multb][ieta]->Fill(recopt,jbkgd,wxs*wcen*wvz);

	hjetbkgd_genm  [nj][nmult-1]->Fill(jbkgd,wxs*wcen*wvz);
	hjetptbkgd_genm[nj][nmult-1]->Fill(recopt,jbkgd,wxs*wcen*wvz);
	hjetptbkgd_etab_genm[nj][nmult-1][ieta]->Fill(recopt,jbkgd,wxs*wcen*wvz);

	hPFFraction_genm[nj][multb][0]->Fill(recopt,iJet->photonSum[gj]/rawpt);
	hPFFraction_genm[nj][multb][1]->Fill(recopt,iJet->neutralSum[gj]/rawpt);
	hPFFraction_genm[nj][multb][2]->Fill(recopt,iJet->chargedSum[gj]/rawpt);

	hPFFraction_genm[nj][nmult-1][0]->Fill(recopt,iJet->photonSum[gj]/rawpt);
	hPFFraction_genm[nj][nmult-1][1]->Fill(recopt,iJet->neutralSum[gj]/rawpt);
	hPFFraction_genm[nj][nmult-1][2]->Fill(recopt,iJet->chargedSum[gj]/rawpt);
	
	hNjets_genm [nj][multb]->Fill(refpt);
	hgenpt_genm [nj][multb]->Fill(refpt,wxs*wcen*wvz);
	hrecopt_genm[nj][multb]->Fill(recopt,wxs*wcen*wvz);	  
	hrawpt_genm [nj][multb]->Fill(rawpt,wxs*wcen*wvz);

	hNjets_genm [nj][nmult-1]->Fill(refpt);
	hgenpt_genm [nj][nmult-1]->Fill(refpt,wxs*wcen*wvz);
	hrecopt_genm[nj][nmult-1]->Fill(recopt,wxs*wcen*wvz);	  
	hrawpt_genm [nj][nmult-1]->Fill(rawpt,wxs*wcen*wvz);
	
	//! Very fine bin in ref pt
	hrescrpt_genm[nj][multb]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	hresrrpt_genm[nj][multb]->Fill(refpt,rawpt/refpt,wxs*wcen*wvz);
	hresrcrpt_genm[nj][multb]->Fill(recopt,recopt/rawpt,wxs*wcen*wvz);

	hrescrpt_genm [nj][nmult-1]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	hresrrpt_genm [nj][nmult-1]->Fill(refpt,rawpt/refpt,wxs*wcen*wvz);
	hresrcrpt_genm[nj][nmult-1]->Fill(recopt,recopt/rawpt,wxs*wcen*wvz);

	//! Very fine bin in ref eta
	hrescreta_genm[nj][multb]->Fill(refeta,recopt/refpt,wxs*wcen*wvz);
	hresrreta_genm[nj][multb]->Fill(refeta,rawpt/refpt,wxs*wcen*wvz);
	hresrcreta_genm[nj][multb]->Fill(recoeta,recopt/rawpt,wxs*wcen*wvz);

	hrescreta_genm[nj][nmult-1]->Fill(refeta,recopt/refpt,wxs*wcen*wvz);
	hresrreta_genm[nj][nmult-1]->Fill(refeta,rawpt/refpt,wxs*wcen*wvz);
	hresrcreta_genm[nj][nmult-1]->Fill(recoeta,recopt/rawpt,wxs*wcen*wvz);


        //! Response in different eta and phi bins
        int etabin = GetEtaBin(fabs(refeta));
	int phibin = GetPhiBin(refphi);

        //! Response in eta and phi bins
        if(etabin >= 0 && etabin<maxe){
	  hpteta[nj][multb][etabin]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	  hpteta[nj][nmult-1][etabin]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	}
        if(phibin >= 0 && phibin<maxph){
	  hptphi[nj][multb][phibin]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	  hptphi[nj][nmult-1][phibin]->Fill(refpt,recopt/refpt,wxs*wcen*wvz);
	}
      }//! igen loop
      hNevt[nj][multb]->Fill(vz,wxs*wcen*wvz);
      hNevt[nj][nmult-1]->Fill(vz,wxs*wcen*wvz);

      delete [] ljet;
    }//! nj jet loop

    if(istat){
      hBin->Fill(hiBin,wxs*wcen*wvz);
      hVz->Fill(vz,wxs*wcen*wvz);
      hHF->Fill(hiHF,wxs*wcen*wvz);

      hNTracks->Fill(ntracks,wxs*wcen*wvz);
      hNTracksHF->Fill(ntracks,hiHF,wxs*wcen*wvz);
      hBinHF->Fill(hiBin,hiHF,wxs*wcen*wvz);
      hBinNTracks->Fill(hiBin,ntracks,wxs*wcen*wvz);
      hNTracksHFplusEta4->Fill(ntracks,HFplusEta4,wxs*wcen*wvz);
      iEvent++;
    }
    //std::cout<<"Completed event #  "<<ievt<<std::endl; 
  }//! event loop ends
  
  std::cout<<std::endl;
  std::cout<<std::endl;
  std::cout<<std::endl;

  std::cout<<"Events which passed the pT hat cut : "<<hTotEve->Integral()<<" out of  : "<<hEvt->Integral()
	   <<" efficiency of all the cuts : " <<hTotEve->Integral()/hEvt->Integral()<<std::endl;
  std::cout<<std::endl;


  for(int nj=0;nj<knj;nj++){
    std::cout<<"# of Events for : "<<cjets[nj]<<"\t"<<hNevt[nj][nmult-1]->Integral()<<"\t # of Jets : "<<hNjets_genm[nj][nmult-1]->Integral()<<std::endl;
  }

  //! Write to output file
  fout->cd();
  fout->Write();
  fout->Close();


  //! Check
  timer.Stop();
  double rtime  = timer.RealTime();
  double ctime  = timer.CpuTime();

  std::cout<<std::endl;
  std::cout<<Form("RealTime=%f seconds, CpuTime=%f seconds",rtime,ctime)<<std::endl;
  std::cout<<std::endl;
  std::cout<<"Good bye : " <<"\t"<<std::endl;
  return 1;
}
/*
int GetCentbin(int bin)
{
  int ibin=-1;
  //! centrality is defined as 1 bins of cross section
  //! in 0-100 bins

  if(bin<10)ibin=0; //! 0-10%
  else if(bin>=10 && bin<30)ibin=1;   //! 10-30%
  else if(bin>=30 && bin<50)ibin=2;   //! 30-50%
  else if(bin>=50 && bin<70)ibin=3;   //! 50-70%
  else if(bin>=70 && bin<100)ibin=4;  //! 70-100%

  return ibin;
}
*/
int GetEtaBin(float eta)
{
  int ibin=-1;
  float min =0.0;
  ibin = (eta - min)*maxe/2.;
  return ibin;
}
int GetPhiBin(float phi)
{
  int ibin=-1;
  float min = -pi;
  ibin = (phi - min)*maxph/2*pi;
  return ibin;
}
bool selecJet(Jets *jetc,int it)
{
  bool goodJet = jetc->refpt[it]!=-999 && jetc->refphi[it]!=-999 && jetc->refeta[it]!=-999 && jetc->jtpt[it]!=-999 && jetc->jtphi[it]!=-999 && jetc->jteta[it]!=-999;
  return goodJet;
}
int GetDetEtaBin(float eta)
{
  int ibin=-1;
  if(eta>=0 && eta<1.3)ibin=0;       //! barrel
  else if(eta>=1.3 && eta<2.0)ibin=1;//! endcap+tracks
  else if(eta>=2.0 && eta<3.0)ibin=2;//! endcap-notracks
  else if(eta>=3.0 && eta<5.0)ibin=3;//! forward

  return ibin;
}
int GetPtBin(float pt)
{
  for(int ix=0;ix<bins;ix++){
    if(pt>=ptbins[ix] && pt<ptbins[ix+1]){
      return ix;
    }
  }
  return -1;
}
void FindLeadSubLeadJets(Jets *jetc, int *ljet)
{
  ljet[0]=-1; ljet[1]=-2;ljet[2]=-3;

  float tempt=-9;
  //! Get the leading jet
  for(int ij=0; ij<jetc->nref; ij++){
    if(fabs(jetc->jteta[ij])>2.0 || jetc->jtpt[ij]<30)continue;
    float jetpt = jetc->jtpt[ij];
    if(jetpt > tempt){
      tempt = jetpt;
      ljet[0] = ij;
    }
  }

  tempt=-9;
  for(int ij=0; ij<jetc->nref; ij++){
    if(ij==ljet[0])continue;
    if(fabs(jetc->jteta[ij])>2.0 || jetc->jtpt[ij]<30)continue;
    float jetpt = jetc->jtpt[ij];
    float dphi  = jetc->jtphi[ij] - jetc->jtphi[ljet[0]];
    if (dphi > pi ) dphi = dphi - 2 * pi;
    if (dphi < -pi) dphi = dphi + 2 * pi;
    if(dphi < 2*pi/3.)continue;
    if (jetpt > tempt){
      tempt = jetpt;
      ljet[1] = ij;
    }
  }

  tempt=-9;
  for(int ij=0; ij<jetc->nref; ij++){
    if(ij==ljet[0] || ij==ljet[1])continue;
    if(fabs(jetc->jteta[ij])>2.0|| jetc->jtpt[ij]<30)continue;
    float jetpt = jetc->jtpt[ij];
    if (jetpt > tempt){
      tempt = jetpt;
      ljet[2] = ij;
    }
  }
}
int GetMultBin(int nt)
{
  int ibin=-1;
  if(nt<60)ibin=0;
  else if(nt>=60 && nt<90)ibin=1;
  else if(nt>=90 && nt<110)ibin=2;
  else if(nt>=110&& nt<150)ibin=3;
  else if(nt>=150&& nt<180)ibin=4;
  else if(nt>=180)ibin=5;

  return ibin;
}
int GetHFplusEta4Bin(float nt)
{
    int ibin=-1;
    if(nt<20.)ibin=0;
    else if(nt>=20. && nt<25.)ibin=1;
    else if(nt>=25. && nt<30.)ibin=2;
    else if(nt>=30. && nt<40.)ibin=3;
    else if(nt>=40)ibin=4;
    
    return ibin;
}
