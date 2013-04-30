#include "TROOT.h"
#include "TString.h"
#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h" 
using namespace std;

void weight_PtHatBins_IncJet(const char *infile="input.root", const char *outfile = "output.root",float maxpthat=30., float xSection=1.,bool updateCorrections=true)
{

  //Reset ROOT and connect tree file
  gROOT->Reset();

  TFile *fin = TFile::Open(infile, "readonly");
  cout<<infile<<endl;
  bool IsMC = kFALSE ;
//  const Int_t ndir=10;
  const Int_t ndir=6;
  
  string inputDirName[ndir]= {
  //  "ak2PFJetAnalyzer",
    "ak3PFJetAnalyzer",
    "ak4PFJetAnalyzer",
    "ak5PFJetAnalyzer",
 //   "ak6PFJetAnalyzer",
 //   "ak2CaloJetAnalyzer",
    "ak3CaloJetAnalyzer",
    "ak4CaloJetAnalyzer",
    "ak5CaloJetAnalyzer"
 //   "ak6CaloJetAnalyzer"
  };

  string outputDirName[ndir]= {
  //  "ak2PFJetAnalyzer",
    "ak3PFJetAnalyzer",
    "ak4PFJetAnalyzer",
    "ak5PFJetAnalyzer",
 //   "ak6PFJetAnalyzer",
 //   "ak2CaloJetAnalyzer",
    "ak3CaloJetAnalyzer",
    "ak4CaloJetAnalyzer",
    "ak5CaloJetAnalyzer"
//    "ak6CaloJetAnalyzer"
  };

  string calgo[ndir]= {
//    "ak2PF",
    "ak3PF",
    "ak4PF",
    "ak5PF",
//    "ak6PF",
//    "ak2Calo",
    "ak3Calo",
    "ak4Calo",
    "ak5Calo"
 //   "ak6Calo"
  };
 
  string corrFileName[ndir]= {
 //   "AK2PF",
    "AK3PF",
    "AK4PF",
    "AK5PF",
//    "AK6PF",
//    "AK2Calo",
    "AK3Calo",
    "AK4Calo",
    "AK5Calo"
 //   "AK6Calo"
  };
  /*
  string corrFileName[ndir]= {
    "AK1",
    "AK2",
    "AK3",
    "AK4",
    "AK5",
    "AK6",
    "AK1",
    "AK2",
    "AK3",
    "AK4",
    "AK5",
    "AK6"
  };
*/

  TFile *fout=new TFile(outfile,"RECREATE");

  TTree *tr_in[ndir], *tr_out[ndir];
  
  TTree *tr_ev = (TTree*)fin->Get("hiEvtAnalyzer/HiTree");
  Int_t           hiBin;
  tr_ev->SetBranchAddress("hiBin",&hiBin);
  tr_ev->SetBranchStatus("*",0,0);
  tr_ev->SetBranchStatus("hiBin",1,0);
  //  float fentries = (float)tr_ev->GetEntries();

//  for(Int_t idir=1;idir<ndir;idir++){
  for(Int_t idir=0;idir<ndir;idir++){
    cout <<"idir =" <<idir <<" JetName ="<< inputDirName[idir].c_str() <<endl ;
    fin->cd(inputDirName[idir].c_str());
  //  fin->cd(inputDirName[idir].Data());
    tr_in[idir] = (TTree*)gDirectory->Get("t");
    
    //Declaration of leaves types
    int   nref;
    float pthat;
    float corrpt[100];
    float jtpt[100];
    float rawpt[100];
    float jteta[100];
    float jtphi[100];
    float jty[100];
    float jtpu[100];
    float refpt[100];
    float refeta[100];
    float refphi[100];
    float refdphijt[100];       
    float refdrjt[100];       
    float refparton_pt[100];
    int refparton_flavor[100];
    
    float fentries = (float)tr_in[idir]->GetEntries();
    float weight = 1.;
   if(IsMC) weight = xSection/(fentries/1000.);
    else weight = 1. ;
   // float weight = xSection/fentries;
    cout<<" weight "<<weight<<endl;
    
    tr_in[idir]->SetBranchAddress("nref",&nref);
    tr_in[idir]->SetBranchAddress("pthat",&pthat);
    tr_in[idir]->SetBranchAddress("rawpt",rawpt);
    tr_in[idir]->SetBranchAddress("jtpt",jtpt);
    tr_in[idir]->SetBranchAddress("jteta",jteta);
    tr_in[idir]->SetBranchAddress("jty",jty);
    tr_in[idir]->SetBranchAddress("jtphi",jtphi);
    tr_in[idir]->SetBranchAddress("jtpu",jtpu);
    tr_in[idir]->SetBranchAddress("refpt",refpt);
    tr_in[idir]->SetBranchAddress("refphi",refphi);
    tr_in[idir]->SetBranchAddress("refeta",refeta);
    tr_in[idir]->SetBranchAddress("refdphijt",refdphijt);
    tr_in[idir]->SetBranchAddress("refdrjt",refdrjt);
    tr_in[idir]->SetBranchAddress("refparton_pt",refparton_pt);
    tr_in[idir]->SetBranchAddress("refparton_flavor",refparton_flavor);
      cout<<"get jet trees!!! "<<endl;
    //! Add Friends to the TTree
    tr_in[idir]->AddFriend(tr_ev);
    fout->mkdir(outputDirName[idir].c_str());
    fout->cd(outputDirName[idir].c_str());
    
    tr_out[idir] = new TTree("t","Jet  Response Analyzer");
    tr_out[idir]->SetMaxTreeSize(4200000000);
    
    // Set output branch addresses.
    tr_out[idir]->Branch("hiBin",&hiBin,"hiBin/I");
    tr_out[idir]->Branch("nref",&nref,"nref/I");
    tr_out[idir]->Branch("pthat",&pthat,"pthat/F");
    tr_out[idir]->Branch("weight",&weight,"weight/F");
    
    tr_out[idir]->Branch("jtpt",jtpt,"jtpt[nref]/F");
    tr_out[idir]->Branch("rawpt",rawpt,"rawpt[nref]/F");
    tr_out[idir]->Branch("corrpt",corrpt,"corrpt[nref]/F");
    tr_out[idir]->Branch("jtpu",jtpu,"jtpu[nref]/F");
    tr_out[idir]->Branch("jteta",jteta,"jteta[nref]/F");
    tr_out[idir]->Branch("jty",jty,"jty[nref]/F");
    tr_out[idir]->Branch("jtphi",jtphi,"jtphi[nref]/F");
    tr_out[idir]->Branch("refpt",refpt,"refpt[nref]/F");
    tr_out[idir]->Branch("refeta",refeta,"refeta[nref]/F");
    tr_out[idir]->Branch("refphi",refphi,"refphi[nref]/F");
    tr_out[idir]->Branch("refdphijt",refdphijt,"refdphijt[nref]/F");
    tr_out[idir]->Branch("refdrjt",refdrjt,"refdrjt[nref]/F");
    tr_out[idir]->Branch("refparton_pt",refparton_pt,"refparton_pt[nref]/F");
    tr_out[idir]->Branch("refparton_flavor",refparton_flavor,"refparton_flavor[nref]/I");
    //! grab the JEC's
    string L2Name, L3Name;
    JetCorrectorParameters* parHI44x_l2, * parHI44x_l3;
    vector<JetCorrectorParameters> vpar_HI44x;   
  //  FactorizedJetCorrector *_JEC_HI44X;
    FactorizedJetCorrector *_JEC_HI44X = new FactorizedJetCorrector(vpar_HI44x);//JR
  //FactorizedJetCorrector *_JEC_HI44X= new JetCorrectorParameters("YYYY.txt"); 
    if(updateCorrections){   

     // L2Name = "/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/JEC_dijet_L2Relative_"+corrFileName[idir]+".txt";
    //  L2Name = "/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/JEC_"+calgo[idir]+"_dijet_L2Relative_"+corrFileName[idir]+".txt";
      L2Name = "/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/txt/pp276/JEC_dijet_L2Relative_"+corrFileName[idir]+".txt";
      cout<<"**** ++++++  L2Name "<<L2Name<<endl;
     // L3Name = "/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/JEC_dijet_L3Absolute_"+corrFileName[idir]+".txt";
      L3Name = "/net/hidsk0001/d00/scratch/maoyx/pPb/CMSSW_5_3_8_HI_patch2/src/MNguyen/JEC/txt/pp276/JEC_dijet_L3Absolute_"+corrFileName[idir]+".txt";
      cout<<"**** ++++++  L3Name "<<L3Name<<endl;

      cout<<" L2Name "<<L2Name<<endl;
      
      parHI44x_l2 = new JetCorrectorParameters(L2Name.c_str());
      parHI44x_l3 = new JetCorrectorParameters(L3Name.c_str());
      
      vpar_HI44x.push_back(*parHI44x_l2);
      vpar_HI44x.push_back(*parHI44x_l3);
      _JEC_HI44X = new FactorizedJetCorrector(vpar_HI44x);
    }     
    
    Long64_t nentries = tr_in[idir]->GetEntries();
    Long64_t nbytes = 0;
    for (Long64_t i=0; i<nentries;i++) {
      nbytes += tr_in[idir]->GetEntry(i);
      if(pthat>maxpthat)continue;	 
      
      //! jet loop
      for(int iref=0;iref<nref;iref++){
	corrpt[iref] = jtpt[iref];
	jtpt[iref]   = rawpt[iref];
	
	if(updateCorrections){
	  _JEC_HI44X->setJetEta(jteta[iref]);
	  _JEC_HI44X->setJetPt(rawpt[iref]);
	  corrpt[iref] = rawpt[iref]*_JEC_HI44X->getCorrection(); 
      //  if(TMath::Abs(corrpt[iref]-(rawpt[iref]*_JEC_HI44X->getCorrection()))>1. && rawpt[iref]>10.) cout<<"*****!!!!!"<<endl;
      //    if(TMath::Abs(_JEC_HI44X->getCorrection()-(corrpt[iref]/rawpt[iref]))>0.2 && rawpt[iref]>10.) cout<<"iref: "<<iref<<" correction is: "<<(_JEC_HI44X->getCorrection())<<" ratio of corr/raw: "<<(corrpt[iref]/rawpt[iref])<<" for rawpt: "<<rawpt[iref]<<" corrpt: "<<corrpt[iref]<<" and rawpT*correction: "<<(rawpt[iref]*_JEC_HI44X->getCorrection())<<" and correction*rawpt: "<<(_JEC_HI44X->getCorrection()*rawpt[iref])<<" and eta: "<<jteta[iref]<<endl;
	}
      }//ref for loop
      tr_out[idir]->Fill();
    }//nentresi loop
  tr_out[idir]->Write();
    cout <<"finish the writing"<<endl;
    if(updateCorrections){
      delete parHI44x_l2;
      delete parHI44x_l3;
      delete _JEC_HI44X;
    }
  //  if (updateCorrections==true)delete _JEC_HI44X;
  }
   cout <<"finish the code!!"<<endl ;
  fout->Close();
   cout<<"XXXX"<<endl;
}
