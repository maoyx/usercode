#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF1.h"
#include <TH1D.h>
#include <TNtuple.h>
#include "TChain.h"
#include <iostream>
#include <TString.h>
#include <TCut.h>
#include "../SmearingFactors.h"
#include "hiForest.h"
#include "commonSetup.h"
#include <stdlib.h>
using namespace std;

//define the kinematics cuts for the analysis
const double conesize = 0.3 ;
const double deltacone = 0.05 ;

//const double radius = 0.025;

//const int nrbin = 10 ;
//double rbin[nrbin+1]={0.0,0.05,0.10, 0.15, 0.20,0.25,0.30,0.35,0.40,0.45,0.50};

const int nrbin = 6 ;
double rbin[nrbin+1]={0.0,0.05,0.10, 0.15, 0.20,0.25,0.30};

const double etacut = 2.0 ;
//const double dphicut = TMath::TwoPi()/3. ;
const double dphicut = 7*TMath::Pi()/8. ;
const double leadingjetcut = 100. ;
const double trackcut = 1.;
const double tracketacut = 2.4; //make sure the tracks around jet cone inside eta range

////for etalimit = 0.3, to remove the possible overlap region by ER method 
const double etalimit = 0.3 ; 

int pthat ; //= 120 ; //=300 ; //30 and 80 for pp; 30,50,80,120,170,200 for PbPb MC 
int ptmax ; //= 170 ;
TString coll = "PP";
//if it is pp, no centrality bins, only one
const int nbin = 1 ;
const int centr[] ={0,100};
////for HI centrality bin
//const int nbin = 4 ;
//const int centr[] ={0,10,30,50,100};

TString JetChoose ="Ref" ; //setting to analyze Gen jet or Ref jet in hiForest
TString TrkChoose ="Gen" ; //setting to analyze SimTrk or hiGenParticle 

const bool DoSumPtNorm = kFALSE ;

TString intputFile ;

//TString dataPath="/Users/ymao/group/CMS/hiForest";
TString dataPath ;

const int nptbin = 1 ;
const double pt[]={100., 300.};
//const int ntrkptbin = 6 ;
//const double trkpt[]={1., 2., 4., 8., 16., 32., 500.};

const int ntrkptbin = 3 ;
const double trkpt[]={1.,4., 16., 300.};

const Double_t jetPtBin[] = {30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,400};
//const Double_t jetPtBin[] = {100, 110, 120, 130, 140, 150, 160, 170, 180, 200, 240, 300};
//const Double_t jetPtBin[] = {100, 110, 120, 130, 140, 150, 160, 180, 200, 240, 300, 500};
const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;

const int nPtBin = 27;
double TrkBin[nPtBin+1] = {0, 0.5, 1, 1.203915, 1.449412, 1.74497, 2.100796, 2.529181, 3.04492, 3.665826, 4.413344, 5.313293, 6.396755, 7.701152, 9.271536, 11.16214, 13.43828, 16.17855, 19.47761, 23.44939, 28.23108, 33.98783, 40.91848, 49.26238, 59.30774, 71.40151, 85.96137, 103.4902}; 

const int netabin = 4 ;
const double deta[]={0.0,0.5,1.0,1.5,2.0} ;

class hist_class {
public:
    hist_class(TString the_desc);
    void Delete();
    void Write();
    
    TString desc;
    
    TH1F * NEvents[nbin];
    TH1F * NJets[nbin];
    TH1F * Ntrack[nbin];
    TH1F * jetpt[nbin];
    TH1F * qjetpt[nbin]; //flavor = 1-6
    TH1F * gjetpt[nbin]; //flavor = 21
    TH2F * jeteta[nbin];
    TH2F * jetphi[nbin];
    TH2F * trackpt[nbin];
    TH2F * trackphi[nbin];
    TH2F * tracketa[nbin];
    TH2F * trackdr[nbin][nptbin];
    TH2F * bkgtrackdr[nbin][nptbin];
    TH2F * jetaxisRes[nbin];
    
    TH2F * jetfrag[nbin] ;
    TH2F * jetbkgfrag[nbin] ;
    TH2F * SumptJetPt[nbin] ;
    TH2F * JetBkgCone[nbin] ;
    TH2F * IncTrackDphi[nbin];
    TH2F * JetTrackDphi[nbin];
    TH2F * JetBkgTrackDphi[nbin];
    TH2F * JetTrackDeta[nbin];
    TH2F * JetBkgTrackDeta[nbin];
    TH2F * bkgtrackpt[nbin];

    TH2F * jetdrSumTrkPt[nbin][nptbin][ntrkptbin];

    TH2F * jetPtSumTrk[nbin][ntrkptbin];

    TH2F * bkgjetdrSumTrkPt[nbin][nptbin][ntrkptbin];
    TH2F * bkgjetPtSumTrk[nbin][ntrkptbin];

    TH1F * deltaR;
    TH1D * CenBin;
    TH1D * RareEvt;
    TH1F * Vertex ;

    //   TH1F * ptbin[6];
    //for jet shape variables
    TH2F * ChargePt[nbin][nrbin];
    TH2F * SumChPt[nbin][nrbin];
    TH2F * ChargeMult[nbin][nrbin];   
    TH2F * DiffJS[nbin][nrbin]; //differential jet shapes (pho(r)) hist.
    TH2F * IntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    //  TH2F * dphiJS[nbin]; //phi diviration from jet phi 
    //   TH2F * detaJS[nbin] ; //eta diviration from jet eta 
    TH1F * DiffJSPt80_100[nbin]; //differential jet shapes (pho(r)) hist.
    TH1F * bkgDiffJSPt80_100[nbin]; //differential jet shapes (pho(r)) hist.

    TH1F * DiffJSPt80_100Rbin[nbin][nrbin]; //differential jet shapes (pho(r)) hist.
    TH1F * bkgDiffJSPt80_100Rbin[nbin][nrbin]; //differential jet shapes (pho(r)) hist.

    //for jet shape background study
    TH2F * bkgChargePt[nbin][nrbin];
    TH2F * bkgSumChPt[nbin][nrbin];
    TH2F * bkgChargeMult[nbin][nrbin];   
    TH2F * bkgDiffJS[nbin][nrbin]; //differential jet shapes (pho(r)) hist.
    TH2F * bkgIntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.

    TH2F * GenAxisDiffJS[nbin][nrbin]; //differential jet shapes (pho(r)) hist.
    TH2F * GenAxisbkgDiffJS[nbin][nrbin]; //differential jet shapes (pho(r)) hist.
    
    TH2F * DiffEtaBinJS[nbin][netabin][nrbin]; //differential jet shapes (pho(r)) hist.
    TH2F * bkgDiffEtaBinJS[nbin][netabin][nrbin]; //differential jet shapes (pho(r)) hist.
};

hist_class::hist_class(TString the_desc)
{
    
    desc = the_desc;
    for(int ibin = 0 ; ibin <nbin; ibin++){
        NEvents[ibin] = new TH1F(Form("Nevents_%d-%d%%",centr[ibin],centr[ibin+1]), Form("Nevents_%d-%d%%",centr[ibin],centr[ibin+1]), 100, 0, 2.);
        NJets[ibin] = new TH1F(Form("NJets_%d-%d%%",centr[ibin],centr[ibin+1]), Form("NJets_%d-%d%%",centr[ibin],centr[ibin+1]), 100, -0.5, 99.5);
        Ntrack[ibin] = new TH1F(Form("Ntracks_%d-%d%%",centr[ibin],centr[ibin+1]), Form("Ntracks_%d-%d%%",centr[ibin],centr[ibin+1]), 800, -1., 799);
        jetpt[ibin] = new TH1F(Form("jetpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500.);
        jetpt[ibin]->Sumw2();
        qjetpt[ibin] = new TH1F(Form("quarkjetpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("quarkjetpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500.);
        qjetpt[ibin]->Sumw2();
        gjetpt[ibin] = new TH1F(Form("gluonjetpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("gluonjetpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500.);
        gjetpt[ibin]->Sumw2();
        jeteta[ibin] = new TH2F(Form("jeteta_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jeteta_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500,100, -5.05, 4.95);
        jeteta[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        jeteta[ibin]->GetYaxis()->SetTitle("#eta");
        jeteta[ibin]->Sumw2();
        jetphi[ibin] = new TH2F(Form("jetphi_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetphi_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500,100, -5.05, 4.95);
        jetphi[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        jetphi[ibin]->GetYaxis()->SetTitle("#phi");
        jetphi[ibin]->Sumw2();
        trackpt[ibin] = new TH2F(Form("trackpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("trackpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, nPtBin, TrkBin); 
//        trackpt[ibin] = new TH2F(Form("trackpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("trackpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, 200, 0., 200); 
        trackpt[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        trackpt[ibin]->GetYaxis()->SetTitle("p_{T}^{trk} (GeV/c)");   
        trackpt[ibin]->Sumw2();
        trackphi[ibin] = new TH2F(Form("trackphi_%d-%d%%",centr[ibin],centr[ibin+1]), Form("trackphi_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, 100, -5.05, 4.95); 
        trackphi[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        trackphi[ibin]->GetYaxis()->SetTitle("#phi^{trk}");   
        trackphi[ibin]->Sumw2();
        tracketa[ibin] = new TH2F(Form("tracketa_%d-%d%%",centr[ibin],centr[ibin+1]), Form("tracketa_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, 100, -5.05, 4.95); 
        tracketa[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        tracketa[ibin]->GetYaxis()->SetTitle("#eta^{trk}");   
        tracketa[ibin]->Sumw2();

        jetfrag[ibin] = new TH2F(Form("FFleading_%d-%d%%",centr[ibin],centr[ibin+1]), Form("FFleading_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, 100, 0., 10.);
        jetfrag[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        jetfrag[ibin]->GetYaxis()->SetTitle("#xi = ln(p_{T}^{jet}/p_{T}^{h})");
        jetfrag[ibin]->Sumw2();
        jetbkgfrag[ibin] = new TH2F(Form("bkgFFleading_%d-%d%%",centr[ibin],centr[ibin+1]), Form("bkgFFleading_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, 100, 0., 10.);
        jetbkgfrag[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        jetbkgfrag[ibin]->GetYaxis()->SetTitle("#xi = ln(p_{T}^{jet}/p_{T}^{h})");
        jetbkgfrag[ibin]->Sumw2();

        for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
            trackdr[ibin][ipt] = new TH2F(Form("recoTrackJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("recoTrackJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize/deltacone), 0., conesize); 
            trackdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{trk} (GeV/c)");
            trackdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");   
            trackdr[ibin][ipt]->Sumw2();
            bkgtrackdr[ibin][ipt] = new TH2F(Form("recobkgTrackJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("recobkgTrackJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize/deltacone), 0., conesize);
            bkgtrackdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{trk} (GeV/c)");
            bkgtrackdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");
            bkgtrackdr[ibin][ipt]->Sumw2();
            for(Int_t itr = 0 ; itr <ntrkptbin ; itr++){
                jetdrSumTrkPt[ibin][ipt][itr]=new TH2F(Form("JetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1],trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), Form("JetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1],trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), (Int_t)(conesize/deltacone), 0., conesize, 100, 0., 1.);
                jetdrSumTrkPt[ibin][ipt][itr]->GetXaxis()->SetTitle("radius r"); 
                jetdrSumTrkPt[ibin][ipt][itr]->GetYaxis()->SetTitle("#Sigma p_{T}^{trk}/p_{T}^{jet}"); 
                jetdrSumTrkPt[ibin][ipt][itr]->Sumw2();
   
                bkgjetdrSumTrkPt[ibin][ipt][itr]=new TH2F(Form("bkgJetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1],trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), Form("bkgJetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1],trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), (Int_t)(conesize/deltacone), 0., conesize, 100, 0., 1.);
                bkgjetdrSumTrkPt[ibin][ipt][itr]->GetXaxis()->SetTitle("radius r"); 
                bkgjetdrSumTrkPt[ibin][ipt][itr]->GetYaxis()->SetTitle("#Sigma p_{T}^{trk}/p_{T}^{jet}"); 
                bkgjetdrSumTrkPt[ibin][ipt][itr]->Sumw2();
            }
        }
        SumptJetPt[ibin] = new TH2F(Form("SumptJetPtRatio_%d-%d%%",centr[ibin],centr[ibin+1]), Form("SumptJetPtRatio_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500.,100, 0., 1.);  
        SumptJetPt[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        SumptJetPt[ibin]->GetYaxis()->SetTitle("#Sigma p_{T}^{trk}/p_{T}^{jet}");
        SumptJetPt[ibin]->Sumw2();
        JetBkgCone[ibin] = new TH2F(Form("JetBkgAxisDifference_%d-%d%%",centr[ibin],centr[ibin+1]), Form("JetBkgAxisDifference_%d-%d%%",centr[ibin],centr[ibin+1]), 140, -7., 7.,100, -5.05, 4.95);
        JetBkgCone[ibin]->GetXaxis()->SetTitle("#phi_{bkg}-#phi_{jet}");
        JetBkgCone[ibin]->GetYaxis()->SetTitle("#eta_{bkg}-#eta_{jet}");
        JetBkgCone[ibin]->Sumw2();
        IncTrackDphi[ibin] = new TH2F(Form("AllTrackDphi_%d-%d%%",centr[ibin],centr[ibin+1]), Form("AllTrackDphi_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500,120, -4.05, 7.95);
        IncTrackDphi[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        IncTrackDphi[ibin]->GetYaxis()->SetTitle("#phi_{trk}-#phi_{jet}");
        IncTrackDphi[ibin]->Sumw2();
        JetTrackDphi[ibin] = new TH2F(Form("JetTrackDphi_%d-%d%%",centr[ibin],centr[ibin+1]), Form("JetTrackDphi_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500,120, -4.05, 7.95);
        JetTrackDphi[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        JetTrackDphi[ibin]->GetYaxis()->SetTitle("#phi_{trk}-#phi_{jet}");
        JetTrackDphi[ibin]->Sumw2();
        JetBkgTrackDphi[ibin] = new TH2F(Form("JetBkgTrackDphi_%d-%d%%",centr[ibin],centr[ibin+1]), Form("JetBkgTrackDphi_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500,120, -4.05, 7.95);
        JetBkgTrackDphi[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        JetBkgTrackDphi[ibin]->GetYaxis()->SetTitle("#phi_{trk}-#phi_{jet}");
        JetBkgTrackDphi[ibin]->Sumw2();
        JetTrackDeta[ibin] = new TH2F(Form("JetTrackDeta_%d-%d%%",centr[ibin],centr[ibin+1]), Form("JetTrackDeta_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500,100, -5.05, 4.95);
        JetTrackDeta[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        JetTrackDeta[ibin]->GetYaxis()->SetTitle("#eta_{trk}-#eta_{jet}");
        JetTrackDeta[ibin]->Sumw2();
        JetBkgTrackDeta[ibin] = new TH2F(Form("JetBkgTrackDeta_%d-%d%%",centr[ibin],centr[ibin+1]), Form("JetBkgTrackDeta_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500,100, -5.05, 4.95);
        JetBkgTrackDeta[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        JetBkgTrackDeta[ibin]->GetYaxis()->SetTitle("#eta_{trk}-#eta_{jet}");
        JetBkgTrackDeta[ibin]->Sumw2();
        bkgtrackpt[ibin] = new TH2F(Form("bkgtrackpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("bkgtrackpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, nPtBin, TrkBin); 
//        bkgtrackpt[ibin] = new TH2F(Form("bkgtrackpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("bkgtrackpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, 200, 0., 200); 
        bkgtrackpt[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        bkgtrackpt[ibin]->GetYaxis()->SetTitle("p_{T}^{trk} (GeV/c)");   
        bkgtrackpt[ibin]->Sumw2();
        for(Int_t itr = 0 ; itr <ntrkptbin ; itr++){
            jetPtSumTrk[ibin][itr]=new TH2F(Form("JetSumTrkPt%.f_%.f_%d-%d%%",trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), Form("JetSumTrkPt%.f_%.f_%d-%d%%",trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), 500, 0., 500, 500, 0., 500.);
            jetPtSumTrk[ibin][itr]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)"); 
            jetPtSumTrk[ibin][itr]->GetYaxis()->SetTitle("#Sigma p_{T}^{trk}"); 
            jetPtSumTrk[ibin][itr]->Sumw2();
            bkgjetPtSumTrk[ibin][itr]=new TH2F(Form("bkgJetSumTrkPt%.f_%.f_%d-%d%%",trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), Form("bkgJetSumTrkPt%.f_%.f_%d-%d%%",trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), 500, 0., 500, 500, 0., 500.);
            bkgjetPtSumTrk[ibin][itr]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)"); 
            bkgjetPtSumTrk[ibin][itr]->GetYaxis()->SetTitle("#Sigma p_{T}^{trk}"); 
            bkgjetPtSumTrk[ibin][itr]->Sumw2();
 
        }

        DiffJSPt80_100[ibin] = new TH1F(Form("diffJSPt80_100_Cen%d-%d%%",centr[ibin],centr[ibin+1]), Form("diffJSPt80_100JS_Cen%d-%d%%",centr[ibin],centr[ibin+1]), 8, -0.05, 0.35);
        DiffJSPt80_100[ibin]->GetXaxis()->SetTitle("#rho (r)");
        DiffJSPt80_100[ibin]->Sumw2();
        bkgDiffJSPt80_100[ibin] = new TH1F(Form("bkgdiffJSPt80_100_Cen%d-%d%%",centr[ibin],centr[ibin+1]), Form("bkgdiffJSPt80_100JS_Cen%d-%d%%",centr[ibin],centr[ibin+1]), 8, -0.05, 0.35);
        bkgDiffJSPt80_100[ibin]->GetXaxis()->SetTitle("#rho (r)");
        bkgDiffJSPt80_100[ibin]->Sumw2();

        for(int ir = 0 ; ir <nrbin; ir++){
            ChargePt[ibin][ir] = new TH2F(Form("chargeptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("chargeptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 500, 0., 500);
            ChargePt[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            ChargePt[ibin][ir]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
            ChargePt[ibin][ir]->Sumw2();
            SumChPt[ibin][ir] = new TH2F(Form("sumchptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("sumchargeptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
            SumChPt[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            SumChPt[ibin][ir]->GetYaxis()->SetTitle("#Sigma p_{T}^{h^{#pm}}/p_{T}^{jet}");
            SumChPt[ibin][ir]->Sumw2();
            ChargeMult[ibin][ir] = new TH2F(Form("ChMultdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("ChMultdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 100, -0.5, 99.5);
            ChargeMult[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            ChargeMult[ibin][ir]->GetYaxis()->SetTitle("# of Charge");   
            ChargeMult[ibin][ir]->Sumw2();
            DiffJS[ibin][ir] = new TH2F(Form("differentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("differentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
            DiffJS[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            DiffJS[ibin][ir]->GetYaxis()->SetTitle("#rho (r)");
            DiffJS[ibin][ir]->Sumw2();
            IntJS[ibin][ir] = new TH2F(Form("IntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("IntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
            IntJS[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            IntJS[ibin][ir]->GetYaxis()->SetTitle(Form("#psi(r=%1.2f)",deltacone*(ir+1)));
            IntJS[ibin][ir]->Sumw2();
            DiffJSPt80_100Rbin[ibin][ir] = new TH1F(Form("diffJSPt80_100dR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("diffJSPt80_100JSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 150, -0.005, 1.495);
            DiffJSPt80_100Rbin[ibin][ir]->GetXaxis()->SetTitle("#rho (r)");
            DiffJSPt80_100Rbin[ibin][ir]->Sumw2();
            bkgDiffJSPt80_100Rbin[ibin][ir] = new TH1F(Form("bkgdiffJSPt80_100dR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("bkgdiffJSPt80_100JSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 150, -0.005, 1.495);
            bkgDiffJSPt80_100Rbin[ibin][ir]->GetXaxis()->SetTitle("#rho (r)");
            bkgDiffJSPt80_100Rbin[ibin][ir]->Sumw2();
            
            //for bkg histos.
            bkgChargePt[ibin][ir] = new TH2F(Form("bkgchargeptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("bkgchargeptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 500, 0., 500);
            bkgChargePt[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            bkgChargePt[ibin][ir]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
            bkgChargePt[ibin][ir]->Sumw2();
            bkgSumChPt[ibin][ir] = new TH2F(Form("bkgsumchptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("bkgsumchargeptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
            bkgSumChPt[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            bkgSumChPt[ibin][ir]->GetYaxis()->SetTitle("#Sigma p_{T}^{h^{#pm}}/p_{T}^{jet}");    
            bkgSumChPt[ibin][ir]->Sumw2();
            bkgChargeMult[ibin][ir] = new TH2F(Form("bkgChMultdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("bkgChMultdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 100, -0.5, 99.5);
            bkgChargeMult[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            bkgChargeMult[ibin][ir]->GetYaxis()->SetTitle("# of Charge");
            bkgChargeMult[ibin][ir]->Sumw2();
            bkgDiffJS[ibin][ir] = new TH2F(Form("bkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("bkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
            bkgDiffJS[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            bkgDiffJS[ibin][ir]->GetYaxis()->SetTitle("#rho (r)");
            bkgDiffJS[ibin][ir]->Sumw2();
            bkgIntJS[ibin][ir] = new TH2F(Form("bkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("bkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
            bkgIntJS[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            bkgIntJS[ibin][ir]->GetYaxis()->SetTitle(Form("#psi(r=%1.2f)",deltacone*(ir+1)));
            bkgIntJS[ibin][ir]->Sumw2();
        }
            for(Int_t ieta = 0 ; ieta <netabin ; ieta++){
                for(int ir = 0 ; ir <nrbin; ir++){
                DiffEtaBinJS[ibin][ieta][ir] = new TH2F(Form("diffJSdEtaBin%.f_%.fdR%.f_%.f_Cen%d-%d%%", deta[ieta]*10,deta[ieta+1]*10, rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("diffJSdEtaBin%.f_%.fdR%.f_%.f_Cen%d-%d%%", deta[ieta]*10,deta[ieta+1]*10,rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
                DiffEtaBinJS[ibin][ieta][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                DiffEtaBinJS[ibin][ieta][ir]->GetYaxis()->SetTitle("#rho (r)");
                DiffEtaBinJS[ibin][ieta][ir]->Sumw2();
                bkgDiffEtaBinJS[ibin][ieta][ir] = new TH2F(Form("bkgdiffJSdEtaBin%.f_%.fdR%.f_%.f_Cen%d-%d%%", deta[ieta]*10,deta[ieta+1]*10, rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("bkgdiffJSdEtaBin%.f_%.fdR%.f_%.f_Cen%d-%d%%", deta[ieta]*10,deta[ieta+1]*10,rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
                bkgDiffEtaBinJS[ibin][ieta][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                bkgDiffEtaBinJS[ibin][ieta][ir]->GetYaxis()->SetTitle("#rho (r)");
                bkgDiffEtaBinJS[ibin][ieta][ir]->Sumw2();
                }
  
            }
        }  //centrality bins loop
        
    //    for(int i = 0 ; i < 6 ; i++){
    //     ptbin[i] = new TH1F(Form("ptbin_%.f-%.f",pt[i], pt[i+1]), Form("ptbin_%.f-%.f",pt[i], pt[i+1]), 500, 0., 500.);
    //    }
    // if(!deltaR)
    deltaR = new TH1F((TString) (desc + "_deltaR"), "", 100, 0., 10.);
    CenBin = new TH1D((TString) (desc + "_Cent"), "", 40, 0, 40);
    CenBin->Sumw2();
    Vertex = new TH1F((TString) (desc + "_Vz"), "", 400, -20., 20.);
    Vertex->Sumw2();
    RareEvt = new TH1D((TString) (desc + "_RareEvt"), "", 40, 0, 40);
}

void hist_class::Delete()
{
    for(int ibin = 0 ; ibin <nbin; ibin++){
        delete NEvents[ibin];
        //        delete TwoTrackTwoJet[ibin];
        //        delete TwoTrackThreeJet[ibin];
        //        delete OneTrackTwoJet[ibin];
        //        delete OneTrackThreeJet[ibin];
        //        delete SumPtFraction[ibin];
        delete NJets[ibin];
        delete Ntrack[ibin];  
        delete jetpt[ibin];
        delete qjetpt[ibin];
        delete gjetpt[ibin];
        delete jeteta[ibin];
        delete jetphi[ibin];
        delete trackpt[ibin];
        delete trackphi[ibin];
        delete tracketa[ibin];
        
        delete jetfrag[ibin];
        delete jetbkgfrag[ibin];
        for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
            for(Int_t itr = 0 ; itr <ntrkptbin ; itr++){
             delete trackdr[ibin][ipt];
             delete bkgtrackdr[ibin][ipt];
              delete jetdrSumTrkPt[ibin][ipt][itr];
                delete bkgjetdrSumTrkPt[ibin][ipt][itr];
                
            }
        }
        delete SumptJetPt[ibin];
        delete JetBkgCone[ibin];
        delete IncTrackDphi[ibin];
        delete JetTrackDphi[ibin];
        delete JetBkgTrackDphi[ibin];
        delete JetTrackDeta[ibin];
        delete JetBkgTrackDeta[ibin];
        delete bkgtrackpt[ibin];
        for(Int_t itr = 0 ; itr <ntrkptbin ; itr++){
            delete jetPtSumTrk[ibin][itr];
            delete bkgjetPtSumTrk[ibin][itr];
            
        }
        delete DiffJSPt80_100[ibin];        
        delete bkgDiffJSPt80_100[ibin];  
        for(int ir = 0 ; ir <nrbin; ir++){
            delete ChargePt[ibin][ir];
            delete SumChPt[ibin][ir];
            delete ChargeMult[ibin][ir];
            delete DiffJS[ibin][ir];
            delete IntJS[ibin][ir];
            delete bkgChargePt[ibin][ir];      
            delete bkgSumChPt[ibin][ir];
            delete bkgChargeMult[ibin][ir];
            delete bkgDiffJS[ibin][ir];
            delete bkgIntJS[ibin][ir];
            delete DiffJSPt80_100Rbin[ibin][ir];        
            delete bkgDiffJSPt80_100Rbin[ibin][ir];  
         }
        for(Int_t ieta = 0 ; ieta <netabin ; ieta++){
            for(int ir = 0 ; ir <nrbin; ir++){
                delete DiffEtaBinJS[ibin][ieta][ir];
                delete bkgDiffEtaBinJS[ibin][ieta][ir];
            }
        }
    } //centrality loop
    //   for(int i=0; i<6;i++) delete ptbin[i];
    delete deltaR;
    delete CenBin;
    delete RareEvt;
    delete Vertex;
}

void hist_class::Write()
{
    TString out_name ;
    TString dataType ; 
    dataType="MC" ;
    TString met ;
    if(etalimit==etacut) met = "RC";
    else met = "ER";
    TString Norm ;
    if(DoSumPtNorm)Norm="NormSum";
    else Norm ="NormJet";
    out_name=Form("%s%s_%sJetPt%.f_%sTrk%.fEtaCut%.fLimit%.f_%sbkgJS%sCone%.f_CenBin%d_Nrbin%d_LJbin%dpthat%d_%s",dataType.Data(),coll.Data(),JetChoose.Data(), leadingjetcut,TrkChoose.Data(),trackcut,etacut*10, etalimit*10,met.Data(), Norm.Data(),conesize*10, nbin,nrbin,nptbin, pthat, intputFile.Data());       
 
       TFile *out_file = new TFile(Form("/net/hidsk0001/d00/scratch/maoyx/JetShape/CMSSW_4_4_2_patch5/src/HiForest_V2_02_16/%s",out_name.Data()), "RECREATE");
    for(int ibin = 0 ; ibin <nbin; ibin++){
        NEvents[ibin]->Write();
        NJets[ibin]->Write();
        Ntrack[ibin]->Write();
        jetpt[ibin]->Write();
        qjetpt[ibin]->Write();
        gjetpt[ibin]->Write(); 
        jeteta[ibin]->Write();
        jetphi[ibin]->Write();
        trackpt[ibin]->Write();
        trackphi[ibin]->Write();
        tracketa[ibin]->Write();
        
        jetfrag[ibin]->Write();
        jetbkgfrag[ibin]->Write();
        for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
            trackdr[ibin][ipt]->Write();
            bkgtrackdr[ibin][ipt]->Write();
            for(Int_t itr = 0 ; itr <ntrkptbin ; itr++){
                jetdrSumTrkPt[ibin][ipt][itr]->Write();
                 bkgjetdrSumTrkPt[ibin][ipt][itr]->Write();
                
            }

        }
        SumptJetPt[ibin]->Write();
        JetBkgCone[ibin]->Write();
        IncTrackDphi[ibin]->Write();
        JetTrackDphi[ibin]->Write();
        JetBkgTrackDphi[ibin]->Write();
        JetTrackDeta[ibin]->Write();
        JetBkgTrackDeta[ibin]->Write();
        bkgtrackpt[ibin]->Write();
        for(Int_t itr = 0 ; itr <ntrkptbin ; itr++){
            jetPtSumTrk[ibin][itr]->Write();
            bkgjetPtSumTrk[ibin][itr]->Write();
            
        }

        DiffJSPt80_100[ibin]->Write();
        bkgDiffJSPt80_100[ibin]->Write();
        for(int ir = 0 ; ir <nrbin; ir++){
            ChargePt[ibin][ir]->Write();    
            SumChPt[ibin][ir]->Write(); 
            ChargeMult[ibin][ir]->Write();    
            DiffJS[ibin][ir]->Write();
            IntJS[ibin][ir]->Write();
            bkgChargePt[ibin][ir]->Write();    
            bkgSumChPt[ibin][ir]->Write(); 
            bkgChargeMult[ibin][ir]->Write();   
            bkgDiffJS[ibin][ir]->Write();
            bkgIntJS[ibin][ir]->Write();
            DiffJSPt80_100Rbin[ibin][ir]->Write();
            bkgDiffJSPt80_100Rbin[ibin][ir]->Write();
    
        }
            for(Int_t ieta = 0 ; ieta <netabin ; ieta++){
                for(int ir = 0 ; ir <nrbin; ir++){
                DiffEtaBinJS[ibin][ieta][ir]->Write();
                bkgDiffEtaBinJS[ibin][ieta][ir]->Write();
                }
            }
        } //centrality bins
    
    //  for(int i=0; i<6;i++) ptbin[i]->Write();
    deltaR->Write();
    CenBin->Write();
    RareEvt->Write();
    Vertex->Write();
    out_file->Close();
    cout <<"Output file: " <<Form("%s",out_name.Data()) <<endl ;
    
}


void anaGenJetJS()
{
    hist_class *my_hists = new hist_class("pfjet");

    pthat=atoi(getenv("PTHAT")) ;
    ptmax=atoi(getenv("PTMAX")) ;
    cout <<"pthat = " <<pthat <<"  pthatmax =" <<ptmax <<endl ;

    std::cout << "start working\n";
    //         dataPath= "/Users/ymao/group/CMS/hiForest"; //local analysis
    if(coll=="HI"){ 
            if(pthat==50||pthat==80||pthat==100||pthat==170)
                dataPath= Form("/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/"); //MIT MC normial
            else 
                dataPath= Form("/mnt/hadoop/cms/store/user/yenjie/HiForest_v28/"); //MIT MC normial
    }
    else    
        dataPath= Form("/net/hisrv0001/home/zhukova/scratch/HIHighPt/forest/pthat%d", pthat); //lxplus path for pp

    //MC sample
    if(coll=="HI"){             
        //  intputFile=Form("pyquen%s%d_HYDJET.root", para.Data(), pthat);
        if(pthat==50||pthat==80||pthat==100||pthat==170)
            intputFile=Form("Dijet%d_HydjetDrum_v27_mergedV1.root", pthat);
        else 
            intputFile=Form("Dijet%d_HydjetDrum_v28_mergedV1.root", pthat);
        /*	        if(pthat==80) 
         //  intputFile=Form("Pythia%d_HydjetDrum_mix01_HiForest2_v22_simTrack05.root", pthat);
         intputFile=Form("Pythia%d_HydjetDrum_mix01_HiForest2_v20.root", pthat);
         else if(pthat==30||pthat==50 ||pthat==170)
         intputFile=Form("Pythia%d_HydjetDrum_mix01_HiForest2_v19.root", pthat);  
         else 
         intputFile=Form("Pythia%d_HydjetDrum_mix01_HiForest2_v21_ivan.root", pthat);
         //  intputFile="merged_HydjetDrum.root";
         */  
    }
    else 
        //   intputFile=Form("pp276Dijet%d_merged.root",pthat);  
        intputFile=Form("mergedFile.root");  
    
    TString inname=Form("%s/%s", dataPath.Data(),intputFile.Data());
    // Define the input file and HiForest
    bool isPP =false ;
    if(coll=="PP") isPP =true ;
    HiForest * c = new HiForest(inname,"forest",isPP);
    c->doTrackCorrections = false;
    c->doTrackingSeparateLeadingSubleading = false;
    c->InitTree();
  //  cout << "start working222222\n";
  //  TFile *my_file=TFile::Open(Form("%s/%s", dataPath.Data(),intputFile.Data()));
    cout <<"Input file" << inname<<endl ;
    
//    if(my_file->IsZombie()) {
//    }
    
    // evt tree
//    TTree *evt_tree = (TTree*) c->inf->Get("hiEvtAnalyzer/HiTree");
//    Evts offSel;
//    if (evt_tree) {
//        c->CheckTree(evt_tree,"HiTree");
//        setupEvtTree(evt_tree,offSel,1);
//    }
    Evts * offSel = &(c->evt); 

    //skim tree
//    TTree *skim_tree = (TTree*) c->inf->Get("skimanalysis/HltTree");
//    Skims my_skim;
//    if (skim_tree) {
//        c->CheckTree(skim_tree,"HltTree");
//        setupSkimTree(skim_tree,my_skim,1);
//    }
  //  Skims * my_skim = &(c->skim); 

    //hlt tree
//    TTree *hlt_tree = (TTree*) c->inf->Get("hltanalysis/HltTree");
//     Hlts trigSel;
//    if (hlt_tree) {
//        c->CheckTree(hlt_tree,"HltTree");
//        setupHltTree(hlt_tree,trigSel,1);
//    }
  //  Hlts * trigSel = &(c->hlt); 
    //jet tree
//    TTree *inp_tree = (TTree*) c->inf->Get("akPu3PFJetAnalyzer/t");
//    Jets my_ct;
//    if (inp_tree) {
//        c->CheckTree(inp_tree,"t");
//        SetupJetTree(inp_tree,my_ct,1);
//    }
    
    //    if(coll=="HI") 
   // Jets * my_ct = &(c->akPu3PF); 
    //   else 
       Jets * my_ct = &(c->ak3PF);
    
    //track tree
//    TTree *tr_tree = (TTree*) c->inf->Get("mergedTrack/trackTree");
////    TTree *tr_tree = (TTree*) c->inf->Get("anaTrack/trackTree");
//    Tracks my_tr;
//    if (tr_tree) {
//        c->CheckTree(tr_tree,"trackTree");
//        SetupTrackTree(tr_tree,my_tr,1);
//    }
    Tracks * my_tr = &(c->track);
    
    //GenParticle tree
//    TTree *GenPart_tree = (TTree*) c->inf->Get("HiGenParticleAna/hi");
//    GenParticles my_GenPart;
//    if (GenPart_tree) {
//        c->CheckTree(GenPart_tree,"hi");
//        SetupGenParticleTree(GenPart_tree,my_GenPart,1);
//    }
    GenParticles * my_GenPart = &(c->genparticle);
    cout <<"GenPart" <<my_GenPart->npart<<endl ;

    int curr_bin = nbin-1 ;
    cout <<"Number of events ="<<c->GetEntries()<<endl ;
    for(int evi = 0; evi < c->GetEntries(); evi++) {
        //          for(int evi = 0; evi < 10; evi++) {
        c->GetEntry(evi);
       //cout <<"evt = "<<evi <<endl ; 
      //  int noise_evt = my_skim->pHBHENoiseFilter ;
        //        int ecal_noise = my_skim->phiEcalRecHitSpikeFilter ;
        //        if(ecal_noise==0) continue ;
        
        double vz = offSel->vz ;
        int hiBin = offSel->hiBin ;
        
        //   cout <<"vz =" <<vz <<endl ;
        if(TMath::Abs(vz)>15.) continue ;

        //if there is no jets or no PF candidates, skip the event
        if(my_ct->nref==0) continue ;
        //put the higher pthat cut
        if(my_ct->pthat>ptmax) continue ;
        
        if(coll=="HI"){
            double centrality = hiBin*2.5 ;
            //   my_hists->CenBin->Fill(offSel->hiBin);
            
            for(int ibin = 0 ; ibin <nbin; ibin++){
                if(centrality >=centr[ibin] && centrality<centr[ibin+1]) curr_bin = ibin ;
            }
        }
        else {
            curr_bin=nbin-1 ;
        }
       //    cout << "  cent_bin:" <<curr_bin <<endl ;
        if(evi%10000==1)cout <<" coll = " <<coll <<" evt = " <<evi <<endl ;
        
        my_hists->Ntrack[curr_bin]->Fill(1);
        //cout <<my_hists->NEvents[curr_bin]->GetName() <<endl;
        //cout << "start working222222\n";
        
        
        ///--------------------------------------------------------------------
        /// ---- JETS  ----- ///////////////////////////////////////////
        ///--------------------------------------------------------------------
        
        //define the variables used
        double jet_pt = -999.;
        double jet_phi = -999.;
        double jet_eta = -999. ;
        double bkg_phi = -999.;
        double bkg_eta = -999.;
        TVector3 jet_vec;
        TVector3 bkg_vec;
        TVector3 track_vec;
        jet_vec.SetPtEtaPhi(0, 0, 0);
        bkg_vec.SetPtEtaPhi(0, 0, 0);
        track_vec.SetPtEtaPhi(0, 0, 0);
        
        int njets = 0 ;
        int ptBin = -1 ;

        int dEtaBin = -1 ;
        int nJets = 0 ;
        int flavor = 0 ;
        if(JetChoose=="Ref")
            nJets = my_ct->nref ;
        else 
            nJets = my_ct->ngen ;
        for(int j4i = 0; j4i < nJets ; j4i++) {
            if(JetChoose=="Ref"){
                jet_pt= my_ct->refpt[j4i];
                jet_phi =my_ct->refphi[j4i];
                jet_eta = my_ct->refeta[j4i];  
                flavor = TMath::Abs(my_ct->refparton_flavor[j4i]);
              //  flavor = (my_ct->refparton_flavor[j4i]);
            }
            else {
                jet_pt= my_ct->genpt[j4i];
                jet_phi =my_ct->genphi[j4i];
                jet_eta = my_ct->geneta[j4i];
            }
            if(TMath::Abs(jet_eta)>etacut) continue ;
            if(TMath::Abs(jet_eta)<etalimit) continue ;            
            if(jet_pt<leadingjetcut) continue ;
            if(jet_pt>300.) continue ;
            //remove flucluation with too large pt from small MC pthat sample
          //  if(jet_pt>=5*pthat) continue ; 
            if(my_ct->trackMax[j4i]/jet_pt <=0.01) continue ;            
            //remove too high energy jets flucutulation
            if(jet_pt>pt[nptbin]) continue ;
            if(jet_phi<=-TMath::Pi())jet_phi+=TMath::TwoPi();

            njets++ ;        
            bkg_phi = jet_phi;
            bkg_eta = -jet_eta;
            if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetBkgCone[curr_bin]->Fill(bkg_phi-jet_phi, bkg_eta-jet_eta);
            for(Int_t ieta = 0 ; ieta <netabin; ieta++){
                //   if(TMath::Abs(bkg_eta-jet_eta)>deta[ieta]&&TMath::Abs(bkg_eta-jet_eta)<=deta[ieta+1]) dEtaBin = ieta ; 
                if(TMath::Abs(jet_eta)>deta[ieta]&&TMath::Abs(jet_eta)<=deta[ieta+1]) dEtaBin = ieta ; 
            }

                // if(jet_pt<leadingjetcut) continue ;   
                if(TMath::Abs(jet_eta)>etalimit) {
                  my_hists->jetpt[curr_bin]->Fill(jet_pt);
                  if(flavor>=1&&flavor<=6) my_hists->qjetpt[curr_bin]->Fill(jet_pt);
                  if(flavor==21) my_hists->gjetpt[curr_bin]->Fill(jet_pt);
                  }
                my_hists->jeteta[curr_bin]->Fill(jet_pt, jet_eta);
                my_hists->jetphi[curr_bin]->Fill(jet_pt, jet_phi);  
                 for(Int_t ipt = 0 ; ipt <nptbin; ipt++){
                    if(jet_pt>=pt[ipt] && jet_pt <pt[ipt+1]) ptBin = ipt ;
                }
                if(TMath::Abs(jet_eta)>etalimit)my_hists->NEvents[curr_bin]->Fill(1);
                
               // cout <<"flavor =" <<flavor <<endl ;
                //seperate quark and gluon jet analysis
              //  if(flavor<1 || flavor >6) continue ; //only quark jets selected
		//   if(flavor!=21) continue ;  //only gluon jets selected
                //for track analysis
                int charge[nrbin] ;
                double sumchpt[nrbin];
                double rho[nrbin] ;
                double psi[nrbin] ;
                double sumpt = 0. ;
                double meanphi = 0.;
                double meaneta = 0. ;
                int bkgcharge[nrbin] ;
                double bkgsumchpt[nrbin]  ;
                double bkgrho[nrbin];
                double bkgpsi[nrbin] ;
                double bkgsumpt = 0. ;
                double sumTrkpt[nrbin][ntrkptbin];
                
                double jetsumTrk[ntrkptbin];
                
                double bkgsumTrkpt[nrbin][ntrkptbin];
                
                double bkgjetsumTrk[ntrkptbin];
                
                for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                    jetsumTrk[ipt]=0.;
                    bkgjetsumTrk[ipt]=0.;
                }
                for(int ir = 0 ; ir <nrbin ; ir++){
                    charge[ir]=0 ;
                    sumchpt[ir]=0. ;
                    rho[ir]=0. ;
                    psi[ir]= 0.;
                    bkgcharge[ir]=0 ;
                    bkgsumchpt[ir]=0. ;
                    bkgrho[ir]=0.;
                    bkgpsi[ir]=0. ;    
                    for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                        sumTrkpt[ir][ipt]=0.;
                        bkgsumTrkpt[ir][ipt]=0.;
                    }
                    
                    
                }
            
            jet_vec.SetPtEtaPhi(jet_pt, jet_eta, jet_phi);
            bkg_vec.SetPtEtaPhi(jet_pt, bkg_eta, bkg_phi);

            Int_t trackcount = 0 ;
            Int_t Ntrk = 0 ;
            if(TrkChoose=="Sim") Ntrk = my_tr->nParticle ;
            else  Ntrk = my_GenPart->mult ;
            //cout <<"Ntrk =" <<Ntrk <<endl ;
                for(int itr = 0 ; itr < Ntrk ; itr++){
                    double dr = 0.;
                    double bkg_dr =0.;
                    double tr_pt = 0. ;
                    double tr_phi = 0. ;
                    double tr_eta = 0. ;
                    if(TrkChoose == "Sim"){
                        tr_pt = my_tr->pPt[itr];
                        tr_phi = my_tr->pPhi[itr];
                        tr_eta = my_tr->pEta[itr];
                      //  if(my_tr->pStatus[itr]!=1) continue ;
                    }
                    else {
                        tr_pt = my_GenPart->pt[itr];
                        tr_phi = my_GenPart->phi[itr];
                        tr_eta = my_GenPart->eta[itr]; 
                        int chg = my_GenPart->chg[itr];
                        if(chg==0) continue ;        
                    }
                    if(TMath::Abs(tr_eta)>2.4) continue ;
                    if(tr_pt<trackcut) continue ;
                    
                    track_vec.SetPtEtaPhi(tr_pt, tr_eta, tr_phi);
                    dr = jet_vec.DeltaR(track_vec); 
                    bkg_dr = bkg_vec.DeltaR(track_vec); 
                    
//                        dr =TMath::Sqrt((tr_phi-jet_phi)*(tr_phi-jet_phi)+(tr_eta-jet_eta)*(tr_eta-jet_eta));
//                        bkg_dr = TMath::Sqrt((tr_phi-bkg_phi)*(tr_phi-bkg_phi)+(tr_eta-bkg_eta)*(tr_eta-bkg_eta));
                        
                        if(TMath::Abs(jet_eta)>=etalimit)my_hists->deltaR->Fill(dr); 
                        double jtrdphi = tr_phi-jet_phi ;

                        if(TMath::Abs(jet_eta)>=etalimit)my_hists->IncTrackDphi[curr_bin]->Fill(jet_pt,TMath::Abs(jtrdphi));
                        if(dr<=conesize){  //for leading jet shape study                        
                            trackcount++;
                            if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetTrackDphi[curr_bin]->Fill(jet_pt,TMath::Abs(jtrdphi));                
                            if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetTrackDeta[curr_bin]->Fill(jet_pt,tr_eta-jet_eta);                
                            if(TMath::Abs(jet_eta)>=etalimit)my_hists->jetfrag[curr_bin]->Fill(jet_pt, TMath::Log(jet_pt/tr_pt));
                            
                            if(TMath::Abs(jet_eta)>=etalimit)my_hists->trackpt[curr_bin]->Fill(jet_pt, tr_pt);
                            if(TMath::Abs(jet_eta)>=etalimit)my_hists->trackphi[curr_bin]->Fill(jet_pt, tr_phi);
                            if(TMath::Abs(jet_eta)>=etalimit)my_hists->tracketa[curr_bin]->Fill(jet_pt, tr_eta);
                            if(TMath::Abs(jet_eta)>=etalimit)my_hists->trackdr[curr_bin][ptBin]->Fill(tr_pt,dr);
                            sumpt+=tr_pt ;
                            meaneta+=tr_pt*tr_eta;
                            meanphi+=tr_pt*tr_phi;  
                            for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                              //  if(tr_pt>trkpt[ipt] && tr_pt<=trkpt[ipt+1])jetsumTrk[ipt]+=tr_pt ;
                                if(tr_pt<=trkpt[ipt+1])jetsumTrk[ipt]+=tr_pt ;
                            }
                            
                            for(int ir = 0 ; ir <nrbin; ir++){
                                if(dr>rbin[ir]&&dr<=rbin[ir+1]){
                                    rho[ir]+=tr_pt ;
                                    //     cout<<" !!!! track index =" <<itr << " trk pt =" <<tr_pt << " tr phi =" <<tr_phi << " tr eta =" <<tr_eta  <<" dr =" <<dr <<" rho =" <<rho[ir]<<endl ;
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->ChargePt[curr_bin][ir]->Fill(jet_pt, tr_pt);
                                    sumchpt[ir]+=tr_pt ;
                                    charge[ir]++ ;
                                    for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                      //  if(tr_pt>trkpt[ipt] && tr_pt<=trkpt[ipt+1])sumTrkpt[ir][ipt]+=tr_pt ;
                                        if(tr_pt<=trkpt[ipt+1])sumTrkpt[ir][ipt]+=tr_pt ;
                                    }
                                    
                                }
                                if(dr<=rbin[ir+1]) {
                                    psi[ir]+=tr_pt ; 
                                }
                            }  //radius loop for rho calculation
                        } //inside leading jet cone
                        else { //outside leading and subleading jet cone
                            if(bkg_dr<=conesize){
                                if(TMath::Abs(jet_eta)>=etalimit)my_hists->bkgtrackpt[curr_bin]->Fill(jet_pt, tr_pt);
                                if(TMath::Abs(jet_eta)>=etalimit)my_hists->bkgtrackdr[curr_bin][ptBin]->Fill(tr_pt,dr);
                                bkgsumpt+=tr_pt ;
                                for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                  //  if(tr_pt>trkpt[ipt] && tr_pt<=trkpt[ipt+1])bkgjetsumTrk[ipt]+=tr_pt ;
                                    if(tr_pt<=trkpt[ipt+1])bkgjetsumTrk[ipt]+=tr_pt ;
                                }
                                if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetBkgTrackDphi[curr_bin]->Fill(jet_pt,TMath::Abs(jtrdphi));
                                if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetBkgTrackDeta[curr_bin]->Fill(jet_pt,tr_eta-jet_eta);
                                if(TMath::Abs(jet_eta)>=etalimit)my_hists->jetbkgfrag[curr_bin]->Fill(jet_pt, TMath::Log(jet_pt/tr_pt));                            for(int ir = 0 ; ir <nrbin; ir++){
                                    if(bkg_dr>rbin[ir]&&bkg_dr<=rbin[ir+1]){
                                        if(TMath::Abs(jet_eta)>=etalimit)my_hists->bkgChargePt[curr_bin][ir]->Fill(jet_pt, tr_pt);
                                        bkgrho[ir]+=tr_pt ;
                                        bkgsumchpt[ir]+=tr_pt ;
                                        bkgcharge[ir]++ ;  
                                        for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                          //  if(tr_pt>trkpt[ipt] && tr_pt<=trkpt[ipt+1])bkgsumTrkpt[ir][ipt]+=tr_pt ;
                                            if(tr_pt<=trkpt[ipt+1])bkgsumTrkpt[ir][ipt]+=tr_pt ;
                                        }
                                        
                                    }
                                    if(bkg_dr<=rbin[ir+1]) {
                                        bkgpsi[ir]+=tr_pt ;                                       
                                    } //fill bkg subleading JS
                                }  //radius loop for leading jet backgound 
                            }  //bkg cone for leading jet bkg
                        }
                    } //track loop
                if(TMath::Abs(jet_eta)>=etalimit){ //remove overlap region for ER bkg
                    my_hists->SumptJetPt[curr_bin]->Fill(jet_pt,sumpt/jet_pt);
                    if(DoSumPtNorm){
                        if(sumpt==0.) continue ;
                    }
                    for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                        my_hists->jetPtSumTrk[curr_bin][ipt]->Fill(jet_pt, jetsumTrk[ipt]);
                        my_hists->bkgjetPtSumTrk[curr_bin][ipt]->Fill(jet_pt, bkgjetsumTrk[ipt]);
                        
                    }
                    //fill jet shape histogram 
                    for(int ir = 0 ; ir <nrbin; ir++){
                        if(DoSumPtNorm){
                            psi[ir]/=sumpt;
                            rho[ir]/=sumpt;
                            bkgpsi[ir]/=sumpt;
                            bkgrho[ir]/=sumpt;
                            sumchpt[ir]/=sumpt ;
                            bkgsumchpt[ir]/=sumpt ;
                            for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                sumTrkpt[ir][ipt]/=sumpt ;
                                bkgsumTrkpt[ir][ipt]/=sumpt ;
                            }
                            
                        }
                        else {
                            psi[ir]/=jet_pt;
                            rho[ir]/=jet_pt;
                            bkgpsi[ir]/=jet_pt;
                            bkgrho[ir]/=jet_pt;
                            sumchpt[ir]/=jet_pt ;
                            bkgsumchpt[ir]/=jet_pt ;
                            for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                sumTrkpt[ir][ipt]/=jet_pt ;
                                bkgsumTrkpt[ir][ipt]/=jet_pt ;
                            }
                            
                            
                        }
                        for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                            my_hists->jetdrSumTrkPt[curr_bin][ptBin][ipt]->Fill((rbin[ir]+deltacone/2.), sumTrkpt[ir][ipt]);
                            my_hists->bkgjetdrSumTrkPt[curr_bin][ptBin][ipt]->Fill((rbin[ir]+deltacone/2.), bkgsumTrkpt[ir][ipt]);
                        }
                        
                        my_hists->ChargeMult[curr_bin][ir]->Fill(jet_pt, charge[ir]);
                        my_hists->SumChPt[curr_bin][ir]->Fill(jet_pt, sumchpt[ir]); 
                        my_hists->bkgChargeMult[curr_bin][ir]->Fill(jet_pt, bkgcharge[ir]);
                        my_hists->bkgSumChPt[curr_bin][ir]->Fill(jet_pt, bkgsumchpt[ir]); 
                        my_hists->DiffJS[curr_bin][ir]->Fill(jet_pt, rho[ir]);
                        my_hists->IntJS[curr_bin][ir]->Fill(jet_pt, psi[ir]);
                        // rho/=deltacone ;
                        my_hists->bkgDiffJS[curr_bin][ir]->Fill(jet_pt, bkgrho[ir]);
                        //   my_hists->IntJS[curr_bin]->Fill(jet_pt, 1-psi);
                        my_hists->bkgIntJS[curr_bin][ir]->Fill(jet_pt, bkgpsi[ir]);

                        if(jet_pt>=100 &&jet_pt<120) {            
                            //   my_hists->DiffJSPt80_100[curr_bin][ir]->Fill(rho[ir]);
                            my_hists->DiffJSPt80_100[curr_bin]->Fill(rbin[ir]+0.025, rho[ir]);
                            my_hists->bkgDiffJSPt80_100[curr_bin]->Fill(rbin[ir]+0.025, bkgrho[ir]);
                            my_hists->DiffJSPt80_100Rbin[curr_bin][ir]->Fill(rho[ir]);
                            my_hists->bkgDiffJSPt80_100Rbin[curr_bin][ir]->Fill(bkgrho[ir]);
                        }
                        if(dEtaBin!=-1){        
                            my_hists->DiffEtaBinJS[curr_bin][dEtaBin][ir]->Fill(jet_pt, rho[ir]);
                            my_hists->bkgDiffEtaBinJS[curr_bin][dEtaBin][ir]->Fill(jet_pt, bkgrho[ir]);
                        }
                    } //radius bin loop            
                    my_hists->NJets[curr_bin]->Fill(1);   
                    if(bkgsumpt>sumpt)my_hists->RareEvt->Fill(offSel->hiBin);
                }  //remove overlape eta regions from jet cone and ER bkg cone

        
        } //jet loop
        //  cout << "still working222222\n";
        my_hists->Vertex->Fill(offSel->vz);
        my_hists->CenBin->Fill(offSel->hiBin); 

    }  ///event loop
    
    my_hists->Write();
    //   my_hists->Delete();
    //  delete my_hists;
    std::cout << "working done\n";
}




