#include <Riostream.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGrid.h>
#include <TFileMerger.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TF1.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>

using namespace std;

const int nrbin = 6 ;
double rbin[nrbin+1]={0.0,0.05,0.10, 0.15, 0.20,0.25,0.30};
const int nbin = 4 ;
const int centr[] ={0,10,30,50,100}; 
const double deltacone = 0.05 ;

void  LoadParameters();
//! Smearing uncertainty
TF1 *fsmear[nbin], *ftrkeff[nbin], *fpptrkeff;
Double_t smear[nbin][2]={{9.99573e-01,1.23108e-02},{9.96382e-01,6.84402e-02},{9.92243e-01,1.43513e-01},{9.99105e-01,1.64548e-02}};
Double_t trkeff[nbin][2]={{9.78076e-01,3.93014e-01},{9.80138e-01,3.87776e-01},{9.93786e-01,1.20619e-01},{9.97038e-01,8.12540e-02}};
Double_t pptrkeff[2]={9.97869e-01,3.12048e-02};

bool DoSmoothF = kFALSE ;
enum Combine_t {kRAA, kRcp, kSpect};

void  CombineSys(Combine_t com = kRAA){
    double rad[nrbin];    
    for(int i = 0 ; i <nrbin; i++){
        rad[i]=deltacone/2.+i*deltacone ;
    }
    TString smooth;
    if(DoSmoothF)smooth="SmoothLinear";
    else smooth="";
//    //QM priminary results systematics 
//    //smearing factor scaling up and down for PbPb/pp ratio
//    Double_t smear[nbin][nrbin] = {{0.117009, 0.160462, 0.355374, 0.329601, 0.631783, 0.26382},
//        {0.19069, 0.346523, 0.28749, 1.20045, 0.254487, 0.209863},
//        {0.0457928, 0.189465, 0.184973, 0.257686, 0.874293, 0.307497},
//        {0.0551462, 0.182227, 0.205381, 0.464287, 0.404984, 0.150108}
//};
//    //JEC scaling uncertainty for PbPb
//    Double_t jetscale[nbin][nrbin] = {{0.980605,2.0665,3.1524,4.2383,5.3242,6.41009},
//        {0.941142,1.61833,2.29553,2.97272,3.64991,4.3271},
//        {0.706215,1.14375,1.58129,2.01882,2.45636,2.8939},
//        {0.725743,1.23578,1.74582,2.25586,2.76589,3.27593}
//    };
//    
//    //tracking efficiency non-closure for PbPb, from myself, diffenertial JS
//    Double_t trkEffHI[nbin][nrbin] = {{1.20628, 2.68967, 4.17306, 5.65645, 7.13984, 8.62322},
//        {1.10088, 2.4446025, 3.788314, 5.1320275, 6.47574, 7.8194525},
//        {0.4065, 0.952143, 1.497787, 2.04343, 2.589073, 3.134716},
//        {0.675835, 1.527328, 2.378821, 3.230314, 4.12504, 4.9333}
//    };
//    
//    Double_t pptrkEff[nrbin] = {0.0255594, 0.59507752, 1.16459564, 1.73411376, 0.675226, 2.30363188};
//    //different bkg model, QM time
//        Double_t Bkg[nbin][nrbin] = {{0.7015095, 1.9864061, 3.2713027,4.5561993,5.8410959,7.1259925},
//            {0.7015095, 1.9864061, 3.2713027,4.5561993,5.8410959,7.1259925},
//            {0.7015095, 1.9864061, 3.2713027,4.5561993,5.8410959,7.1259925},
//            {0.7015095, 1.9864061, 3.2713027,4.5561993,5.8410959,7.1259925}
//        };
//    
//    //different correction table
//    Double_t diffCorr[nbin][nrbin] = {{0.6203, 0.220166, 2.36552, 3.14131, 3.79164, 0.140103},
//        {0.660071, 0.423419, 1.98414, 2.69331, 3.90305, 2.59888},
//        {0.0507467, 0.313536, 0.00823412, 0.247674, 0.403449, 2.58265},
//        {0.1717, 0.508969, 0.0954958, 0.311635, 0.0470193, 1.19699}
//    };
//    // end of QM priminary results systematics 
//    

//    //tracking efficiency non-closure for PbPb, from KK, weighted tracks
//    Double_t trkEffHI[nbin][nrbin] = {{2.3169, 4.04826, 5.77962, 7.51098, 9.24234, 10.9737},
//        {2.751, 3.470526, 4.190052, 4.909578, 5.629104, 6.34863},
//        {3.1635, 4.91518, 6.66686, 8.41854, 10.170220, 11.9219},
//        {3.7909, 4.338102, 4.885304, 5.432506, 5.979708, 6.52691}
//    };
    //tracking efficiency non-closure double ratio for PbPb/pp
    Double_t trkCorr[nbin][nrbin] ;
//    = {{4.82406, 5.1858, 5.5476, 5.9094, 6.2712, 6.63297},
//        {4.42025, 3.8926, 3.3651, 2.8375, 2.3099, 1.78239},
//        {4.03654, 4.7547, 5.4729, 6.1910, 6.9092, 7.62735},
//        {3.45293, 3.1562, 2.8595, 2.5628, 2.2661, 1.96936}
//    };
    //tracking efficiency non-closure central/pripheral PbPb
    Double_t trkCorrRcp[nbin][nrbin] ;
//    = {{1.42016, 2.08759, 2.75502, 3.42245, 4.089879, 4.75731},
//        {1.00191, 0.8396817, 0.6774436, 0.5152054, 0.3529672, 0.190729},
//        {0.604485, 1.63792, 2.671355, 3.70479, 4.738225, 5.77166},
//        {0,0,0,0,0,0}
//    };

    Double_t ppjetalgo[nrbin] = {0.929141, 1.10686, 2.37452, 2.26061, 2.20016, 2.10552};

    //different bkg model, new with full statistics
    Double_t Bkg[nbin][nrbin] = {{0.450306,   2.4229447, 4.3955836,   6.3682223,  8.3408612,  10.3135},
        {0.40175,    2.46172, 4.52169,   6.58166,  8.64163,   10.7016},
        {0.0821024,  0.793663, 1.5052254, 2.216787, 2.92834848,  3.63991},
        {0.0412596,  0.3822465, 0.7232336,  1.06422,  1.4052076, 1.7461946}
    };
    Double_t BkgRatio[nbin][nrbin];
    if(DoSmoothF){
    LoadParameters();
        Double_t smear[nbin][nrbin] ;
      //  Double_t jetscale[nbin][nrbin];
        Double_t trkEffHI[nbin][nrbin];
        Double_t pptrkEff[nrbin];
        for(int ibin = 0 ; ibin <nbin; ibin++){ 
            for(int ir =0 ; ir <nrbin; ir++){
                smear[ibin][ir]=TMath::Abs(fsmear[ibin]->Eval(rad[ir])-1)*100;
                trkEffHI[ibin][ir]=TMath::Abs(ftrkeff[ibin]->Eval(rad[ir])-1)*100;
                cout << "ibin =" <<ibin<<"ir =" <<ir<<"smear =" <<smear[ibin][ir]<<" eff =" <<trkEffHI[ibin][ir]<<endl;

            }
        }
        for(int ir =0 ; ir <nrbin; ir++){
            pptrkEff[ir]=TMath::Abs(fpptrkeff->Eval(rad[ir])-1)*100;  
        }
        Double_t jetscale[nbin][nrbin] = {{2.31602,3.29938,5.21023,5.58148,5.988975,6.39647},
            {2.29559,3.89185,5.01919,5.160485,5.30178,5.81605},
            {1.77396,3.39562,3.57282,4.143595,4.71437,5.285145},
            {1.9815,3.63538,4.18799,4.66353,4.901295,5.13906}
        };
    }
    else {
    Double_t smear[nbin][nrbin] ={{0.0640153, 0.243896, 1.05567, 2.32071, 0.832067,3.67566},
        {0.0274903, 0.802691, 0.889982, 2.36008, 0.61518, 1.19731},
        {0.333834, 0.342016, 2.54903, 0.339365,3.8902,2.07496},
        {0.0469137, 0.266133, 0.0272689, 0.877288, 0.698775, 1.22649}};
//    Double_t pptrkEff[nrbin] = {9.5, 8.04, 7.08, 6.12, 5.16, 4.2};
//    Double_t trkEffHI[nbin][nrbin] ={{5.87003, 1.59031, 1.7146, 1.33531, 0.910102, 2.59896},
//        {6.58632, 3.74656, 0.852732, 1.41239, 1.52152, 2.76232},
//        {6.28413, 3.20206, 0.904124, 1.92622, 1.27265, 2.11952},
//        {6.88574, 4.82951, 1.3063, 2.4074, 3.93293, 5.33402}};

        Double_t pptrkEff[nrbin] = {5, 5, 5, 5, 5, 5};
//        Double_t pptrkEff[nrbin] = {9.5, 9.5, 9.5, 9.5, 9.5, 9.5};
        Double_t trkEffHI[nbin][nrbin] ={{6.88574, 6.88574, 6.88574, 6.88574, 6.88574, 6.88574},
            {6.88574, 6.88574, 6.88574, 6.88574, 6.88574, 6.88574},
            {6.88574, 6.88574, 6.88574, 6.88574, 6.88574, 6.88574},
            {6.88574, 6.88574, 6.88574, 6.88574, 6.88574, 6.88574}};
        
  Double_t jetscale[nbin][nrbin] = {{2.31602,3.29938,5.21023,5.58148,3.94677,6.39647},
            {2.29559,3.89185,5.01919,3.88394,5.30178,5.81605},
            {1.77396,3.39562,3.57282,3.24563,4.71437,2.15097},
            {1.9815,3.63538,4.18799,4.66353,3.99957,1.08497}
        };

    }

   for(int ibin = 0 ; ibin <nbin; ibin++){ 
       for(int ir =0 ; ir <nrbin; ir++){
           if(pptrkEff[ir]<=1) trkCorr[ibin][ir] = trkEffHI[ibin][ir] ;
               else 
                   trkCorr[ibin][ir] = TMath::Abs((trkEffHI[ibin][ir]/pptrkEff[ir]));
           
           if(trkEffHI[nbin-1][ir]<=1) trkCorrRcp[ibin][ir] = trkEffHI[ibin][ir] ;
           else 
               trkCorrRcp[ibin][ir] = TMath::Abs(trkEffHI[ibin][ir]/trkEffHI[nbin-1][ir]);
           BkgRatio[ibin][ir] = TMath::Abs(Bkg[ibin][ir]/Bkg[nbin-1][ir]-1);           
           cout << "ibin =" <<ibin<<"ir =" <<ir<<"Bkg =" <<BkgRatio[ibin][ir]<<endl;
       }
   }
    
    //average the centrality dependent values if no centrality dep. found
    Double_t AverSys[nrbin]={0,0,0,0,0,0};
    for(int ir =0 ; ir <nrbin; ir++){
            for(int ibin = 0 ; ibin <nbin; ibin++){ 
                AverSys[ir]+=trkCorr[ibin][ir];
       //         cout << "ibin =" <<ibin<<"Bkg =" <<BkgRatio[ibin][ir]<<endl;
            }
        AverSys[ir]/=nbin ;
        cout <<"ir =" <<ir <<"aver =" <<AverSys[ir]<<endl ;
    }
    //using linear fitting to the data points to remove fluctuations
    TF1 * fLinear = new TF1("flinear","[0]+[1]*x", 0, 0.3);
    Double_t SmoothSys[nbin][nrbin]={{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}};
    for(int ibin = 0 ; ibin <nbin; ibin++){  
        Double_t a = (trkEffHI[ibin][nrbin-1]-trkEffHI[ibin][0])/5. ;
        Double_t b = trkEffHI[ibin][0]-a ;        
        fLinear->SetParameters(b, a);
        for(int ir =0 ; ir <nrbin; ir++){
            SmoothSys[ibin][ir]=fLinear->Eval(rad[ir]);
            cout <<"smooth =" <<SmoothSys[ibin][ir]<<endl ;
          //  trkEffHI[ibin][ir]=SmoothSys[ibin][ir] ;
            
        }
    } //finish smooth the systematics
    
    
    Double_t total[nbin][nrbin];
       for(int ibin = 0 ; ibin <nbin; ibin++){  
            for(int ir =0 ; ir <nrbin; ir++){
               if(com==kRAA)total[ibin][ir]=TMath::Sqrt(TMath::Power(smear[ibin][ir],2)+TMath::Power(trkCorr[ibin][ir],2)+TMath::Power(Bkg[ibin][ir],2)+TMath::Power(ppjetalgo[ir],2));
                else if(com==kSpect)total[ibin][ir]=TMath::Sqrt(TMath::Power(jetscale[ibin][ir],2)+TMath::Power(trkEffHI[ibin][ir],2)+TMath::Power(Bkg[ibin][ir],2));
               else 
                  // total[ibin][ir]=TMath::Sqrt(TMath::Power(trkCorrRcp[ibin][ir],2));
                 //  total[ibin][ir]=TMath::Sqrt(TMath::Power(trkCorrRcp[ibin][ir],2)+TMath::Power(BkgRatio[ibin][ir],2));
                total[ibin][ir]=TMath::Sqrt(TMath::Power(trkCorrRcp[ibin][ir],2)+TMath::Power(Bkg[ibin][ir],2));
               cout <<"ibin =" <<ibin <<"ir ="<<ir <<"sys= " <<total[ibin][ir]<<endl ;
           }
       }
    
    ofstream myfile;
    if(com==kRAA)myfile.open(Form("%ssysFixedRaaTotal.txt",smooth.Data()));
    else if(com==kSpect) myfile.open(Form("%ssysFixedProfile.txt",smooth.Data()));
    else   myfile.open(Form("%ssysRcpTotal.txt",smooth.Data()));
    for(int ibin = 0 ; ibin <nbin; ibin++){
        myfile << total[ibin][0]<<", "<< total[ibin][1]<<", "<< total[ibin][2]<<", "<< total[ibin][3]<<", "<< total[ibin][4]<<", "<< total[ibin][5]<<"\n";
    }
    myfile.close();

}

void LoadParameters()
{
    fpptrkeff = new TF1(Form("fpptrkeff"),"[0] + [1]*x",0,0.3);
    fpptrkeff->SetParameters(pptrkeff[0],pptrkeff[1]);
        for(int i=0;i<nbin;i++){
            fsmear[i]  = new TF1(Form("fsmear_%d",i),"[0] + [1]*x",0,0.3);
            ftrkeff[i] = new TF1(Form("ftrkeff_%d",i),"[0] + [1]*x",0,0.3);
            for(int im=0;im<2;im++){
                fsmear[i]->SetParameter(im,smear[i][im]);
                ftrkeff[i]->SetParameter(im,trkeff[i][im]);
            }
        }

    
    //! Reweight factor for 0-10% and 50-100%
    fReWe = new TF1("fReWe","[0] + [1]/x + [2]/x/x",30,400);
}
