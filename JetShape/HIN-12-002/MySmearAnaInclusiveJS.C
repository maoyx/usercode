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

const int pthat =50 ; //30 and 80 for pp; 30,50,80,120,170,200 for PbPb MC 
TString coll = "PP";
TString para ="Full" ; //PYQUEN parameters setting: Wide or Full
const bool DoGenAna=kFALSE ; //should be only be true for MC sample

//do tracking corrections
const bool doTrackCorrections = kTRUE; //for tracking efficiency correction
TString corrMet = "Hist" ; //Hist table from Frank, or Para from Krisztian
TrackingParam *trackCorrFromParam;

vector<TrackingCorrections*> trackCorrections;

//for jet pt smearing and reweight
const bool DoSmear = kTRUE ;
const bool DoJetPtWeight = kTRUE ; //!only be true when DoSmear TRUE 

const bool DoSumPtNorm = kFALSE ;
//weight factor from vertex and centrality
double weight =1. ;

TString intputFile ;

//TString dataPath="/Users/ymao/group/CMS/hiForest";
TString dataPath ;

//if it is pp, no centrality bins, only one
//const int nbin = 1 ;
//const int centr[] ={0,100};
////for HI centrality bin
const int nbin = 4 ;
const int centr[] ={0,10,30,50,100};
//const int nbin = 6 ;
//const int centr[] ={0,5, 10, 30,50, 70, 90};

const int ntrkptbin = 6 ;
//const double pt[]={100., 120., 140., 160., 200., 300., 500.};
const int nptbin = 1 ;
const double pt[]={100., 500.};
const double trkpt[]={1., 2., 4., 8., 16., 32., 500.};

//const Double_t jetPtBin[] = {30, 50, 70, 90, 120, 160, 200, 240, 280, 340, 400, 500};
const Double_t jetPtBin[] = {30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300,310,320,330,340,350,400};
//const Double_t jetPtBin[] = {100, 110, 120, 130, 140, 150, 160, 180, 200, 240, 300, 500};
const int nJetPtBin = sizeof(jetPtBin)/sizeof(Double_t)-1 ;

const int netabin = 4 ;
const double deta[]={0.0,0.5,1.0,1.5,2.0} ;

// ! for akPu3PF jet in pp
//const double ppmeanshift[11] ={0.0453817, 0.018024, 0.0144651, 0.0120364, 0.00819245, 0.00722936, 0.00682016, 0.00933281, 0.00360392, 0.00949609, 0.0123716};
// ! for ak3PF jet in pp
const float ppmeanshift[33] ={-0.0039463, -0.000293245, -0.00401251, -0.0025241, -0.0021653, -0.00387832, -0.00134047, -0.00417586, -0.000891294, -0.00501452, -0.00120124, -0.0052245, -0.000536276, -0.000288925, -0.000840595, -0.00920611, 0.00173915, 0.000697007, -0.00404328, 0.00269357, 0.00571458, 0.00285188, 0.0093249, 0.00309667, -0.00478382, -4.78632e-05, -0.00323687, -0.000163158, -0.00345717, -0.00303542, -0.0165555, 0.00465018, 0.00618863} ;
/*//! for akPu3PF jet JES meanshift
const double jesms[4][11]={{0.0334522, 0.0301365, 0.0270331, 0.0234527, 0.0189977, 0.0147029, 0.0112578, 0.00866234, 0.00661296, 0.00590621, 0.00793783},
                                     {0.0481923, 0.0399863, 0.0326162, 0.0245789, 0.0155205, 0.00830237, 0.00442739, 0.00389558, 0.00793215, 0.0196713, 0.0470246},
                                     {0.0457836, 0.0381061, 0.0312445, 0.0238145, 0.0155538, 0.00917213, 0.00605338, 0.0061976, 0.0109664, 0.0234188, 0.0514425},
                                     {0.0396879, 0.0321911, 0.0254954, 0.0182522, 0.0102147, 0.00403293, 0.00105551, 0.00128242, 0.00607214, 0.0184287, 0.0461193}};
// ! for akPu3PF jet JER smearing
const double jersm[4][11]={{6.88406, 8.16419, 9.08794, 9.89106, 10.509, 10.6005, 10.0495, 8.73549, 5.17097, 0, 0},
                                     {4.82393, 6.14146, 7.04037, 7.79262, 8.33827, 8.35724, 7.73099, 6.26905, 0, 0, 0},
                                     {3.00269, 4.01536, 4.66889, 5.18458, 5.49874, 5.35759, 4.63093, 2.90884, 0, 0, 0},
                                     {0.911816, 1.58212, 1.71039, 1.30248, 0, 0, 0, 0, 0, 0, 0}}; 
*/
/* //using real pp_mean value
//! for ak3PF jet JES meanshift
const float jesms[4][33]={{0.0350974, 0.0349506, 0.0374512, 0.0335472, 0.0305995, 0.0298662, 0.0251113, 0.0259674, 0.0209231, 0.023481, 0.0182716, 0.021045, 0.0152332, 0.0139716, 0.013604, 0.0211329, 0.00942357, 0.00976531, 0.0138615, 10.0713, 0.00295951, 0.00531208, -0.00163544, 0.00415035, 0.0116174, 0.00649416, 0.00931976, 0.00590439, 0.00887658, 0.00815119, 0.0213844, -9.29126e-05, -0.00236523},
                           {0.0673801, 0.0484925, 0.0431067, 0.0355834, 0.0309388, 0.0294542, 0.0244407, 0.0253035, 0.0204108, 0.023198, 0.0182574, 0.0213168, 0.015795, 0.0148197, 0.0147302, 0.0225267, 0.0110731, 0.0116583, 0.0159857, 8.98523, 0.00551096, 0.00806049, 0.00129962, 0.00726236, 0.0148972, 0.00993334, 0.0129103, 0.00963879, 0.0127479, 0.0121529, 0.0255103, 0.00415157, 0.00220406},
                           {0.068861, 0.0466669, 0.040367, 0.0327531, 0.0283239, 0.0271651, 0.0225073, 0.0237218, 0.0191625, 0.02226, 0.0176052, 0.020927, 0.0156458, 0.0148912, 0.0150044, 0.0229875, 0.0117059, 0.0124501, 0.0169248, 6.59249, 0.00671391, 0.00938203, 0.00273197, 0.00879846, 0.0165306, 0.0116582, 0.0147212, 0.0115308, 0.0147166, 0.014194, 0.0276201, 0.00632647, 0.00455551},
                           {0.0648125, 0.0409796, 0.034174, 0.0264548, 0.0220757, 0.0210278, 0.0165021, 0.0178523, 0.0134245, 0.0166459, 0.0121065, 0.015535, 0.0103522, 0.00968821, 0.00988498, 0.0179451, 0.0067348, 0.00754502, 0.0120809, 0, 0.00198009, 0.00469776, -0.00190597, 0.00420394, 0.0119769, 0.00714277, 0.0102419, 0.00708562, 0.0103035, 0.00981147, 0.0232664, 0.00200017, 0.000303577}
                        };
// ! for ak3PF jet JER smearing
const float jersm[4][33]={{7.91141, 8.26654, 8.58311, 8.86524, 9.11615, 9.33833, 9.53381, 9.7042, 9.8508, 9.97465, 10.0766, 10.1573, 10.2173, 10.2569, 10.2763, 10.2757, 10.2551, 10.2144, 10.1532, 10.0713, 9.96807, 9.84293, 9.69498, 9.52318, 9.3262, 9.1024, 8.84976, 8.56572, 8.24705, 7.88953, 7.48761, 7.03368, 5.2183},
                           {6.38946, 6.59582, 6.79117, 6.97645, 7.15244, 7.3198, 7.47913, 7.63091, 7.7756, 7.91358, 8.0452, 8.17077, 8.29055, 8.4048, 8.51375, 8.61758, 8.71648, 8.81063, 8.90017, 8.98523, 9.06595, 9.14244, 9.21481, 9.28315, 9.34754, 9.40808, 9.46483, 9.51786, 9.56724, 9.61302, 9.65525, 9.69398, 9.78956}, 
                           {5.34573, 5.26427, 5.19917, 5.15105, 5.12038, 5.10748, 5.11248, 5.13533, 5.1758, 5.23347, 5.30778, 5.39806, 5.5035, 5.62326, 5.75645, 5.90215, 6.05947, 6.22752, 6.40546, 6.59249, 6.78785, 6.99085, 7.20084, 7.41724, 7.63948, 7.86708, 8.0996, 8.3366, 8.57773, 8.82265, 9.07104, 9.32263, 10.0942},
                           {4.60364, 4.24156, 3.86577, 3.47183, 3.0527, 2.59622, 2.07793, 1.43221, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
                       };
*/
/*//! for ak3PF jet JES meanshift, set pp JES ==1 all the time
const float jesms[4][33]={{0.0311511, 0.0346574, 0.0334387, 0.0310231, 0.0284342, 0.0259879, 0.0237709, 0.0217915, 0.0200318, 0.0184665, 0.0170704, 0.0158205, 0.0146969, 0.0136827, 0.0127634, 0.0119268, 0.0111627, 0.0104623, 0.00981817, 10.0713, 0.00867409, 0.00816396, 0.00768946, 0.00724702, 0.00683353, 0.00644629, 0.0060829, 0.00574123, 0.00541941, 0.00511578, 0.00482884, 0.00455727, 0.00382339},
                           {0.0634338, 0.0481992, 0.0390942, 0.0330593, 0.0287735, 0.0255758, 0.0231002, 0.0211277, 0.0195195, 0.0181835, 0.0170562, 0.0160923, 0.0152587, 0.0145308, 0.0138896, 0.0133206, 0.0128122, 0.0123553, 0.0119424, 8.98523, 0.0112255, 0.0109124, 0.0106245, 0.010359, 0.0101134, 0.00988548, 0.00967342, 0.00947563, 0.00929072, 0.00911746, 0.00895478, 0.00880176, 0.00839268},
                           {0.0649147, 0.0463736, 0.0363545, 0.030229, 0.0261586, 0.0232867, 0.0211668, 0.019546, 0.0182713, 0.0172454, 0.016404, 0.0157025, 0.0151095, 0.0146023, 0.0141638, 0.0137814, 0.013445, 0.0131471, 0.0128815, 6.59249, 0.0124285, 0.0122339, 0.0120569, 0.0118951, 0.0117468, 0.0116103, 0.0114843, 0.0113677, 0.0112594, 0.0111586, 0.0110646, 0.0109767, 0.0107441},
                           {0.0608662, 0.0406863, 0.0301615, 0.0239307, 0.0199104, 0.0171494, 0.0151616, 0.0136764, 0.0125332, 0.0116314, 0.0109053, 0.0103105, 0.00981588, 0.00939928, 0.00904439, 0.00873903, 0.00847396, 0.00824203, 0.00803764, 0, 0.00769467, 0.00754964, 0.00741893, 0.00730061, 0.00719305, 0.0070949, 0.00700503, 0.00692246, 0.00684637, 0.00677605, 0.00671089, 0.00665035, 0.0064922}
                        };
// ! for ak3PF jet JER smearing
const float jersm[4][33]={{7.91141, 8.26654, 8.58311, 8.86524, 9.11615, 9.33833, 9.53381, 9.7042, 9.8508, 9.97465, 10.0766, 10.1573, 10.2173, 10.2569, 10.2763, 10.2757, 10.2551, 10.2144, 10.1532, 10.0713, 9.96807, 9.84293, 9.69498, 9.52318, 9.3262, 9.1024, 8.84976, 8.56572, 8.24705, 7.88953, 7.48761, 7.03368, 5.2183},
                           {6.38946, 6.59582, 6.79117, 6.97645, 7.15244, 7.3198, 7.47913, 7.63091, 7.7756, 7.91358, 8.0452, 8.17077, 8.29055, 8.4048, 8.51375, 8.61758, 8.71648, 8.81063, 8.90017, 8.98523, 9.06595, 9.14244, 9.21481, 9.28315, 9.34754, 9.40808, 9.46483, 9.51786, 9.56724, 9.61302, 9.65525, 9.69398, 9.78956},
                           {5.34573, 5.26427, 5.19917, 5.15105, 5.12038, 5.10748, 5.11248, 5.13533, 5.1758, 5.23347, 5.30778, 5.39806, 5.5035, 5.62326, 5.75645, 5.90215, 6.05947, 6.22752, 6.40546, 6.59249, 6.78785, 6.99085, 7.20084, 7.41724, 7.63948, 7.86708, 8.0996, 8.3366, 8.57773, 8.82265, 9.07104, 9.32263, 10.0942},
                           {4.60364, 4.24156, 3.86577, 3.47183, 3.0527, 2.59622, 2.07793, 1.43221, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
                       };
*/
//! for ak3PF jet JES meanshift, set pp JES ==1 all the time, using Pawan fitting
const float jesms[4][33]={{-0.0190023, 0.0103056, 0.0209182, 0.0246031, 0.0254011, 0.024924, 0.0238817, 0.0226078, 0.0212652, 0.0199344, 0.0186546, 0.0174433, 0.016307, 0.0152457, 0.0142569, 0.0133362, 0.0124789, 0.0116799, 0.0109345, 10.4618, 0.00958689, 0.00897664, 0.00840402, 0.00786589, 0.0073594, 0.00688199, 0.00643133, 0.00600532, 0.00560207, 0.00521987, 0.00485715, 0.0045125, 0.0035747},
                           {0.02251, 0.0311241, 0.0323047, 0.0311624, 0.0293057, 0.0273041, 0.0253677, 0.02357, 0.0219296, 0.020443, 0.0190985, 0.0178816, 0.0167782, 0.0157751, 0.0148604, 0.0140239, 0.0132565, 0.0125505, 0.0118991, 8.39175, 0.0107373, 0.0102173, 0.00973275, 0.00928005, 0.00885626, 0.00845875, 0.00808519, 0.0077335, 0.00740185, 0.00708859, 0.00679225, 0.0065115, 0.00575154},
                           {0.019068, 0.0287287, 0.0304181, 0.0295565, 0.0278684, 0.0259744, 0.0241099, 0.022362, 0.0207571, 0.0192964, 0.017971, 0.0167685, 0.0156761, 0.0146813, 0.0137731, 0.0129416, 0.0121782, 0.0114752, 0.0108261, 3.3134, 0.00966744, 0.00914855, 0.00866468, 0.00821249, 0.00778903, 0.0073917, 0.00701819, 0.00666646, 0.00633468, 0.00602123, 0.00572464, 0.0054436, 0.00468258},
                           {0.0232598, 0.0274357, 0.0271138, 0.0255336, 0.0236688, 0.0218399, 0.0201493, 0.0186212, 0.0172514, 0.0160254, 0.0149269, 0.0139399, 0.0130501, 0.012245, 0.0115137, 0.0108471, 0.0102374, 0.00967779, 0.00916256, 0, 0.00824614, 0.00783702, 0.0074562, 0.00710088, 0.00676862, 0.00645729, 0.00616498, 0.00589002, 0.00563093, 0.00538639, 0.00515521, 0.00493633, 0.00434451}
                        };
// ! for ak3PF jet JER smearing, using Pawan fitting
const float jersm[4][33]={{7.61535, 8.0656, 8.46349, 8.81611, 9.12871, 9.40529, 9.64893, 9.86209, 10.0467, 10.2043, 10.3362, 10.4432, 10.5262, 10.5858, 10.6223, 10.636, 10.6269, 10.595, 10.5401, 10.4618, 10.3596, 10.2328, 10.0805, 9.90141, 9.69413, 9.45679, 9.18705, 8.88198, 8.53777, 8.14947, 7.71042, 7.21137, 5.17322},
                          {5.90295, 6.24922, 6.55752, 6.83301, 7.0795, 7.29995, 7.49664, 7.67141, 7.82573, 7.96078, 8.07753, 8.17676, 8.25911, 8.32507, 8.37504, 8.4093, 8.42804, 8.43137, 8.41929, 8.39175, 8.3486, 8.28958, 8.21436, 8.12249, 8.01339, 7.88635, 7.74049, 7.57472, 7.38769, 7.17775, 6.94283, 6.68027, 5.68671},
                           {3.92274, 4.23445, 4.49054, 4.70011, 4.86917, 5.00183, 5.10092, 5.16839, 5.20546, 5.21278, 5.19047, 5.13815, 5.05489, 4.93912, 4.78849, 4.59957, 4.3674, 4.0846, 3.73971, 3.3134, 2.7683, 2.00998, 0.329073, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                           {3.65654, 3.5538, 3.43881, 3.31028, 3.16657, 3.00551, 2.82411, 2.61817, 2.38132, 2.10314, 1.7642, 1.31842, 0.548396, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
                       };
class hist_class {
public:
    hist_class(TString the_desc);
    void Delete();
    void Write();
    
    TString desc;
    bool IsMC ;
    
    TH1F * NEvents[nbin];
    TH1F * NJets[nbin];
    TH1F * Ntrack[nbin];
    TH1F * jetpt[nbin];
    TH2F * jeteta[nbin];
    TH2F * jetphi[nbin];
    TH2F * jetsmf[nbin]; //jet smear factor
    TH2F * jetms[nbin]; //mean shift
    
     TH2F * trackpt[nbin];
    TH2F * trackphi[nbin];
    TH2F * tracketa[nbin];
    TH2F * trackdr[nbin][nptbin];
    TH2F * bkgtrackdr[nbin][nptbin];
    TH2F * wttrackdr[nbin][nptbin];
     TH2F * bkgwttrackdr[nbin][nptbin];
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

    TH1F * jetdrSumTrkPt[nbin][nptbin][ntrkptbin];

    TH2F * jetPtSumTrk[nbin][ntrkptbin];

    TH1F * bkgjetdrSumTrkPt[nbin][nptbin][ntrkptbin];
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
    //   if(IsMC){
    TH2F * genjetpt[nbin];
    TH2F * genjeteta[nbin];
    TH2F * genjetphi[nbin];
    TH2F * gendphi[nbin];
     TH2F * genptRatio[nbin];
    TH2F * genphiRatio[nbin] ;
    TH2F * genetaRatio[nbin] ; 
    TH2F * genpartpt[nbin] ;
    TH2F * genpartphi[nbin] ;
    TH2F * genparteta[nbin] ;
    TH2F * bkgpartpt[nbin];
    TH2F * bkgparteta[nbin];
    TH2F * genpartdr[nbin][nptbin];
    TH2F * matchpartdr[nbin][nptbin];
    TH2F * genbkgpartdr[nbin][nptbin];
    TH2F * matchbkgpartdr[nbin][nptbin];
    
    //weighted pt distribution
    TH2F * wtgenpartdr[nbin][nptbin];
    TH2F * wtmatchpartdr[nbin][nptbin];
    TH2F * wtgenbkgpartdr[nbin][nptbin];
    TH2F * wtmatchbkgpartdr[nbin][nptbin];

    TH1F * genjetdrSumTrkPt[nbin][nptbin][ntrkptbin];
    TH1F * bkggenjetdrSumTrkPt[nbin][nptbin][ntrkptbin];

    TH2F * genSumPt[nbin] ;
    TH2F * genChargePt[nbin][nrbin];
    TH2F * genSumChPt[nbin][nrbin];
    TH2F * genChargeMult[nbin][nrbin];   
    TH2F * genDiffJS[nbin][nrbin]; //differential jet shapes (pho(r)) hist.
    TH2F * genIntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.
    //  TH2F * dphiJS[nbin]; //phi diviration from jet phi 
    //   TH2F * detaJS[nbin] ; //eta diviration from jet eta 
    
    //for jet shape background study  
    TH2F * genbkgChargePt[nbin][nrbin];
    TH2F * genbkgSumChPt[nbin][nrbin];
    TH2F * genbkgChargeMult[nbin][nrbin];   
    TH2F * genbkgDiffJS[nbin][nrbin]; //differential jet shapes (pho(r)) hist.
    TH2F * genbkgIntJS[nbin][nrbin]; //Integrated jet shapes (psi(r)) hist.

//    //for subleading jet shape, fill according leading jet bins
//    TH2F * gensubDiffJS[nbin][nptbin][nrbin]; //differential sub-jet shapes (pho(r)) hist.
//    TH2F * gensubIntJS[nbin][nptbin][nrbin]; //Integrated sub-jet shapes (psi(r)) hist.
//    TH2F * genbkgsubDiffJS[nbin][nptbin][nrbin]; //differential sub-jet shapes (pho(r)) hist.
//    TH2F * genbkgsubIntJS[nbin][nptbin][nrbin]; //Integrated sub-jet shapes (psi(r)) hist. 

    //   }
};

hist_class::hist_class(TString the_desc)
{
    
    desc = the_desc;
    IsMC =kFALSE ;
    for(int ibin = 0 ; ibin <nbin; ibin++){
        NEvents[ibin] = new TH1F(Form("Nevents_%d-%d%%",centr[ibin],centr[ibin+1]), Form("Nevents_%d-%d%%",centr[ibin],centr[ibin+1]), 100, 0, 2.);
        NJets[ibin] = new TH1F(Form("NJets_%d-%d%%",centr[ibin],centr[ibin+1]), Form("NJets_%d-%d%%",centr[ibin],centr[ibin+1]), 100, -0.5, 99.5);
        Ntrack[ibin] = new TH1F(Form("Ntracks_%d-%d%%",centr[ibin],centr[ibin+1]), Form("Ntracks_%d-%d%%",centr[ibin],centr[ibin+1]), 800, -1., 799);
        if(IsMC)jetpt[ibin] = new TH1F(Form("jetpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500.);	
        else jetpt[ibin] = new TH1F(Form("jetpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetpt_%d-%d%%",centr[ibin],centr[ibin+1]), nJetPtBin, jetPtBin);
        jetpt[ibin]->Sumw2();
      //  jeteta[ibin] = new TH2F(Form("jeteta_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jeteta_%d-%d%%",centr[ibin],centr[ibin+1]), nJetPtBin, jetPtBin,100, -5.05, 4.95);
        jeteta[ibin] = new TH2F(Form("jeteta_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jeteta_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500.,100, -5.05, 4.95);
        jeteta[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        jeteta[ibin]->GetYaxis()->SetTitle("#eta");
        jeteta[ibin]->Sumw2();
      //  jetphi[ibin] = new TH2F(Form("jetphi_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetphi_%d-%d%%",centr[ibin],centr[ibin+1]), nJetPtBin, jetPtBin,100, -5.05, 4.95);
        jetphi[ibin] = new TH2F(Form("jetphi_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetphi_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500.,100, -5.05, 4.95);
        jetphi[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        jetphi[ibin]->GetYaxis()->SetTitle("#phi");
        jetphi[ibin]->Sumw2();
      //  jetsmf[ibin] = new TH2F(Form("jetSmearFactor_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetSmearFactor_%d-%d%%",centr[ibin],centr[ibin+1]), nJetPtBin, jetPtBin,100, 0., 100);
        jetsmf[ibin] = new TH2F(Form("jetSmearFactor_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetSmearFactor_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500.,100, 0., 100);
        jetsmf[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        jetsmf[ibin]->GetYaxis()->SetTitle("Smear Factor");
        jetsmf[ibin]->Sumw2();
      //  jetms[ibin] = new TH2F(Form("jetMeanShift_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetMeanShift_%d-%d%%",centr[ibin],centr[ibin+1]), nJetPtBin, jetPtBin,100, -50., 50.);
        jetms[ibin] = new TH2F(Form("jetMeanShift_%d-%d%%",centr[ibin],centr[ibin+1]), Form("jetMeanShift_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500.,100, -50., 50.);
        jetms[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
        jetms[ibin]->GetYaxis()->SetTitle("Mean Shift");
        jetms[ibin]->Sumw2();

        trackpt[ibin] = new TH2F(Form("trackpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("trackpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, 200, 0., 200); 
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
      //  jetfrag[ibin]->GetYaxis()->SetTitle("#xi = ln(p_{T}^{jet}/p_{T}^{h})");
        jetfrag[ibin]->GetYaxis()->SetTitle("#xi = ln(p^{jet}/p_{||}^{trk})");
        jetfrag[ibin]->Sumw2();
        jetbkgfrag[ibin] = new TH2F(Form("bkgFFleading_%d-%d%%",centr[ibin],centr[ibin+1]), Form("bkgFFleading_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, 100, 0., 10.);
        jetbkgfrag[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
      //  jetbkgfrag[ibin]->GetYaxis()->SetTitle("#xi = ln(p_{T}^{jet}/p_{T}^{h})");
        jetbkgfrag[ibin]->GetYaxis()->SetTitle("#xi = ln(p^{jet}/p_{||}^{trk})");
        jetbkgfrag[ibin]->Sumw2();
        //jet pt bins loop
        for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
	  trackdr[ibin][ipt] = new TH2F(Form("recoTrackJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("recoTrackJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize*100/(deltacone*100)), 0., conesize); 
            trackdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{trk} (GeV/c)");
            trackdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");   
            trackdr[ibin][ipt]->Sumw2();
            bkgtrackdr[ibin][ipt] = new TH2F(Form("recobkgTrackJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("recobkgTrackJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize*100/(deltacone*100)), 0., conesize);
            bkgtrackdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{trk} (GeV/c)");
            bkgtrackdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");
            bkgtrackdr[ibin][ipt]->Sumw2();
            wttrackdr[ibin][ipt] = new TH2F(Form("recowtTrackJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("recowtTrackJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize*100/(deltacone*100)), 0., conesize); 
            wttrackdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{trk} (GeV/c)");
            wttrackdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");   
            wttrackdr[ibin][ipt]->Sumw2();
            bkgwttrackdr[ibin][ipt] = new TH2F(Form("recobkgwtTrackJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("recobkgwtTrackJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize*100/(deltacone*100)), 0., conesize);
            bkgwttrackdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{trk} (GeV/c)");
            bkgwttrackdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");
            bkgwttrackdr[ibin][ipt]->Sumw2();            
            for(Int_t itr = 0 ; itr <ntrkptbin ; itr++){ //track bins loop
                jetdrSumTrkPt[ibin][ipt][itr]=new TH1F(Form("JetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1],trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), Form("JetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1],trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), (Int_t)(conesize*100/(deltacone*100)), 0., conesize);
                jetdrSumTrkPt[ibin][ipt][itr]->GetXaxis()->SetTitle("radius r"); 
                jetdrSumTrkPt[ibin][ipt][itr]->GetYaxis()->SetTitle("#Sigma p_{T}^{trk}/p_{T}^{jet}"); 
                jetdrSumTrkPt[ibin][ipt][itr]->Sumw2();
   
                bkgjetdrSumTrkPt[ibin][ipt][itr]=new TH1F(Form("bkgJetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1],trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), Form("bkgJetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1],trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), (Int_t)(conesize*100/(deltacone*100)), 0., conesize);
                bkgjetdrSumTrkPt[ibin][ipt][itr]->GetXaxis()->SetTitle("radius r"); 
                bkgjetdrSumTrkPt[ibin][ipt][itr]->GetYaxis()->SetTitle("#Sigma p_{T}^{trk}/p_{T}^{jet}"); 
                bkgjetdrSumTrkPt[ibin][ipt][itr]->Sumw2();
            }
        }
        SumptJetPt[ibin] = new TH2F(Form("SumptJetPtRatio_%d-%d%%",centr[ibin],centr[ibin+1]), Form("SumptJetPtRatio_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500,100, 0., 1.);  
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
        bkgtrackpt[ibin] = new TH2F(Form("bkgtrackpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("bkgtrackpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500, 200, 0., 200); 
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
        if(IsMC&&DoGenAna){
            genjetpt[ibin] = new TH2F(Form("genjetpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genjetpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500., 500, 0., 500.);
            genjetpt[ibin]->GetXaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
            genjetpt[ibin]->GetYaxis()->SetTitle("p_{T}^{gen} (GeV/c)");   
            genjetpt[ibin]->Sumw2();
            genjeteta[ibin] = new TH2F(Form("genjeteta_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genjeteta_%d-%d%%",centr[ibin],centr[ibin+1]), 100, -5.05, 4.95,100, -5.05, 4.95);
            genjeteta[ibin]->GetXaxis()->SetTitle("#eta^{rec}");
            genjeteta[ibin]->GetYaxis()->SetTitle("#eta^{gen}");
            genjeteta[ibin]->Sumw2();
            genjetphi[ibin] = new TH2F(Form("genjetphi_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genjetphi_%d-%d%%",centr[ibin],centr[ibin+1]), 100, -5.05, 4.95,100, -5.05, 4.95);
            genjetphi[ibin]->GetXaxis()->SetTitle("#phi^{rec}");
            genjetphi[ibin]->GetYaxis()->SetTitle("#phi^{gen}");
            genjetphi[ibin]->Sumw2();
            genptRatio[ibin] = new TH2F(Form("genptRatio_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genptRatio_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0.,500., 51, -0.05,5.05);
            genptRatio[ibin]->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
            genptRatio[ibin]->GetYaxis()->SetTitle("p_{T}^{rec}/p_{T}^{gen}"); 
            genptRatio[ibin]->Sumw2();
            genphiRatio[ibin] = new TH2F(Form("genphiRatio_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genphiRatio_%d-%d%%",centr[ibin],centr[ibin+1]), 500,0., 500.,100, -1.01,0.99);
            genphiRatio[ibin]->GetXaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
            genphiRatio[ibin]->GetYaxis()->SetTitle("#phi_{gen}-#phi_{rec}"); 
            genphiRatio[ibin]->Sumw2();
            genetaRatio[ibin] = new TH2F(Form("genetaRatio_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genetaRatio_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500.,100, -1.01,0.99);
            genetaRatio[ibin]->Sumw2();
            genetaRatio[ibin]->GetXaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
            genetaRatio[ibin]->GetYaxis()->SetTitle("#eta_{gen}-#eta_{rec}"); 
            genpartpt[ibin] = new TH2F(Form("genpartpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genpartpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500., 200, 0., 200);
            genpartpt[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            genpartpt[ibin]->GetYaxis()->SetTitle("p_{T}^{part} (GeV/c)");   
            genpartpt[ibin]->Sumw2();
            genpartphi[ibin] = new TH2F(Form("genpartphi_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genpartphi_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500., 100, -5.05, 4.95);
            genpartphi[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            genpartphi[ibin]->GetYaxis()->SetTitle("#phi^{part}");   
            genpartphi[ibin]->Sumw2();
            genparteta[ibin] = new TH2F(Form("genparteta_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genparteta_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500., 100, -5.05, 4.95);
            genparteta[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            genparteta[ibin]->GetYaxis()->SetTitle("#eta^{part} ");   
            genparteta[ibin]->Sumw2();
            bkgpartpt[ibin] = new TH2F(Form("bkgpartpt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("bkgpartpt_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500., 200, 0., 200); 
            bkgpartpt[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            bkgpartpt[ibin]->GetYaxis()->SetTitle("p_{T}^{part} (GeV/c)");   
            bkgpartpt[ibin]->Sumw2();
            bkgparteta[ibin] = new TH2F(Form("bkgparteta_%d-%d%%",centr[ibin],centr[ibin+1]), Form("bkgparteta_%d-%d%%",centr[ibin],centr[ibin+1]), 200, 0., 200, 100, -5.05, 4.95); 
            bkgparteta[ibin]->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
            bkgparteta[ibin]->GetYaxis()->SetTitle("#eta^{part} (GeV/c)");   
            bkgparteta[ibin]->Sumw2();

            genSumPt[ibin] = new TH2F(Form("genSumPt_%d-%d%%",centr[ibin],centr[ibin+1]), Form("genSumPto_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500.,500, 0., 500.);  
            genSumPt[ibin]->GetXaxis()->SetTitle("#Sum p_{T}^{rec}(R) (GeV/c)");
            genSumPt[ibin]->GetYaxis()->SetTitle("#Sum p_{T}^{gen}(R) (GeV/c)");
            genSumPt[ibin]->Sumw2();
            jetaxisRes[ibin] = new TH2F(Form("JetAxisRes_%d-%d%%",centr[ibin],centr[ibin+1]), Form("JetAxisRes_%d-%d%%",centr[ibin],centr[ibin+1]), 500, 0., 500., 200, -1., 1.);
            jetaxisRes[ibin]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
            jetaxisRes[ibin]->GetYaxis()->SetTitle("#Delta r=r_{genAxis}-r_{recoAxis}");
            jetaxisRes[ibin]->Sumw2();
            for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
                genpartdr[ibin][ipt] = new TH2F(Form("genPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("genPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize*100/(deltacone*100)), 0., conesize); 
                genpartdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
                genpartdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");   
                genpartdr[ibin][ipt]->Sumw2();
                matchpartdr[ibin][ipt] = new TH2F(Form("matchPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("matchPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize*100/(deltacone*100)), 0., conesize); 
                matchpartdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
                matchpartdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");   
                matchpartdr[ibin][ipt]->Sumw2();
                genbkgpartdr[ibin][ipt] = new TH2F(Form("genbkgPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("genbkgPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize*100/(deltacone*100)), 0., conesize); 
                genbkgpartdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
                genbkgpartdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");   
                genbkgpartdr[ibin][ipt]->Sumw2();
                matchbkgpartdr[ibin][ipt] = new TH2F(Form("matchbkgPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("matchbkgPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize*100/(deltacone*100)), 0., conesize); 
                matchbkgpartdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
                matchbkgpartdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");   
                matchbkgpartdr[ibin][ipt]->Sumw2();

                wtgenpartdr[ibin][ipt] = new TH2F(Form("wtgenPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("wtgenPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize*100/(deltacone*100)), 0., conesize); 
                wtgenpartdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
                wtgenpartdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");   
                wtgenpartdr[ibin][ipt]->Sumw2();
                wtmatchpartdr[ibin][ipt] = new TH2F(Form("wtmatchPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("wtmatchPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize*100/(deltacone*100)), 0., conesize); 
                wtmatchpartdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
                wtmatchpartdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");   
                wtmatchpartdr[ibin][ipt]->Sumw2();
                wtgenbkgpartdr[ibin][ipt] = new TH2F(Form("wtgenbkgPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("wtgenbkgPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize*100/(deltacone*100)), 0., conesize); 
                wtgenbkgpartdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
                wtgenbkgpartdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");   
                wtgenbkgpartdr[ibin][ipt]->Sumw2();
                wtmatchbkgpartdr[ibin][ipt] = new TH2F(Form("wtmatchbkgPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), Form("wtmatchbkgPartJetPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1], centr[ibin],centr[ibin+1]), 200, 0., 200, (Int_t)(conesize*100/(deltacone*100)), 0., conesize); 
                wtmatchbkgpartdr[ibin][ipt]->GetXaxis()->SetTitle("p_{T}^{part} (GeV/c)");
                wtmatchbkgpartdr[ibin][ipt]->GetYaxis()->SetTitle("radius r");   
                wtmatchbkgpartdr[ibin][ipt]->Sumw2();

                for(Int_t itr = 0 ; itr <ntrkptbin ; itr++){ //track bins loop
                    genjetdrSumTrkPt[ibin][ipt][itr]=new TH1F(Form("genJetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1],trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), Form("genJetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1],trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), (Int_t)(conesize*100/(deltacone*100)), 0., conesize);
                    genjetdrSumTrkPt[ibin][ipt][itr]->GetXaxis()->SetTitle("radius r"); 
                    genjetdrSumTrkPt[ibin][ipt][itr]->GetYaxis()->SetTitle("#Sigma p_{T}^{trk}/p_{T}^{jet}"); 
                    genjetdrSumTrkPt[ibin][ipt][itr]->Sumw2();
                    
                    bkggenjetdrSumTrkPt[ibin][ipt][itr]=new TH1F(Form("bkggenJetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1],trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), Form("bkggenJetPt%.f_%.fSumTrkPt%.f_%.f_%d-%d%%",pt[ipt],pt[ipt+1],trkpt[itr],trkpt[itr+1], centr[ibin],centr[ibin+1]), (Int_t)(conesize*100/(deltacone*100)), 0., conesize);
                    bkggenjetdrSumTrkPt[ibin][ipt][itr]->GetXaxis()->SetTitle("radius r"); 
                    bkggenjetdrSumTrkPt[ibin][ipt][itr]->GetYaxis()->SetTitle("#Sigma p_{T}^{trk}/p_{T}^{jet}"); 
                    bkggenjetdrSumTrkPt[ibin][ipt][itr]->Sumw2();
                } //track bins loop

            }  //jet bins loop
        } //MC histogram
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
        if(IsMC&&DoGenAna){
            for(int ir = 0 ; ir <nrbin; ir++){
                genChargePt[ibin][ir] = new TH2F(Form("genchargeptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("genchargeptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 500, 0., 500);
                genChargePt[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                genChargePt[ibin][ir]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
                genChargePt[ibin][ir]->Sumw2();
                genSumChPt[ibin][ir] = new TH2F(Form("gensumchptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("gensumchargeptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
                genSumChPt[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                genSumChPt[ibin][ir]->GetYaxis()->SetTitle("#Sigma p_{T}^{h^{#pm}}/p_{T}^{jet}");
                genSumChPt[ibin][ir]->Sumw2();
                genChargeMult[ibin][ir] = new TH2F(Form("genChMultdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("genChMultdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 100, -0.5, 99.5);
                genChargeMult[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                genChargeMult[ibin][ir]->GetYaxis()->SetTitle("# of Charge");   
                genChargeMult[ibin][ir]->Sumw2();
                
                genDiffJS[ibin][ir] = new TH2F(Form("gendifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("gendifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
                genDiffJS[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                genDiffJS[ibin][ir]->GetYaxis()->SetTitle("#rho (r)");
                genDiffJS[ibin][ir]->Sumw2();
                genIntJS[ibin][ir] = new TH2F(Form("genIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("genIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
                genIntJS[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                genIntJS[ibin][ir]->GetYaxis()->SetTitle(Form("#psi(r=%1.2f)",deltacone*(ir+1)));
                genIntJS[ibin][ir]->Sumw2();
                //for bkg histos.
                genbkgChargePt[ibin][ir] = new TH2F(Form("genbkgchargeptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("genbkgchargeptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 500, 0., 500);
                genbkgChargePt[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                genbkgChargePt[ibin][ir]->GetYaxis()->SetTitle("p_{T} (GeV/c)");
                genbkgChargePt[ibin][ir]->Sumw2();
                genbkgSumChPt[ibin][ir] = new TH2F(Form("genbkgsumchptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("genbkgsumchargeptdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
                genbkgSumChPt[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                genbkgSumChPt[ibin][ir]->GetYaxis()->SetTitle("#Sigma p_{T}^{h^{#pm}}/p_{T}^{jet}");    
                genbkgSumChPt[ibin][ir]->Sumw2();
                genbkgChargeMult[ibin][ir] = new TH2F(Form("genbkgChMultdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("genbkgChMultdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 100, -0.5, 99.5);
                genbkgChargeMult[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                genbkgChargeMult[ibin][ir]->GetYaxis()->SetTitle("# of Charge");
                genbkgChargeMult[ibin][ir]->Sumw2();
                
                genbkgDiffJS[ibin][ir] = new TH2F(Form("genbkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("genbkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
                genbkgDiffJS[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                genbkgDiffJS[ibin][ir]->GetYaxis()->SetTitle("#rho (r)");
                genbkgDiffJS[ibin][ir]->Sumw2();
                genbkgIntJS[ibin][ir] = new TH2F(Form("genbkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("genbkgIntegratedJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
                genbkgIntJS[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                genbkgIntJS[ibin][ir]->GetYaxis()->SetTitle(Form("#psi(r=%1.2f)",deltacone*(ir+1)));
                genbkgIntJS[ibin][ir]->Sumw2();
                GenAxisDiffJS[ibin][ir] = new TH2F(Form("GenAxisdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("GenAxisdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
                GenAxisDiffJS[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                GenAxisDiffJS[ibin][ir]->GetYaxis()->SetTitle("#rho (r)");
                GenAxisDiffJS[ibin][ir]->Sumw2();
                GenAxisbkgDiffJS[ibin][ir] = new TH2F(Form("GenAxisbkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), Form("GenAxisbkgdifferentialJSdR%.f_%.f_Cen%d-%d%%",rbin[ir]*100,rbin[ir+1]*100,centr[ibin],centr[ibin+1]), 500, 0., 500, 150, -0.005, 1.495);
                GenAxisbkgDiffJS[ibin][ir]->GetXaxis()->SetTitle("p_{T}^{jet} (GeV/c)");
                GenAxisbkgDiffJS[ibin][ir]->GetYaxis()->SetTitle("#rho (r)");
                GenAxisbkgDiffJS[ibin][ir]->Sumw2();
                                
            }

        } //for MC histograms
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
        delete jeteta[ibin];
        delete jetphi[ibin];
        if(DoSmear){
            delete jetsmf[ibin];
            delete jetms[ibin];
        }      
        delete trackpt[ibin];
        delete trackphi[ibin];
        delete tracketa[ibin];
        
        delete jetfrag[ibin];
        delete jetbkgfrag[ibin];
        for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
            delete trackdr[ibin][ipt];
            delete bkgtrackdr[ibin][ipt];
            delete wttrackdr[ibin][ipt];
            delete bkgwttrackdr[ibin][ipt];
            for(Int_t itr = 0 ; itr <ntrkptbin ; itr++){
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
        if(IsMC&&DoGenAna){
            delete genjetpt[ibin];
            delete genjeteta[ibin];
            delete genjetphi[ibin];
            delete genptRatio[ibin];
            delete genphiRatio[ibin];
            delete genetaRatio[ibin];
            delete genpartpt[ibin];
            delete genpartphi[ibin];
            delete genparteta[ibin];
            delete bkgpartpt[ibin];
            delete bkgparteta[ibin];
            delete genSumPt[ibin];
            delete jetaxisRes[ibin];
            for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
                delete genpartdr[ibin][ipt];
                delete matchpartdr[ibin][ipt];
                delete genbkgpartdr[ibin][ipt];
                delete matchbkgpartdr[ibin][ipt];
                delete wtgenpartdr[ibin][ipt];
                delete wtmatchpartdr[ibin][ipt];
                delete wtgenbkgpartdr[ibin][ipt];
                delete wtmatchbkgpartdr[ibin][ipt];
                for(Int_t itr = 0 ; itr <ntrkptbin ; itr++){
                    delete genjetdrSumTrkPt[ibin][ipt][itr];
                    delete bkggenjetdrSumTrkPt[ibin][ipt][itr];
                }
            }
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
        if(IsMC&&DoGenAna){
            for(int ir = 0 ; ir <nrbin; ir++){
                delete genChargePt[ibin][ir];
                delete genSumChPt[ibin][ir];
                delete genChargeMult[ibin][ir];
                delete genDiffJS[ibin][ir];
                delete genIntJS[ibin][ir];
                delete genbkgChargePt[ibin][ir];      
                delete genbkgSumChPt[ibin][ir];
                delete genbkgChargeMult[ibin][ir];
                delete genbkgDiffJS[ibin][ir];
                delete genbkgIntJS[ibin][ir];
                delete GenAxisDiffJS[ibin][ir];
                delete GenAxisbkgDiffJS[ibin][ir];
 
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
    if(IsMC) dataType="MC" ;
    else dataType="DATA";
    TString met ;
    if(etalimit==etacut) met = "RC";
    else met = "ER";
    TString level ;
    if(DoJetPtWeight) level = "wt";
    else level = "";
    TString Norm ;
    if(DoSumPtNorm)Norm="NormSum";
    else Norm ="NormJet";
    if(doTrackCorrections){
     if(IsMC==kTRUE && coll=="PP")   
          out_name=Form("%s%s_MySmearPawanFitAk3PF%sIncJetPt%.f_%sCorrTrk%.fEtaCut%.fLimit%.f_%sbkgJS%sCone%.f_CenBin%d_Nrbin%d_LJbin%dpthat%d_%s",dataType.Data(),coll.Data(),level.Data(),leadingjetcut,corrMet.Data(), trackcut,etacut*10, etalimit*10,met.Data(), Norm.Data(),conesize*10, nbin,nrbin,nptbin,pthat, intputFile.Data());
       else 
         out_name=Form("%s%s_MySmearPawanFitAk3PF%sIncJetPt%.f_%sCorrTrk%.fEtaCut%.fLimit%.f_%sbkgJS%sCone%.f_CenBin%d_Nrbin%d_LJbin%d_%s",dataType.Data(),coll.Data(),level.Data(),leadingjetcut,corrMet.Data(), trackcut,etacut*10, etalimit*10,met.Data(), Norm.Data(),conesize*10, nbin,nrbin,nptbin, intputFile.Data());
      }
       else {
    if(IsMC==kTRUE && coll=="PP")
            out_name=Form("%s%s_AkPu3PF%sIncJetPt%.f_Trk%.fEtaCut%.fLimit%.f_%sbkgJS%sCone%.f_CenBin%d_Nrbin%d_LJbin%dpthat%d_%s",dataType.Data(),coll.Data(),level.Data(),leadingjetcut,trackcut,etacut*10, etalimit*10,met.Data(), Norm.Data(),conesize*10, nbin,nrbin,nptbin,pthat, intputFile.Data());
    else        
            out_name=Form("%s%s_AkPu3PF%sIncJetPt%.f_Trk%.fEtaCut%.fLimit%.f_%sbkgJS%sCone%.f_CenBin%d_Nrbin%d_LJbin%d_%s",dataType.Data(),coll.Data(),level.Data(),leadingjetcut,trackcut,etacut*10, etalimit*10,met.Data(), Norm.Data(),conesize*10, nbin,nrbin,nptbin, intputFile.Data());       
    }
     
       TFile *out_file = new TFile(Form("/net/hidsk0001/d00/scratch/maoyx/JetShape/CMSSW_4_4_2_patch5/src/HiForest_V2_02_16/%s",out_name.Data()), "RECREATE");
    for(int ibin = 0 ; ibin <nbin; ibin++){
        NEvents[ibin]->Write();
        NJets[ibin]->Write();
        Ntrack[ibin]->Write();
        jetpt[ibin]->Write();
        jeteta[ibin]->Write();
        jetphi[ibin]->Write();
        if(DoSmear){
                jetsmf[ibin]->Write();
                jetms[ibin]->Write();
        }
         trackpt[ibin]->Write();
        trackphi[ibin]->Write();
        tracketa[ibin]->Write();
        
        jetfrag[ibin]->Write();
        jetbkgfrag[ibin]->Write();
        for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
            trackdr[ibin][ipt]->Write();
            bkgtrackdr[ibin][ipt]->Write();
            wttrackdr[ibin][ipt]->Write();
            bkgwttrackdr[ibin][ipt]->Write();
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
        if(IsMC&&DoGenAna){
            genjetpt[ibin]->Write();
            genjeteta[ibin]->Write();
            genjetphi[ibin]->Write();
             genptRatio[ibin]->Write();
            genphiRatio[ibin]->Write();
            genetaRatio[ibin]->Write();
            genpartpt[ibin]->Write();
            genpartphi[ibin]->Write();
            genparteta[ibin]->Write();
            bkgpartpt[ibin]->Write();
            bkgparteta[ibin]->Write();
            genSumPt[ibin]->Write();
            jetaxisRes[ibin]->Write();
            for(Int_t ipt = 0 ; ipt <nptbin ; ipt++){
                genpartdr[ibin][ipt]->Write();
                matchpartdr[ibin][ipt]->Write();
                genbkgpartdr[ibin][ipt]->Write();
                matchbkgpartdr[ibin][ipt]->Write();
                wtgenpartdr[ibin][ipt]->Write();
                wtmatchpartdr[ibin][ipt]->Write();
                wtgenbkgpartdr[ibin][ipt]->Write();
                wtmatchbkgpartdr[ibin][ipt]->Write();
                for(Int_t itr = 0 ; itr <ntrkptbin ; itr++){
                    genjetdrSumTrkPt[ibin][ipt][itr]->Write();
                    bkggenjetdrSumTrkPt[ibin][ipt][itr]->Write();
                }
            }    
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
            if(IsMC&&DoGenAna){
                for(int ir = 0 ; ir <nrbin; ir++){
                genChargePt[ibin][ir]->Write();    
                genSumChPt[ibin][ir]->Write(); 
                genChargeMult[ibin][ir]->Write();    
                genDiffJS[ibin][ir]->Write();
                genIntJS[ibin][ir]->Write();
                genbkgChargePt[ibin][ir]->Write();    
                genbkgSumChPt[ibin][ir]->Write(); 
                genbkgChargeMult[ibin][ir]->Write();   
                genbkgDiffJS[ibin][ir]->Write();
                genbkgIntJS[ibin][ir]->Write();
                GenAxisDiffJS[ibin][ir]->Write();
                GenAxisbkgDiffJS[ibin][ir]->Write();

                }
            } //only write if it MC
        } //centrality bins
    
    //  for(int i=0; i<6;i++) ptbin[i]->Write();
    deltaR->Write();
    CenBin->Write();
    RareEvt->Write();
    Vertex->Write();

    out_file->Close();
    cout <<"Output file: " <<Form("%s",out_name.Data()) <<endl ;
    
}


void MySmearAnaInclusiveJS()
{
    // gROOT->ForceStyle(0);
    //for centrality reweight parameterization
    TF1 *fcen = new TF1("fcen","[0]*exp([1]+[2]*x+[3]*x*x+[4]*x*x*x)",0,40);
   // fcen->SetParameters(2.10653e-02,5.61607,-1.41493e-01,1.00586e-03,-1.32625e-04); //from Pawn/Pelin
    fcen->SetParameters(1.98261e-02,5.55963e+00,-1.34951e-01,1.70895e-03,-9.28386e-05); //fit by myself
   
  //  TF1 * fcen = new TF1("fcen","[0]+[1]*x+[2]*x*x+[3]*TMath::Power(x,3)", 0, 40);
  //  fcen->SetParameters(3.55169, -0.246077, 0.00541725,-3.64594e-05);
  //  fcen->SetParErrors(0.143397, 0.0166513, 0.000639065,7.94585e-06);
    //for vertex reweight parameterization
    TF1 * fVz = new TF1("fVx","[0]+[1]*x+[2]*TMath::Power(x,2)+[3]*TMath::Power(x,3)+[4]*TMath::Power(x,4)", -15., 15.);
    if(coll=="HI")
      //  fVz->SetParameters(0.803816, -0.0179222, 0.00716733, -0.000165785, 7.30741e-05); //PbPb vertex reweighting
        fVz->SetParameters(7.62788e-01,-1.13228e-02,5.85199e-03,-3.04550e-04,4.43440e-05); 
    else 
     //    fVz->SetParameters(0.86946, -0.0353677, 0.00497902, -0.00016535, 4.53564e-05); //p vertex reweighting
      fVz->SetParameters(8.41684e-01,-2.58609e-02,4.86550e-03,-3.10581e-04,2.07918e-05); 
 //   fVz->SetParErrors(0.0218623, 0.00619651,  0.00144498,0.000105209, 1.3927e-05); //for PbPb vertext reweight parameters error
//   fVz->SetParErrors( 0.0280385, 0.00792399,  0.00173743,0.000122401, 1.59032e-05); //for PbPb vertext reweight parameters error
        
    hist_class *my_hists = new hist_class("pfjet");
    cout <<my_hists->IsMC<<endl ;

    //! Edited by Pawan
    //! Load the Smearing factors
    if(coll=="PP"&& DoSmear)LoadParameters();
    
//for pt jet reweight after smearing
    TH1F * ptwt[nbin];
    TFile * ptwtfile ;    
    if(DoJetPtWeight && DoSmear && coll=="PP"){
     //     TFile * ptwtfile = TFile::Open("/Users/ymao/group/CMS/macros/DataJetPtReweight.root", "readonly");
  //      TFile * ptwtfile = TFile::Open("/home/group/CMS/WorkSpace/analysis/PFjet/DataJetPtReweight.root", "readonly");
        //  TFile * ptwtfile = TFile::Open("/afs/cern.ch/user/y/ymao/scratch0/JetShape/DataSubjetPtReweight.root", "readonly");     
     //   TFile * ptwtfile = TFile::Open("/net/hidsk0001/d00/scratch/maoyx/JetShape/CMSSW_4_4_2_patch5/src/DataJetPtReweight.root", "readonly");
        if(my_hists->IsMC==kTRUE)  
            ptwtfile = TFile::Open("/net/hidsk0001/d00/scratch/maoyx/JetShape/CMSSW_4_4_2_patch5/src/MCJetPtReweightRebin.root", "readonly");
        else  
         //  ptwtfile = TFile::Open("/net/hidsk0001/d00/scratch/maoyx/JetShape/CMSSW_4_4_2_patch5/src/DataJetPtReweightRebin.root", "readonly");
       //    ptwtfile = TFile::Open("/net/hidsk0001/d00/scratch/maoyx/JetShape/CMSSW_4_4_2_patch5/src/DataJetPtReweight.root", "readonly");
           ptwtfile = TFile::Open("/net/hidsk0001/d00/scratch/maoyx/JetShape/CMSSW_4_4_2_patch5/src/DataMySmearPawanFitAk3PFJetPtReweightRebin.root", "readonly");
        for(int ibin = 0 ; ibin <nbin; ibin++){
          //  ptwt[ibin] = (TH1F*)ptwtfile->Get(Form("WeightFactor_JetPtCen%d-%d%%", centr[ibin],centr[ibin+1]));
            ptwt[ibin] = (TH1F*)ptwtfile->Get(Form("SmearedRatioJetPtCen%d-%d%%", centr[ibin],centr[ibin+1]));

        }
        //   ptwtfile->Close();
    }

    std::cout << "start working\n";
    if(my_hists->IsMC!=kTRUE){ 
      //    dataPath="/Users/ymao/group/CMS/hiForest";
   //  dataPath="/home/group/CMS/WorkSpace/data/HI2011/forest" ; //lxplus data path
      if(coll=="HI")
        //  dataPath="/net/hisrv0001/home/yenjie/slocal/merge/pbpb" ; //mit PbPb data for QM analysis
      //    dataPath="/net/hisrv0001/home/yenjie/slocal/merge/pbpbDijet" ; //mit PbPb data full dataset
        dataPath="/net/hidsk0001/d00/scratch/yjlee/merge/pbpbDijet_v20" ;//mit PbPb data path
   else 
        dataPath="/net/hisrv0001/home/yenjie/scratch/hiForest/prod/productionPP/CMSSW_4_4_2_patch6/test/ppForest2";
    }
    else {
    //         dataPath= "/Users/ymao/group/CMS/hiForest"; //local analysis
        if(coll=="HI") {
        //  dataPath= Form("/net/hisrv0001/home/yenjie/slocal/merge/v28/pthat%d",pthat); //lxplus MC normial
          if(pthat==50||pthat==80||pthat==100||pthat==170)
             dataPath= Form("/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/"); //MIT MC normial
           else 
                dataPath= Form("/mnt/hadoop/cms/store/user/yenjie/HiForest_v28/"); //MIT MC normial
        //  dataPath= Form("/net/hisrv0001/home/yenjie/slocal/merge/pyquen/pyquen%s%d",para.Data(), pthat); //lxplus MC pquen
       }
        else    
         dataPath= Form("/net/hisrv0001/home/zhukova/scratch/HIHighPt/forest/pthat%d", pthat); //lxplus path for pp
       //  dataPath= Form("/mnt/hadoop/cms/store/user/icali/Pythia/Z2/ppDijet_merged"); //lxplus path for pp
    }
    if(my_hists->IsMC!=kTRUE){  //real data analysis
        if(coll=="HI")             
          //  intputFile="mergedFile.root" ;  //QM analysis results
            intputFile="promptskim-hihighpt-hltjet80-pt90-v20.root" ; //full dataset
         //  intputFile="castor_frank_PbPb_dijetforest2_HiForest-promptskim-hihighpt-hltjet80-pt90-v2_v3_part2.root";
	    //  intputFile="HiForest-promptskim-hihighpt-hltjet80-pt90-v3_part2.root";
       // intputFile="promptskim-hihighpt-hltjet80-pt90-v3-forest-v2.root";
        else 
	      intputFile="pp_merged_full.root";
     //    intputFile="forest2_2760GeV_ppJet_full_July14.root";
	  //   intputFile="castor_frank_forest2_HiForest-ppskim-hihighpt-pt90-v1_v3_part.root";
    }
    else { //MC sample
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
    }
    
    TString inname=Form("%s/%s", dataPath.Data(),intputFile.Data());
    // Define the input file and HiForest
    bool isPP =false ;
    if(coll=="PP") isPP =true ;
    HiForest * c = new HiForest(inname,"forest",isPP);
   // HiForest *c = new HiForest(inname);
    c->doTrackCorrections = 1;
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
    Skims * my_skim = &(c->skim); 

    //hlt tree
//    TTree *hlt_tree = (TTree*) c->inf->Get("hltanalysis/HltTree");
//     Hlts trigSel;
//    if (hlt_tree) {
//        c->CheckTree(hlt_tree,"HltTree");
//        setupHltTree(hlt_tree,trigSel,1);
//    }
    Hlts * trigSel = &(c->hlt); 
    //jet tree
//    TTree *inp_tree = (TTree*) c->inf->Get("akPu3PFJetAnalyzer/t");
//    Jets my_ct;
//    if (inp_tree) {
//        c->CheckTree(inp_tree,"t");
//        SetupJetTree(inp_tree,my_ct,1);
//    }
    
  //  Jets * my_ct = &(c->akPu3PF); 
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

    int curr_bin = nbin-1 ;
    cout <<"Number of events ="<<c->GetEntries()<<endl ;
    for(int evi = 0; evi < c->GetEntries(); evi++) {
        //          for(int evi = 0; evi < 10; evi++) {
        c->GetEntry(evi);
       //  if(evi%2==1) continue ; 
       //cout <<"evt = "<<evi <<endl ; 
        int noise_evt = my_skim->pHBHENoiseFilter ;
        //        int ecal_noise = my_skim->phiEcalRecHitSpikeFilter ;
        //        if(ecal_noise==0) continue ;
        
        double vz = offSel->vz ;
        int hiBin = offSel->hiBin ;
        weight = 1. ;       
        //   cout <<"vz =" <<vz <<endl ;
        if(TMath::Abs(vz)>15.) continue ;
        if(my_hists->IsMC==kTRUE) weight*=fVz->Eval(vz);
        else weight=1. ;
        
 //  if(coll=="HI") {
       if(my_hists->IsMC!=kTRUE){
            int evt_sel = my_skim->pcollisionEventSelection ;
            if(evt_sel==0) continue ;
        }
        if(my_hists->IsMC!=kTRUE){
            if(noise_evt==0) continue ;
              if(coll=="HI"){
	      int jetTr2 = trigSel->HLT_HIJet80_v1 ;
            //    else 
            //            int jetTr1 = trigSel->HLT_Jet40_v1 ;
	           	  //    int jetTr2 = trigSel->HLT_Jet60_v1 ;
            //                
            if(jetTr2==0) continue ;
	}
        }
                
        //if there is no jets or no PF candidates, skip the event
        if(my_ct->nref==0) continue ;
        
        if(coll=="HI"){
            if(my_hists->IsMC==kTRUE) 
                // weight = wt[hiBin]; 
                weight*=fcen->Eval(hiBin);
            // weight = 1.;
            else 
                weight = 1. ;
            double centrality = hiBin*2.5 ;
            //   my_hists->CenBin->Fill(offSel->hiBin);
            
            for(int ibin = 0 ; ibin <nbin; ibin++){
                if(centrality >=centr[ibin] && centrality<centr[ibin+1]) curr_bin = ibin ;
            }
        }
        else {
            curr_bin=nbin-1 ;
            //   weight = 1. ;
        }
       //   weight = 1. ;

       //    cout << "  cent_bin:" <<curr_bin <<endl ;
        if(evi%10000==1)cout <<" coll = " <<coll <<" weight = " <<weight <<" evt = " <<evi <<endl ;
        
        my_hists->Ntrack[curr_bin]->Fill(1, weight);
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
        
        int njets = 0 ;
        int ptBin = -1 ;

        int dEtaBin = -1. ;
        
        //   if(DoGenAna){
        //for MC gen jet info
        double genjet_pt = -999. ;
        double genjet_phi = -999.;
        double genjet_eta = -999. ;
        double genbkg_phi = -999. ;
        double genbkg_eta = -999. ;
    
        int gencharge[nrbin] ;
        double gensumchpt[nrbin];
        double genrho[nrbin] ;
        double genpsi[nrbin] ;
        double gensumpt = 0. ;
        int genbkgcharge[nrbin] ;
        double genbkgsumchpt[nrbin]  ;
        double genbkgrho[nrbin];
        double genbkgpsi[nrbin] ;
        double genbkgsumpt = 0. ;
        double gensumTrkpt[nrbin][ntrkptbin];
        double genbkgsumTrkpt[nrbin][ntrkptbin];
        
        //    }
        
        //jet pt smearing variables according to PbPb resolution
        double smpt[nbin];
        double smbkgjetsumTrk[nbin][ntrkptbin];
        double smjetsumTrk[nbin][ntrkptbin];
        double smsumTrkpt[nbin][nrbin][ntrkptbin];
        double smbkgsumTrkpt[nbin][nrbin][ntrkptbin];
        
        double smsumpt[nbin] ;
        double smbkgsumpt[nbin] ;
        int smcharge[nbin][nrbin] ;
        int smbkgcharge[nbin][nrbin] ;
        double smsumchpt[nbin][nrbin];
        double smrho[nbin][nrbin] ;
        double smpsi[nbin][nrbin] ;
        double smbkgsumchpt[nbin][nrbin]  ;
        double smbkgrho[nbin][nrbin];
        double smbkgpsi[nbin][nrbin] ;
        double smGenAxisrho[nbin][nrbin] ;
        double smGenAxisbkgrho[nbin][nrbin] ;
        double flavor = 0 ;
        for(int j4i = 0; j4i < my_ct->nref ; j4i++) {
            
            jet_pt= my_ct->jtpt[j4i];
            if(TMath::Abs(my_ct->jteta[j4i])>etacut) continue ;
            if(!DoSmear && jet_pt>500.) continue ;
            if(DoSmear && jet_pt<50.) continue ;
          //  if(TMath::Abs(my_ct->jteta[j4i])<etalimit) continue ;            
          //  if(!DoSmear && jet_pt<leadingjetcut) continue ;
          //  if(jet_pt<leadingjetcut) continue ;
            //remove flucluation with too large pt from small MC pthat sample
            if(my_hists->IsMC==kTRUE && jet_pt>=5*pthat) continue ; 
            if(my_ct->trackMax[j4i]/jet_pt <=0.01) continue ;            
            //remove too high energy jets flucutulation
            if(jet_pt>pt[nptbin]) continue ;
            jet_phi =my_ct->jtphi[j4i];
            if(jet_phi<=-TMath::Pi())jet_phi+=TMath::TwoPi();
            jet_eta = my_ct->jteta[j4i];
            if(my_hists->IsMC==kTRUE) flavor = TMath::Abs(my_ct->refparton_flavor[j4i]);
             //seperate quark and gluon jet analysis
          //  if(my_hists->IsMC==kTRUE) {
          //   if(flavor<1 || flavor >6) continue ; //only quark jets selected
          //  if(flavor!=21) continue ;  //only gluon jets selected
           //  }
            njets++ ;        
            bkg_phi = jet_phi;
            bkg_eta = -jet_eta;
            if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetBkgCone[curr_bin]->Fill(bkg_phi-jet_phi, bkg_eta-jet_eta, weight);
            for(Int_t ieta = 0 ; ieta <netabin; ieta++){
                //   if(TMath::Abs(bkg_eta-jet_eta)>deta[ieta]&&TMath::Abs(bkg_eta-jet_eta)<=deta[ieta+1]) dEtaBin = ieta ; 
                if(TMath::Abs(jet_eta)>deta[ieta]&&TMath::Abs(jet_eta)<=deta[ieta+1]) dEtaBin = ieta ; 
            }
            
            if(my_hists->IsMC==kTRUE){
//                //search for leading at the generator level with all matched jets
//                for(int j4i = 0; j4i < my_ct->nref ; j4i++){
//                    
//                    double refjet_pt= my_ct->refpt[j4i];
//                    if(TMath::Abs(my_ct->refeta[j4i])>etacut) continue ;
//                    if(refjet_pt<leadingjetcut) continue ;
//                    //remove flucluation with too large pt from small MC pthat sample
//                    if(my_hists->IsMC==kTRUE && refjet_pt>5*pthat) continue ; 
//                    if(refjet_pt >genjet_pt){
//                        genjet_pt=refjet_pt;
//                        genj4i=j4i;
//                    }
//                } //search for generator leading jet loop
//                if(genj4i<0) continue ; //no leading found
//                genjet_pt =my_ct->refpt[genj4i];  
//                genjet_phi =my_ct->refphi[genj4i]; 
//                if(genjet_phi<=-TMath::Pi())genjet_phi+=TMath::TwoPi();            
//                genjet_eta =my_ct->refeta[genj4i];  
                
                    genjet_pt =my_ct->refpt[j4i];  
                    genjet_phi =my_ct->refphi[j4i]; 
                 if(genjet_phi<=-TMath::Pi())genjet_phi+=TMath::TwoPi();            
                    genjet_eta =my_ct->refeta[j4i];  
               // cout <<"gen jet =" <<genjet_pt <<endl ;
                if(genjet_pt<0.) continue ; 
               genbkg_phi = genjet_phi;
               genbkg_eta = -genjet_eta;
            }            
            if(coll=="PP"&&DoSmear){
                for(Int_t ibin = 0 ; ibin <nbin ; ibin++){
                    smpt[ibin] = jet_pt ;
                    smsumpt[ibin]=0.;
                    smbkgsumpt[ibin]=0.;
                    for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                        smjetsumTrk[ibin][ipt]=0.;
                        smbkgjetsumTrk[ibin][ipt]=0.;
                    }
                    for(int ir = 0 ; ir <nrbin; ir++){
                        smcharge[ibin][ir]=0 ;
                        smbkgcharge[ibin][ir]=0 ;
                        smpsi[ibin][ir]=0.0;
                        smrho[ibin][ir]=0.0;
                        smbkgpsi[ibin][ir]=0.0;
                        smbkgrho[ibin][ir]=0.0;
                        smsumchpt[ibin][ir]=0.0;
                        smbkgsumchpt[ibin][ir]=0.0;     
                        smGenAxisrho[ibin][ir]=0. ;
                        smGenAxisbkgrho[ibin][ir]=0. ;
                        for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                            smsumTrkpt[ibin][ir][ipt]=0.;
                            smbkgsumTrkpt[ibin][ir][ipt]=0.;
                        }
                    }                    
                }
                Int_t JetBin = -1 ;
                for(Int_t ijet = 0 ; ijet < nJetPtBin ; ijet++){
                  if(jet_pt>jetPtBin[ijet] && jet_pt <=jetPtBin[ijet+1]) JetBin = ijet ; 
                 }
                for(int nic=0;nic<4;nic++){
              //  for(int nic=0;nic<6;nic++){
                    double ptReFactor = 1. ;
                    float fsmpt ;
                    float shiftpt ;
                    double fjetsmf ;
                    double jetsms ;
                    //! Smearing function
                    //      float fsmpt = GetSmearedPtData(3,nic,jet_pt);
                  //  if(my_hists->IsMC==kTRUE){
                       //  shiftpt = jet_pt + ppmeanshift[JetBin]*jet_pt-jesms[nic][JetBin]*jet_pt ;
                         shiftpt = jet_pt-jesms[nic][JetBin]*jet_pt ;
                         fsmpt = gRandom->Gaus(shiftpt,jersm[nic][JetBin]);
                         fjetsmf = jersm[nic][JetBin] ;
                         jetsms = jesms[nic][JetBin]-ppmeanshift[JetBin] ;
                      //  fsmpt = GetSmearedPtMC(2,nic,jet_pt,genjet_pt);
                      //  fsmpt = GetSmearedPtMC_NoMeanShift(2,nic,jet_pt,genjet_pt);
                      //  fjetsmf = GetSmFactor(2, nic, genjet_pt) ;
                      //  jetsms = GetMeanShift(2,nic, genjet_pt) ;
		//	}	
                /*    else {
                        fsmpt = GetSmearedPtData(2,nic,my_ct->jtpt[j4i],0, "no");
                    //  fsmpt = GetSmearedPtData(2,nic,jet_pt,5, "low");
                        fjetsmf = GetSmFactor(2, nic, my_ct->jtpt[j4i]) ;
                        jetsms = GetMeanShift(2,nic, my_ct->jtpt[j4i]) ;
			}*/
                  //  if(fsmpt<leadingjetcut) continue ;
                    //      cout <<"starting smearing"<<endl ;
                    //	   cout <<"sm factor =" <<fsmpt <<"lead _pt" <<jet_pt<<endl ;
                    //! For different centrality binning                                                                                                                             
                      //  cout <<"nic =" <<nic <<"jet pt ="<<jet_pt<<"smeared ="<<fsmpt<<endl ;
                    int ibin ; 
                  /*  if(nbin==2){
                        if(nic==0 || nic==1 || nic ==3)ibin=0; 
                      //  else if(nic==4 || nic==5)ibin=1;   
                        else ibin =1 ;
                    }
                    else if(nbin==4){
                        if(nic==0 || nic==1)ibin=0; 
                        else if(nic==2)ibin=1; 
                        else if(nic==3)ibin=2;  
                       // else if(nic==4 || nic==5)ibin=3;   
                        else ibin = 3 ;
                    }
                  else */ 
                    ibin = nic ;
                    smpt[ibin]=fsmpt ;
                    if(my_hists->IsMC==kTRUE){
                       my_hists->jetsmf[ibin]->Fill(genjet_pt, fjetsmf,weight);
                       my_hists->jetms[ibin]->Fill(genjet_pt, jetsms, weight);
                        my_hists->genjetpt[ibin]->Fill(jet_pt,genjet_pt, weight);
                        my_hists->genjeteta[ibin]->Fill(jet_eta, genjet_eta, weight);
                        my_hists->genjetphi[ibin]->Fill(jet_phi, genjet_phi, weight);
                        my_hists->genptRatio[ibin]->Fill(genjet_pt,smpt[ibin]/genjet_pt, weight);
                        my_hists->genphiRatio[ibin]->Fill(jet_pt,genjet_phi-jet_phi, weight);
                        my_hists->genetaRatio[ibin]->Fill(jet_pt,genjet_eta-jet_eta, weight);

                    }
                    else {
                        my_hists->jetsmf[ibin]->Fill(my_ct->jtpt[j4i], fjetsmf, weight);
                        my_hists->jetms[ibin]->Fill(my_ct->jtpt[j4i], jetsms, weight);
                    }
                    //cout <<"ibin =" <<ibin <<"jet pt ="<<jet_pt<<"smeared ="<<fsmpt<<endl ;
                      //  cout <<"ibin =" <<ibin << "curr_bin =" << curr_bin<<"jet pt ="<<jet_pt<<"smeared ="<<smpt[ibin]<<endl ;
                     if(smpt[ibin]<leadingjetcut) continue ;
                     if(smpt[ibin]>500.) continue ;
                    if(DoJetPtWeight){
                        if(ptwt[ibin]->GetBinContent(ptwt[ibin]->FindBin(fsmpt)))
                            ptReFactor=ptwt[ibin]->GetBinContent(ptwt[ibin]->FindBin(fsmpt)); 
                        else ptReFactor = 1. ;
                        //     cout <<"ptwt =" <<ptReFactor <<endl ;
                    }
                    // jet_pt*=ptReFactor ;
 
                    my_hists->jetpt[ibin]->Fill(smpt[ibin], weight*ptReFactor);
                    my_hists->jeteta[ibin]->Fill(smpt[ibin], jet_eta, weight*ptReFactor);
                    my_hists->jetphi[ibin]->Fill(smpt[ibin], jet_phi, weight*ptReFactor);
                    for(Int_t ipt = 0 ; ipt <nptbin; ipt++){
                        if(fsmpt>=pt[ipt] && fsmpt <pt[ipt+1]) ptBin = ipt ;
                    }
                    if(TMath::Abs(my_ct->jteta[j4i])>etalimit) { //continue ;
                        my_hists->NEvents[ibin]->Fill(1, weight); 
                        
                        //for tracking correction factors        
                        float trkweight = 1.; 
                        Int_t trackcount = 0 ;
                        for(int itr = 0 ; itr < my_tr->nTrk ; itr++){
                            double dr = 0.;
                            double bkg_dr =0.;
                            double tr_pt = my_tr->trkPt[itr];
                            double tr_phi = my_tr->trkPhi[itr];
                            double tr_eta = my_tr->trkEta[itr];
                            if(TMath::Abs(tr_eta)>2.4) continue ;
                            
                            //  if(my_tr->trkNHit[itr]<8 ||(my_tr->highPurity[itr]&&my_tr->trkAlgo[itr]==4.0)){
                            if(my_tr->trkAlgo[itr]<4 ||(my_tr->highPurity[itr])){
                                //            if(my_tr->highPurity[itr]==0) continue ;
                                //            if(my_tr->trkAlgo[itr]!=4.0) continue ;
                                // if(TMath::Abs(tr_eta)>tracketacut) continue ;
                                if(tr_pt<trackcut) continue ;
                                
                                trkweight=1. ;
                                if(doTrackCorrections){
                                    if(corrMet=="Hist")trkweight = c->getTrackCorrection(itr);  
                                    else trkweight = c->getTrackCorrectionPara(itr); 
                                    //  if(TMath::Abs(trkweight-1.)<0.001)trkweight = c->getTrackCorrection(itr);  
                                }
                                // cout <<"ibin =" <<ibin <<"trkwt ="<<trkweight<<endl ;
                                dr =TMath::Sqrt((tr_phi-jet_phi)*(tr_phi-jet_phi)+(tr_eta-jet_eta)*(tr_eta-jet_eta));
                                bkg_dr = TMath::Sqrt((tr_phi-bkg_phi)*(tr_phi-bkg_phi)+(tr_eta-bkg_eta)*(tr_eta-bkg_eta));
                                double jtrdphi = tr_phi-jet_phi ;
                                //    cout <<"itr =" <<itr <<"jet pt ="<<jet_pt<<"smeared ="<<smpt[ibin]<<endl ;
                                if(dr<=conesize){  //for leading jet shape study
                                    
                                    trackcount++;
                                    smsumpt[ibin]+=tr_pt*trkweight ;
                                    for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                        // if(tr_pt*trkweight>trkpt[ipt] && tr_pt*trkweight<=trkpt[ipt+1])smjetsumTrk[ibin][ipt]+=tr_pt*trkweight ;
                                        if(tr_pt*trkweight<=trkpt[ipt+1])smjetsumTrk[ibin][ipt]+=tr_pt*trkweight ;
                                    }
                                    //   double ptReFactor = 1. ;
                                    if(DoJetPtWeight){
                                        if(ptwt[ibin]->GetBinContent(ptwt[ibin]->FindBin(smpt[ibin])))
                                            ptReFactor=ptwt[ibin]->GetBinContent(ptwt[ibin]->FindBin(smpt[ibin])); 
                                        else ptReFactor = 1. ;
                                    }
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetTrackDphi[ibin]->Fill(smpt[ibin],TMath::Abs(jtrdphi), weight*trkweight*ptReFactor);                
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetTrackDeta[ibin]->Fill(smpt[ibin],tr_eta-jet_eta, weight*trkweight*ptReFactor);                
                                    
                                    //  if(TMath::Abs(jet_eta)>=etalimit)my_hists->jetfrag[ibin]->Fill(smpt[ibin], TMath::Log(smpt[ibin]/tr_pt), weight*trkweight*ptReFactor);
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->jetfrag[ibin]->Fill(smpt[ibin], TMath::Log((smpt[ibin]*TMath::CosH(jet_eta))/(tr_pt*TMath::CosH(tr_eta)*TMath::Cos(dr))), weight*trkweight*ptReFactor);
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->trackpt[ibin]->Fill(smpt[ibin], tr_pt, weight*trkweight*ptReFactor);
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->trackphi[ibin]->Fill(smpt[ibin], tr_phi, weight*trkweight*ptReFactor);
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->tracketa[ibin]->Fill(smpt[ibin], tr_eta, weight*trkweight*ptReFactor);
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->trackdr[ibin][ptBin]->Fill(tr_pt,dr, weight*trkweight);
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->wttrackdr[ibin][ptBin]->Fill(tr_pt,dr, tr_pt*weight*trkweight);
                                    
                                    for(int ir = 0 ; ir <nrbin; ir++){
                                        if(dr>rbin[ir]&&dr<=rbin[ir+1]){
                                            if(TMath::Abs(jet_eta)>=etalimit)my_hists->ChargePt[ibin][ir]->Fill(smpt[ibin], tr_pt, weight*trkweight*ptReFactor);
                                            smrho[ibin][ir]+=tr_pt*trkweight ;
                                            smsumchpt[ibin][ir]+=tr_pt*trkweight ;
                                            smcharge[ibin][ir]++ ;
                                            for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                                //   if(tr_pt*trkweight>trkpt[ipt] && tr_pt*trkweight<=trkpt[ipt+1])smsumTrkpt[ibin][ir][ipt]+=tr_pt*trkweight ;
                                                if(tr_pt*trkweight<=trkpt[ipt+1])smsumTrkpt[ibin][ir][ipt]+=tr_pt*trkweight ;
                                            }
                                        }
                                        if(dr<=rbin[ir+1]) {
                                            smpsi[ibin][ir]+=tr_pt*trkweight ; 
                                        }
                                        //  cout <<" ibin ="<<ibin<<"smpt =" <<smpt[ibin]<<" ir =" <<ir<<"rho =" <<smrho[ibin][ir]<<endl ;
                                    }  //radius loop to fill leading jet JS histogram                             
                                } //inside leading jet cone
                                else { //outside leading and subleading jet cone
                                    if(bkg_dr<=conesize){
                                        if(TMath::Abs(jet_eta)>=etalimit)my_hists->bkgtrackdr[ibin][ptBin]->Fill(tr_pt,bkg_dr, weight*trkweight);
                                        if(TMath::Abs(jet_eta)>=etalimit)my_hists->bkgwttrackdr[ibin][ptBin]->Fill(tr_pt,bkg_dr, tr_pt*weight*trkweight);
                                        if(TMath::Abs(jet_eta)>=etalimit)my_hists->bkgtrackpt[ibin]->Fill(smpt[ibin], tr_pt, weight*trkweight*ptReFactor);
                                        for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                            // if(tr_pt*trkweight>trkpt[ipt] && tr_pt*trkweight<=trkpt[ipt+1])smbkgjetsumTrk[ibin][ipt]+=tr_pt*trkweight ;
                                            if(tr_pt*trkweight<=trkpt[ipt+1])smbkgjetsumTrk[ibin][ipt]+=tr_pt*trkweight ;
                                        }
                                        smbkgsumpt[ibin]+=tr_pt*trkweight ;
                                        if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetBkgTrackDphi[ibin]->Fill(smpt[ibin],TMath::Abs(jtrdphi), weight*trkweight*ptReFactor);
                                        if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetBkgTrackDeta[ibin]->Fill(smpt[ibin],tr_eta-jet_eta, weight*trkweight*ptReFactor);
                                        
                                        //   if(TMath::Abs(jet_eta)>=etalimit)my_hists->jetbkgfrag[ibin]->Fill(smpt[ibin], TMath::Log(smpt[ibin]/tr_pt), weight*trkweight);
                                        if(TMath::Abs(jet_eta)>=etalimit)my_hists->jetbkgfrag[ibin]->Fill(smpt[ibin], TMath::Log((smpt[ibin]*TMath::CosH(jet_eta))/(tr_pt*TMath::CosH(tr_eta)*TMath::Cos(bkg_dr))), weight*trkweight*ptReFactor);
                                        for(int ir = 0 ; ir <nrbin; ir++){
                                            if(bkg_dr>rbin[ir]&&bkg_dr<=rbin[ir+1]){
                                                if(TMath::Abs(jet_eta)>=etalimit) my_hists->bkgChargePt[ibin][ir]->Fill(smpt[ibin], tr_pt, weight*trkweight);
                                                smbkgrho[ibin][ir]+=tr_pt*trkweight ;
                                                smbkgsumchpt[ibin][ir]+=tr_pt*trkweight ;
                                                smbkgcharge[ibin][ir]++ ;  
                                                for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                                    //  if(tr_pt*trkweight>trkpt[ipt] && tr_pt*trkweight<=trkpt[ipt+1])smbkgsumTrkpt[ibin][ir][ipt]+=tr_pt*trkweight ;
                                                    if(tr_pt*trkweight<=trkpt[ipt+1])smbkgsumTrkpt[ibin][ir][ipt]+=tr_pt*trkweight ;
                                                }
                                                
                                            }
                                            if(bkg_dr<=rbin[ir+1]) {
                                                smbkgpsi[ibin][ir]+=tr_pt*trkweight ;                                       
                                            } //fill bkg subleading JS
                                        }
                                    }  //bkg cone for leading jet bkg
                                } //outside leading and subleading jet
                                if(my_hists->IsMC==kTRUE&&DoGenAna){
                                    double GenAxisdr =TMath::Sqrt((tr_phi-genjet_phi)*(tr_phi-genjet_phi)+(tr_eta-genjet_eta)*(tr_eta-genjet_eta));
                                    double GenAxisbkgdr =TMath::Sqrt((tr_phi-genbkg_phi)*(tr_phi-genbkg_phi)+(tr_eta-genbkg_eta)*(tr_eta-genbkg_eta));
                                    
                                    if(dr<=conesize&&TMath::Abs(jet_eta)>=etalimit)my_hists->jetaxisRes[ibin]->Fill(smpt[ibin], (GenAxisdr-dr), weight);
                                    if(GenAxisdr<=conesize){
                                        for(int ir = 0 ; ir <nrbin; ir++){
                                            if(GenAxisdr>rbin[ir]&&GenAxisdr<=rbin[ir+1]){
                                                smGenAxisrho[ibin][ir]+=tr_pt*trkweight ;
                                            }
                                        }  //radius loop to fill leading jet JS histogram
                                    } //tracks inside for leading jet
                                    if(GenAxisbkgdr<=conesize){
                                        //   if(doTrackCorrections)trkweight = c->getTrackCorrection(itr);  
                                        if(doTrackCorrections){
                                            if(corrMet=="Hist")trkweight = c->getTrackCorrection(itr);  
                                            else trkweight = c->getTrackCorrectionPara(itr); 
                                        }
                                        for(int ir = 0 ; ir <nrbin; ir++){
                                            if(GenAxisbkgdr>rbin[ir]&&GenAxisbkgdr<=rbin[ir+1]){
                                                smGenAxisbkgrho[ibin][ir]+=tr_pt*trkweight ;
                                            }
                                        }  //radius loop to fill leading jet bkg JS histogram
                                    } //leading jet with gen axis background JS
                                } //fill only when MC 
                            } //track selection  
                        } //track loop  
                        
                        //     if(TMath::Abs(jet_eta)>=etalimit){ //remove overlap region for ER bkg
                        my_hists->SumptJetPt[ibin]->Fill(smpt[ibin],smsumpt[ibin]/smpt[ibin], weight);
                        if(DoSumPtNorm && smsumpt[ibin]==0.) continue ;
                        //   double ptReFactor = 1. ;
                        //                double smjetpt = smpt[ibin];
                        //                double smsubjetpt = smsubpt[ibin]; 
                        if(DoJetPtWeight){
                            if(ptwt[ibin]->GetBinContent(ptwt[ibin]->FindBin(smpt[ibin])))
                                ptReFactor=ptwt[ibin]->GetBinContent(ptwt[ibin]->FindBin(smpt[ibin])); 
                            else ptReFactor = 1. ;
                        }
                        for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                            my_hists->jetPtSumTrk[ibin][ipt]->Fill(smpt[ibin], smjetsumTrk[ibin][ipt],weight*ptReFactor);
                            my_hists->bkgjetPtSumTrk[ibin][ipt]->Fill(smpt[ibin], smbkgjetsumTrk[ibin][ipt],weight*ptReFactor);
                            
                        }
                        
                        for(int ir = 0 ; ir <nrbin; ir++){
                            //  cout <<" ibin ="<<ibin<<"smpt =" <<smpt[ibin]<<" ir =" <<ir<<"rho =" <<smrho[ibin][ir]<<endl ;
                            if(DoSumPtNorm){
                                smpsi[ibin][ir]/=smsumpt[ibin];
                                smrho[ibin][ir]/=smsumpt[ibin];
                                smbkgpsi[ibin][ir]/=smsumpt[ibin];
                                smbkgrho[ibin][ir]/=smsumpt[ibin];
                                smsumchpt[ibin][ir]/=smsumpt[ibin];
                                smbkgsumchpt[ibin][ir]/=smsumpt[ibin];   
                                for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                    smsumTrkpt[ibin][ir][ipt]/=smsumpt[ibin] ;
                                    smbkgsumTrkpt[ibin][ir][ipt]/=smsumpt[ibin] ;
                                }
                                
                                //    cout <<"dem =" <<smpt[ibin]<<" ir =" <<ir<<"rho =" <<smrho[ibin][ir]<<endl ;
                            }
                            else {
                                smpsi[ibin][ir]/=smpt[ibin];
                                smrho[ibin][ir]/=smpt[ibin];
                                smbkgpsi[ibin][ir]/=smpt[ibin];
                                smbkgrho[ibin][ir]/=smpt[ibin];
                                smsumchpt[ibin][ir]/=smpt[ibin] ;
                                smbkgsumchpt[ibin][ir]/=smpt[ibin] ; 
                                for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                    smsumTrkpt[ibin][ir][ipt]/=smpt[ibin] ;
                                    smbkgsumTrkpt[ibin][ir][ipt]/=smpt[ibin] ;
                                }
                                
                                //   cout <<"dem =" <<smpt[ibin]<<" ir =" <<ir<<"rho =" <<smrho[ibin][ir]<<endl ;
                            }
                            for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                my_hists->jetdrSumTrkPt[ibin][ptBin][ipt]->Fill((rbin[ir]+0.025), smsumTrkpt[ibin][ir][ipt]*weight);
                                my_hists->bkgjetdrSumTrkPt[ibin][ptBin][ipt]->Fill((rbin[ir]+0.025), smbkgsumTrkpt[ibin][ir][ipt]*weight);
                            }
                            
                            //  cout <<"ijet =" <<j4i <<"jet pt ="<<smpt[ibin]<<endl ;
                            my_hists->ChargeMult[ibin][ir]->Fill(smpt[ibin], smcharge[ibin][ir], weight*ptReFactor);
                            my_hists->SumChPt[ibin][ir]->Fill(smpt[ibin], smsumchpt[ibin][ir],weight*ptReFactor); 
                            my_hists->bkgChargeMult[ibin][ir]->Fill(smpt[ibin], smbkgcharge[ibin][ir],weight*ptReFactor);
                            my_hists->bkgSumChPt[ibin][ir]->Fill(smpt[ibin], smbkgsumchpt[ibin][ir],weight*ptReFactor); 
                            my_hists->DiffJS[ibin][ir]->Fill(smpt[ibin], smrho[ibin][ir],weight*ptReFactor);
                            //  my_hists->DiffJS[ibin][ir]->Fill(smpt[ibin], smrho[ibin][ir],weight);
                            my_hists->IntJS[ibin][ir]->Fill(smpt[ibin], smpsi[ibin][ir],weight*ptReFactor);
                            // rho/=deltacone ;
                            my_hists->bkgDiffJS[ibin][ir]->Fill(smpt[ibin], smbkgrho[ibin][ir],weight*ptReFactor);
                            //  my_hists->bkgDiffJS[ibin][ir]->Fill(smpt[ibin], smbkgrho[ibin][ir],weight);
                            //   my_hists->IntJS[curr_bin]->Fill(jet_pt, 1-psi);
                            my_hists->bkgIntJS[ibin][ir]->Fill(smpt[ibin], smbkgpsi[ibin][ir],weight*ptReFactor);
                            if(my_hists->IsMC==kTRUE&&DoGenAna){ 
                                if(DoSumPtNorm){
                                    smGenAxisrho[ibin][ir]/=smsumpt[ibin];
                                    smGenAxisbkgrho[ibin][ir]/=smsumpt[ibin];
                                }
                                else {
                                    smGenAxisrho[ibin][ir]/=smpt[ibin];
                                    smGenAxisbkgrho[ibin][ir]/=smpt[ibin];
                                }
                                my_hists->GenAxisDiffJS[ibin][ir]->Fill(smpt[ibin], smGenAxisrho[ibin][ir],weight);
                                my_hists->GenAxisbkgDiffJS[ibin][ir]->Fill(smpt[ibin], smGenAxisbkgrho[ibin][ir],weight);
                            }
                            if(smpt[ibin]*ptReFactor>=100 &&smpt[ibin]*ptReFactor<120) {            
                                //   my_hists->DiffJSPt80_100[curr_bin][ir]->Fill(rho[ir]);
                                my_hists->DiffJSPt80_100[ibin]->Fill(rbin[ir]+0.025, smrho[ibin][ir]*weight);
                                my_hists->bkgDiffJSPt80_100[ibin]->Fill(rbin[ir]+0.025, smbkgrho[ibin][ir]*weight);
                                my_hists->DiffJSPt80_100Rbin[ibin][ir]->Fill(smrho[ibin][ir],weight);
                                my_hists->bkgDiffJSPt80_100Rbin[ibin][ir]->Fill(smbkgrho[ibin][ir],weight);
                            }
                            
                            if(dEtaBin!=-1){
                                my_hists->DiffEtaBinJS[ibin][dEtaBin][ir]->Fill(smpt[ibin], smrho[ibin][ir],weight);
                                my_hists->bkgDiffEtaBinJS[ibin][dEtaBin][ir]->Fill(smpt[ibin], smbkgrho[ibin][ir],weight);
                            }
                            // cout <<"ir =" <<ir <<"ibin =" <<ibin <<"jet pt ="<<jet_pt <<" smeared =" <<smpt[ibin]<<endl ;
                        }   //radius loop
                        
                    }//remove overlap region between jet and bkg cone
                    my_hists->NJets[ibin]->Fill(1, weight);  

                }  //centrality bins for smearing 
//                for(Int_t ibin = 0 ; ibin <nbin ; ibin++){
//                 } //centrality loop when smeared pp
            } //finish the smearing part        
            else {
                if(my_hists->IsMC==kTRUE){
                my_hists->genjetpt[curr_bin]->Fill(jet_pt,genjet_pt, weight);
                my_hists->genjeteta[curr_bin]->Fill(jet_eta, genjet_eta, weight);
                my_hists->genjetphi[curr_bin]->Fill(jet_phi, genjet_phi, weight);
                my_hists->genptRatio[curr_bin]->Fill(genjet_pt,jet_pt/genjet_pt, weight);
                my_hists->genphiRatio[curr_bin]->Fill(jet_pt,genjet_phi-jet_phi, weight);
                my_hists->genetaRatio[curr_bin]->Fill(jet_pt,genjet_eta-jet_eta, weight);
                }
                if(jet_pt<leadingjetcut) continue ;   
                my_hists->jetpt[curr_bin]->Fill(jet_pt, weight);
                my_hists->jeteta[curr_bin]->Fill(jet_pt, jet_eta, weight);
                my_hists->jetphi[curr_bin]->Fill(jet_pt, jet_phi, weight);

                for(Int_t ipt = 0 ; ipt <nptbin; ipt++){
                    if(jet_pt>=pt[ipt] && jet_pt <pt[ipt+1]) ptBin = ipt ;
                }
                my_hists->NEvents[curr_bin]->Fill(1, weight);
                if(TMath::Abs(my_ct->jteta[j4i])>etalimit) { //continue ;
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
                    double GenAxisrho[nrbin] ;
                    double GenAxisbkgrho[nrbin] ;
                    
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
                        GenAxisrho[ir]=0. ;
                        GenAxisbkgrho[ir]=0. ;
                        for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                            sumTrkpt[ir][ipt]=0.;
                            bkgsumTrkpt[ir][ipt]=0.;
                        }
                        
                        
                    }
                    //for tracking correction factors        
                    float trkweight = 1.; 
                    Int_t trackcount = 0 ;
                    
                    for(int itr = 0 ; itr < my_tr->nTrk ; itr++){
                        double dr = 0.;
                        double bkg_dr =0.;
                        double tr_pt = my_tr->trkPt[itr];
                        double tr_phi = my_tr->trkPhi[itr];
                        double tr_eta = my_tr->trkEta[itr];
                        if(TMath::Abs(tr_eta)>2.4) continue ;
                        
                        //  if(my_tr->trkNHit[itr]<8 ||(my_tr->highPurity[itr]&&my_tr->trkAlgo[itr]==4.0)){
                        if(my_tr->trkAlgo[itr]<4 ||(my_tr->highPurity[itr])){
                            //            if(my_tr->highPurity[itr]==0) continue ;
                            //            if(my_tr->trkAlgo[itr]!=4.0) continue ;
                            // if(TMath::Abs(tr_eta)>tracketacut) continue ;
                            if(tr_pt<trackcut) continue ;
                            
                            trkweight=1. ;
                            if(doTrackCorrections){
                                if(corrMet=="Hist")trkweight = c->getTrackCorrection(itr);  
                                else trkweight = c->getTrackCorrectionPara(itr); 
                                //  if(TMath::Abs(trkweight-1.)<0.001)trkweight = c->getTrackCorrection(itr);  
                            }
                            
                            dr =TMath::Sqrt((tr_phi-jet_phi)*(tr_phi-jet_phi)+(tr_eta-jet_eta)*(tr_eta-jet_eta));
                            bkg_dr = TMath::Sqrt((tr_phi-bkg_phi)*(tr_phi-bkg_phi)+(tr_eta-bkg_eta)*(tr_eta-bkg_eta));
                            
                            if(TMath::Abs(jet_eta)>=etalimit)my_hists->deltaR->Fill(dr, weight); 
                            double jtrdphi = tr_phi-jet_phi ;
                            //            if(jtrdphi>TMath::Pi())jtrdphi-=TMath::TwoPi();
                            //            if(jtrdphi<=-TMath::Pi()) jtrdphi+=TMath::TwoPi();            
                            if(TMath::Abs(jet_eta)>=etalimit)my_hists->IncTrackDphi[curr_bin]->Fill(jet_pt,TMath::Abs(jtrdphi), weight);
                            if(dr<=conesize){  //for leading jet shape study                        
                                trackcount++;
                                if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetTrackDphi[curr_bin]->Fill(jet_pt,TMath::Abs(jtrdphi), weight*trkweight);                
                                if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetTrackDeta[curr_bin]->Fill(jet_pt,tr_eta-jet_eta, weight*trkweight);                
                                
                                //  if(TMath::Abs(jet_eta)>=etalimit)my_hists->jetfrag[curr_bin]->Fill(jet_pt, TMath::Log(jet_pt/tr_pt), weight*trkweight);
                                if(TMath::Abs(jet_eta)>=etalimit)my_hists->jetfrag[curr_bin]->Fill(jet_pt, TMath::Log((jet_pt*TMath::CosH(jet_eta))/(tr_pt*TMath::CosH(tr_eta)*TMath::Cos(dr))), weight*trkweight);
                                if(TMath::Abs(jet_eta)>=etalimit)my_hists->trackpt[curr_bin]->Fill(jet_pt, tr_pt, weight*trkweight);
                                if(TMath::Abs(jet_eta)>=etalimit)my_hists->trackphi[curr_bin]->Fill(jet_pt, tr_phi, weight*trkweight);
                                if(TMath::Abs(jet_eta)>=etalimit)my_hists->tracketa[curr_bin]->Fill(jet_pt, tr_eta, weight*trkweight);
                                if(TMath::Abs(jet_eta)>=etalimit)my_hists->trackdr[curr_bin][ptBin]->Fill(tr_pt,dr, weight*trkweight);
                                if(TMath::Abs(jet_eta)>=etalimit)my_hists->wttrackdr[curr_bin][ptBin]->Fill(tr_pt,dr, tr_pt*weight*trkweight);
                                sumpt+=tr_pt*trkweight ;
                                meaneta+=tr_pt*tr_eta*trkweight;
                                meanphi+=tr_pt*tr_phi*trkweight;  
                                for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                    //  if(tr_pt*trkweight>trkpt[ipt] && tr_pt*trkweight<=trkpt[ipt+1])jetsumTrk[ipt]+=tr_pt*trkweight ;
                                    if(tr_pt*trkweight<=trkpt[ipt+1])jetsumTrk[ipt]+=tr_pt*trkweight ;
                                }
                                
                                for(int ir = 0 ; ir <nrbin; ir++){
                                    if(dr>rbin[ir]&&dr<=rbin[ir+1]){
                                        rho[ir]+=tr_pt*trkweight ;
                                        //     cout<<" !!!! track index =" <<itr << " trk pt =" <<tr_pt << " tr phi =" <<tr_phi << " tr eta =" <<tr_eta  <<" dr =" <<dr <<" rho =" <<rho[ir]<<endl ;
                                        if(TMath::Abs(jet_eta)>=etalimit)my_hists->ChargePt[curr_bin][ir]->Fill(jet_pt, tr_pt, weight*trkweight);
                                        sumchpt[ir]+=tr_pt*trkweight ;
                                        charge[ir]++ ;
                                        for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                            //  if(tr_pt*trkweight>trkpt[ipt] && tr_pt*trkweight<=trkpt[ipt+1])sumTrkpt[ir][ipt]+=tr_pt*trkweight ;
                                            if(tr_pt*trkweight<=trkpt[ipt+1])sumTrkpt[ir][ipt]+=tr_pt*trkweight ;
                                        }
                                        
                                    }
                                    if(dr<=rbin[ir+1]) {
                                        psi[ir]+=tr_pt*trkweight ; 
                                    }
                                }  //radius loop for rho calculation
                            } //inside leading jet cone
                            else { //outside leading and subleading jet cone
                                if(bkg_dr<=conesize){
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->bkgtrackpt[curr_bin]->Fill(jet_pt, tr_pt, weight);
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->bkgtrackdr[curr_bin][ptBin]->Fill(tr_pt,bkg_dr, weight*trkweight);
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->bkgwttrackdr[curr_bin][ptBin]->Fill(tr_pt,bkg_dr, tr_pt*weight*trkweight);
                                    bkgsumpt+=tr_pt*trkweight ;
                                    for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                        //  if(tr_pt*trkweight>trkpt[ipt] && tr_pt*trkweight<=trkpt[ipt+1])bkgjetsumTrk[ipt]+=tr_pt*trkweight ;
                                        if(tr_pt*trkweight<=trkpt[ipt+1])bkgjetsumTrk[ipt]+=tr_pt*trkweight ;
                                    }
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetBkgTrackDphi[curr_bin]->Fill(jet_pt,TMath::Abs(jtrdphi), weight*trkweight);
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->JetBkgTrackDeta[curr_bin]->Fill(jet_pt,tr_eta-jet_eta, weight*trkweight);
                                    //   if(TMath::Abs(jet_eta)>=etalimit)my_hists->jetbkgfrag[curr_bin]->Fill(jet_pt, TMath::Log(jet_pt/tr_pt), weight*trkweight);
                                    if(TMath::Abs(jet_eta)>=etalimit)my_hists->jetbkgfrag[curr_bin]->Fill(jet_pt, TMath::Log((jet_pt*TMath::CosH(jet_eta))/(tr_pt*TMath::CosH(tr_eta)*TMath::Cos(bkg_dr))), weight*trkweight);
                                    for(int ir = 0 ; ir <nrbin; ir++){
                                        if(bkg_dr>rbin[ir]&&bkg_dr<=rbin[ir+1]){
                                            if(TMath::Abs(jet_eta)>=etalimit)my_hists->bkgChargePt[curr_bin][ir]->Fill(jet_pt, tr_pt, weight*trkweight);
                                            bkgrho[ir]+=tr_pt*trkweight ;
                                            bkgsumchpt[ir]+=tr_pt*trkweight ;
                                            bkgcharge[ir]++ ;  
                                            for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                                //  if(tr_pt*trkweight>trkpt[ipt] && tr_pt*trkweight<=trkpt[ipt+1])bkgsumTrkpt[ir][ipt]+=tr_pt*trkweight ;
                                                if(tr_pt*trkweight<=trkpt[ipt+1])bkgsumTrkpt[ir][ipt]+=tr_pt*trkweight ;
                                            }
                                            
                                        }
                                        if(bkg_dr<=rbin[ir+1]) {
                                            bkgpsi[ir]+=tr_pt*trkweight ;                                       
                                        } //fill bkg subleading JS
                                    }  //radius loop for leading jet backgound 
                                }  //bkg cone for leading jet bkg
                            }
                            if(my_hists->IsMC==kTRUE&&DoGenAna){
                                double GenAxisdr =TMath::Sqrt((tr_phi-genjet_phi)*(tr_phi-genjet_phi)+(tr_eta-genjet_eta)*(tr_eta-genjet_eta));
                                double GenAxisbkgdr =TMath::Sqrt((tr_phi-genbkg_phi)*(tr_phi-genbkg_phi)+(tr_eta-genbkg_eta)*(tr_eta-genbkg_eta));
                                
                                if(dr<=conesize&&TMath::Abs(jet_eta)>=etalimit)my_hists->jetaxisRes[curr_bin]->Fill(jet_pt, (GenAxisdr-dr)/dr, weight);
                                if(GenAxisdr<=conesize){
                                    for(int ir = 0 ; ir <nrbin; ir++){
                                        if(GenAxisdr>rbin[ir]&&GenAxisdr<=rbin[ir+1]){
                                            GenAxisrho[ir]+=tr_pt*trkweight ;
                                        }
                                    }  //radius loop to fill leading jet JS histogram
                                } //tracks inside for leading jet
                                if(GenAxisbkgdr<=conesize){
                                    trkweight=1. ;
                                    //                        if(doTrackCorrections)trkweight = c->getTrackCorrection(itr);        
                                    if(doTrackCorrections){
                                        if(corrMet=="Hist")trkweight = c->getTrackCorrection(itr);  
                                        else trkweight = c->getTrackCorrectionPara(itr); 
                                    }
                                    for(int ir = 0 ; ir <nrbin; ir++){
                                        if(GenAxisbkgdr>rbin[ir]&&GenAxisbkgdr<=rbin[ir+1]){
                                            GenAxisbkgrho[ir]+=tr_pt*trkweight ;
                                        }
                                    }  //radius loop to fill leading jet bkg JS histogram
                                } //leading jet with gen axis background JS
                            } //fill only when MC    
                        } //track selection  
                    } //track loop
                    //     if(TMath::Abs(jet_eta)>=etalimit){ //remove overlap region for ER bkg
                    my_hists->SumptJetPt[curr_bin]->Fill(jet_pt,sumpt/jet_pt, weight);
                    if(DoSumPtNorm){
                        if(sumpt==0.) continue ;
                    }
                    for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                        my_hists->jetPtSumTrk[curr_bin][ipt]->Fill(jet_pt, jetsumTrk[ipt], weight);
                        my_hists->bkgjetPtSumTrk[curr_bin][ipt]->Fill(jet_pt, bkgjetsumTrk[ipt], weight);
                        
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
                            my_hists->jetdrSumTrkPt[curr_bin][ptBin][ipt]->Fill(rbin[ir]+0.025, sumTrkpt[ir][ipt]*weight);
                            my_hists->bkgjetdrSumTrkPt[curr_bin][ptBin][ipt]->Fill(rbin[ir]+0.025, bkgsumTrkpt[ir][ipt]*weight);
                        }
                        
                        my_hists->ChargeMult[curr_bin][ir]->Fill(jet_pt, charge[ir], weight);
                        my_hists->SumChPt[curr_bin][ir]->Fill(jet_pt, sumchpt[ir],weight); 
                        my_hists->bkgChargeMult[curr_bin][ir]->Fill(jet_pt, bkgcharge[ir],weight);
                        my_hists->bkgSumChPt[curr_bin][ir]->Fill(jet_pt, bkgsumchpt[ir],weight); 
                        my_hists->DiffJS[curr_bin][ir]->Fill(jet_pt, rho[ir],weight);
                        my_hists->IntJS[curr_bin][ir]->Fill(jet_pt, psi[ir],weight);
                        // rho/=deltacone ;
                        my_hists->bkgDiffJS[curr_bin][ir]->Fill(jet_pt, bkgrho[ir],weight);
                        //   my_hists->IntJS[curr_bin]->Fill(jet_pt, 1-psi);
                        my_hists->bkgIntJS[curr_bin][ir]->Fill(jet_pt, bkgpsi[ir],weight);
                        if(my_hists->IsMC==kTRUE&&DoGenAna){ 
                            if(DoSumPtNorm){
                                GenAxisrho[ir]/=sumpt;
                                GenAxisbkgrho[ir]/=sumpt;
                            }
                            else {
                                GenAxisrho[ir]/=jet_pt;
                                GenAxisbkgrho[ir]/=jet_pt;
                            }
                            my_hists->GenAxisDiffJS[curr_bin][ir]->Fill(jet_pt, GenAxisrho[ir],weight);
                            my_hists->GenAxisbkgDiffJS[curr_bin][ir]->Fill(jet_pt, GenAxisbkgrho[ir],weight);
                        }
                        if(jet_pt>=100 &&jet_pt<120) {            
                            //   my_hists->DiffJSPt80_100[curr_bin][ir]->Fill(rho[ir]);
                            my_hists->DiffJSPt80_100[curr_bin]->Fill(rbin[ir]+0.025, rho[ir]*weight);
                            my_hists->bkgDiffJSPt80_100[curr_bin]->Fill(rbin[ir]+0.025, bkgrho[ir]*weight);
                            my_hists->DiffJSPt80_100Rbin[curr_bin][ir]->Fill(rho[ir],weight);
                            my_hists->bkgDiffJSPt80_100Rbin[curr_bin][ir]->Fill(bkgrho[ir],weight);
                        }
                        if(dEtaBin!=-1){        
                            my_hists->DiffEtaBinJS[curr_bin][dEtaBin][ir]->Fill(jet_pt, rho[ir],weight);
                            my_hists->bkgDiffEtaBinJS[curr_bin][dEtaBin][ir]->Fill(jet_pt, bkgrho[ir],weight);
                        }
                    } //radius bin loop            
                    my_hists->NJets[curr_bin]->Fill(1, weight);   
                    if(bkgsumpt>sumpt)my_hists->RareEvt->Fill(offSel->hiBin, weight);
                }  //remove overlape eta regions from jet cone and ER bkg cone
            } //no smearing done
            if(my_hists->IsMC==kTRUE&&DoGenAna){
                if(TMath::Abs(genjet_eta)>=etalimit){
                    for(int ir = 0 ; ir <nrbin ; ir++){
                        gencharge[ir]=0;
                        gensumchpt[ir]=0. ;
                        genrho[ir]=0. ;
                        genpsi[ir]= 0.;
                        genbkgcharge[ir]=0 ;
                        genbkgsumchpt[ir]=0. ;
                        genbkgrho[ir]=0.;
                        genbkgpsi[ir]=0. ;
                        
                        for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                            gensumTrkpt[ir][ipt]=0.;
                            genbkgsumTrkpt[ir][ipt]=0.;
                        }
                        
                    }
                    //using the sim track to calculate the tracking efficiency 
                    for(int ipart = 0 ; ipart < my_tr->nParticle ; ipart++){ //sim track loop 
                        double gen_pt = my_tr->pPt[ipart];
                        double gen_phi = my_tr->pPhi[ipart];
                        double gen_eta = my_tr->pEta[ipart];  
                        if(gen_pt<trackcut)continue ;
                        if(TMath::Abs(gen_eta)>2.4)continue ;
                        //   cout <<"gen pt =" <<gen_pt <<endl ;
                        double gendr =TMath::Sqrt((gen_phi-jet_phi)*(gen_phi-jet_phi)+(gen_eta-jet_eta)*(gen_eta-jet_eta));
                        double genbkg_dr = TMath::Sqrt((gen_phi-bkg_phi)*(gen_phi-bkg_phi)+(gen_eta-bkg_eta)*(gen_eta-bkg_eta));
                        //   double gendr =TMath::Sqrt((gen_phi-genjet_phi)*(gen_phi-genjet_phi)+(gen_eta-genjet_eta)*(gen_eta-genjet_eta));
                        //    double genbkg_dr = TMath::Sqrt((gen_phi-genbkg_phi)*(gen_phi-genbkg_phi)+(gen_eta-genbkg_eta)*(gen_eta-genbkg_eta));
                        if(gendr<=conesize)my_hists->genpartdr[curr_bin][ptBin]->Fill(gen_pt, gendr, weight);
                        if(genbkg_dr<=conesize)my_hists->genbkgpartdr[curr_bin][ptBin]->Fill(gen_pt, genbkg_dr, weight);
                        if(gendr<=conesize)my_hists->wtgenpartdr[curr_bin][ptBin]->Fill(gen_pt, gendr,gen_pt*weight);
                        if(genbkg_dr<=conesize)my_hists->wtgenbkgpartdr[curr_bin][ptBin]->Fill(gen_pt, genbkg_dr, gen_pt*weight);
                        if(my_tr->pNRec[ipart]>0&&(my_tr->mtrkAlgo[ipart]<4||(my_tr->mtrkQual[ipart]))){
                            //   if(my_tr->pNRec[ipart]>0){  //no quality cuts                  
                            if(gendr<=conesize){
                                my_hists->matchpartdr[curr_bin][ptBin]->Fill(gen_pt, gendr,weight);
                                my_hists->wtmatchpartdr[curr_bin][ptBin]->Fill(gen_pt, gendr, gen_pt*weight);
                                gensumpt+=gen_pt ;
                                if(TMath::Abs(genjet_eta)>=etalimit)my_hists->genpartpt[curr_bin]->Fill(jet_pt,gen_pt, weight);
                                //      my_hists->genpartpt[curr_bin]->Fill(jet_pt,gen_pt, weight);
                                if(TMath::Abs(genjet_eta)>=etalimit)my_hists->genpartphi[curr_bin]->Fill(jet_pt,gen_phi, weight);
                                if(TMath::Abs(genjet_eta)>=etalimit)my_hists->genparteta[curr_bin]->Fill(jet_pt,gen_eta, weight);
                                for(int ir = 0 ; ir <nrbin ; ir++){
                                    if(gendr>rbin[ir]&&gendr<=rbin[ir+1]){
                                        genrho[ir]+=gen_pt ;
                                        if(TMath::Abs(genjet_eta)>=etalimit)my_hists->genChargePt[curr_bin][ir]->Fill(jet_pt, gen_pt, weight);
                                        gensumchpt[ir]+=gen_pt ;
                                        gencharge[ir]++ ;
                                        for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                            if(gen_pt<=trkpt[ipt+1])gensumTrkpt[ir][ipt]+=gen_pt ;
                                        }
                                        
                                    }
                                    if(gendr<=rbin[ir+1]) {
                                        genpsi[ir]+=gen_pt ; 
                                    }
                                }  //radius loop
                                
                            } //inside jet cone
                            if(genbkg_dr<=conesize){
                                my_hists->matchbkgpartdr[curr_bin][ptBin]->Fill(gen_pt, genbkg_dr, weight); 
                                my_hists->wtmatchbkgpartdr[curr_bin][ptBin]->Fill(gen_pt, genbkg_dr, gen_pt*weight); 
                                genbkgsumpt+=gen_pt ;
                                if(TMath::Abs(genjet_eta)>=etalimit)my_hists->bkgpartpt[curr_bin]->Fill(jet_pt,gen_pt, weight);
                                if(TMath::Abs(genjet_eta)>=etalimit)my_hists->bkgparteta[curr_bin]->Fill(gen_pt,gen_eta, weight);
                                for(int ir = 0 ; ir <nrbin ; ir++){
                                    if(genbkg_dr>rbin[ir]&&genbkg_dr<=rbin[ir+1]){
                                        genbkgrho[ir]+=gen_pt ;
                                        if(TMath::Abs(genjet_eta)>=etalimit)my_hists->genbkgChargePt[curr_bin][ir]->Fill(jet_pt, gen_pt,weight);
                                        genbkgsumchpt[ir]+=gen_pt ;
                                        genbkgcharge[ir]++ ;   
                                        for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                                            if(gen_pt<=trkpt[ipt+1])genbkgsumTrkpt[ir][ipt]+=gen_pt ;
                                        }
                                    }
                                    if(genbkg_dr<=rbin[ir+1]) {
                                        genbkgpsi[ir]+=gen_pt ;
                                    }
                                }  //radius loop for JS calculation  
                            } //inside bkg cone
                            
                        } //gen matched track selection
                    }  //sim track loop
                    if(DoSumPtNorm){
                        if(gensumpt==0.) continue ;
                    }
                    my_hists->genSumPt[curr_bin]->Fill(genjet_pt,gensumpt, weight);
                    
                    for(int ir = 0 ; ir <nrbin; ir++){  
                        if(DoSumPtNorm){
                            genpsi[ir]/=gensumpt;
                            genrho[ir]/=gensumpt;
                            genbkgpsi[ir]/=gensumpt;
                            genbkgrho[ir]/=gensumpt;
                            gensumchpt[ir]/=gensumpt ;
                            genbkgsumchpt[ir]/=gensumpt ;
                        }
                        else {
                            /*                        genpsi[ir]/=jet_pt;
                             genrho[ir]/=jet_pt;
                             genbkgpsi[ir]/=jet_pt;
                             genbkgrho[ir]/=jet_pt;
                             gensumchpt[ir]/=jet_pt ;
                             genbkgsumchpt[ir]/=jet_pt ;
                             */
                            genpsi[ir]/=genjet_pt;
                            genrho[ir]/=genjet_pt;
                            genbkgpsi[ir]/=genjet_pt;
                            genbkgrho[ir]/=genjet_pt;
                            gensumchpt[ir]/=genjet_pt ;
                            genbkgsumchpt[ir]/=genjet_pt ;            
                        }
                        for(int ipt = 0 ; ipt <ntrkptbin; ipt++){
                            my_hists->genjetdrSumTrkPt[curr_bin][ptBin][ipt]->Fill(rbin[ir]+0.025, gensumTrkpt[ir][ipt]*weight);
                            my_hists->bkggenjetdrSumTrkPt[curr_bin][ptBin][ipt]->Fill(rbin[ir]+0.025, genbkgsumTrkpt[ir][ipt]*weight);
                        }
                        my_hists->genChargeMult[curr_bin][ir]->Fill(jet_pt, gencharge[ir], weight);
                        my_hists->genSumChPt[curr_bin][ir]->Fill(jet_pt, gensumchpt[ir], weight); 
                        my_hists->genbkgChargeMult[curr_bin][ir]->Fill(jet_pt, genbkgcharge[ir], weight);
                        my_hists->genbkgSumChPt[curr_bin][ir]->Fill(jet_pt, genbkgsumchpt[ir], weight); 
                        
                        my_hists->genDiffJS[curr_bin][ir]->Fill(jet_pt, genrho[ir], weight);
                        //       my_hists->genDiffJS[curr_bin][ir]->Fill(jet_pt, genrho[ir],weight);
                        my_hists->genIntJS[curr_bin][ir]->Fill(jet_pt, genpsi[ir],weight);
                        // rho/=deltacone ;
                        my_hists->genbkgDiffJS[curr_bin][ir]->Fill(jet_pt, genbkgrho[ir],weight);
                        //      my_hists->genbkgDiffJS[curr_bin][ir]->Fill(jet_pt, genbkgrho[ir],weight);
                        //   my_hists->IntJS[curr_bin]->Fill(jet_pt, 1-psi);
                        my_hists->genbkgIntJS[curr_bin][ir]->Fill(jet_pt, genbkgpsi[ir],weight);
                    }  //radius loop in order to fill jet shape histograms   
                } //eta limit cut for ER bkg
                
            }//MC sim track loop for efficiency study   
        
        //   if(sumpt==0) cout<<"evt ="<<evi <<"jet pt =" <<jet_pt<<" without tracks" <<endl ;
        
//            if(my_hists->IsMC==kTRUE&&DoGenAna){ 
//                
//                genbkg_phi = genjet_phi;
//                genbkg_eta = -genjet_eta;
//                for(int ir = 0 ; ir <nrbin ; ir++){
//                    gencharge[ir]=0;
//                    gensumchpt[ir]=0. ;
//                    genrho[ir]=0. ;
//                    genpsi[ir]= 0.;
//                    genbkgcharge[ir]=0 ;
//                    genbkgsumchpt[ir]=0. ;
//                    genbkgrho[ir]=0.;
//                    genbkgpsi[ir]=0. ;
//                }
//                // for(int ipart = 0 ; ipart < my_tr->nParticle ; ipart++){ //sim track loop
//                for(int ipart = 0 ; ipart < my_GenPart->mult ; ipart++){//gen particle loop
//                    double gen_pt = my_GenPart->pt[ipart];
//                    double gen_phi = my_GenPart->phi[ipart];
//                    double gen_eta = my_GenPart->eta[ipart];
//                    int chg = my_GenPart->chg[ipart];
//                    if(chg==0) continue ;
//                    if(TMath::Abs(gen_eta)>2.4) continue ;
//                    //cout <<"charge =" <<chg <<endl ;
//                   //  int sube = my_GenPart->sube[ipart];
//                    // if(sube>0) continue ; //only PYTHIA signale particles
//                    // if(sube>0)cout <<"sube =" <<sube <<endl ;
//                    //                double gen_pt = my_tr->pPt[ipart];
//                    //                double gen_phi = my_tr->pPhi[ipart];
//                    //                double gen_eta = my_tr->pEta[ipart];  
//                    if(gen_pt<trackcut)continue ;
//                    double gendr =TMath::Sqrt((gen_phi-genjet_phi)*(gen_phi-genjet_phi)+(gen_eta-genjet_eta)*(gen_eta-genjet_eta));
//                    //     double gendr =TMath::Sqrt((gen_phi-jet_phi)*(gen_phi-jet_phi)+(gen_eta-jet_eta)*(gen_eta-jet_eta));
//                    //                double gendr =TMath::Sqrt((gen_phi-jet_phi)*(gen_phi-jet_phi)+(gen_eta-jet_eta)*(gen_eta-jet_eta));
//                    //                double gendr2 =TMath::Sqrt((gen_phi-subjet_phi)*(gen_phi-subjet_phi)+(gen_eta-subjet_eta)*(gen_eta-subjet_eta));
//                    if(gendr <=conesize){
//                        gensumpt+=gen_pt ;
//                        if(TMath::Abs(genjet_eta)>=etalimit)my_hists->genpartpt[curr_bin]->Fill(jet_pt,gen_pt, weight);
//                        //      my_hists->genpartpt[curr_bin]->Fill(jet_pt,gen_pt, weight);
//                        if(TMath::Abs(genjet_eta)>=etalimit)my_hists->genpartphi[curr_bin]->Fill(jet_pt,gen_phi, weight);
//                        if(TMath::Abs(genjet_eta)>=etalimit)my_hists->genparteta[curr_bin]->Fill(jet_pt,gen_eta, weight);
//                        for(int ir = 0 ; ir <nrbin ; ir++){
//                            if(gendr>rbin[ir]&&gendr<=rbin[ir+1]){
//                                genrho[ir]+=gen_pt ;
//                                if(TMath::Abs(genjet_eta)>=etalimit)my_hists->genChargePt[curr_bin][ir]->Fill(jet_pt, gen_pt, weight);
//                                gensumchpt[ir]+=gen_pt ;
//                                gencharge[ir]++ ;
//                                
//                            }
//                            if(gendr<=rbin[ir+1]) {
//                                genpsi[ir]+=gen_pt ; 
//                            }
//                        }  //radius loop
//                        
//                    } //tracks inside for jet
//                    else {
//                        double genbkg_dr = TMath::Sqrt((gen_phi-genbkg_phi)*(gen_phi-genbkg_phi)+(gen_eta-genbkg_eta)*(gen_eta-genbkg_eta));  
//                        //                        double genbkg_dr = TMath::Sqrt((gen_phi-bkg_phi)*(gen_phi-bkg_phi)+(gen_eta-bkg_eta)*(gen_eta-bkg_eta));  
//                        if(genbkg_dr<=conesize){
//                            genbkgsumpt+=gen_pt ;
//                            if(TMath::Abs(genjet_eta)>=etalimit)my_hists->bkgpartpt[curr_bin]->Fill(jet_pt,gen_pt, weight);
//                            if(TMath::Abs(genjet_eta)>=etalimit)my_hists->bkgparteta[curr_bin]->Fill(gen_pt,gen_eta, weight);
//                            for(int ir = 0 ; ir <nrbin ; ir++){
//                                if(genbkg_dr>rbin[ir]&&genbkg_dr<=rbin[ir+1]){
//                                    genbkgrho[ir]+=gen_pt ;
//                                    if(TMath::Abs(genjet_eta)>=etalimit)my_hists->genbkgChargePt[curr_bin][ir]->Fill(jet_pt, gen_pt, weight);
//                                    genbkgsumchpt[ir]+=gen_pt ;
//                                    genbkgcharge[ir]++ ;   
//                                }
//                                if(genbkg_dr<=rbin[ir+1]) {
//                                    genbkgpsi[ir]+=gen_pt ;
//                                }
//                            }  //radius loop for JS calculation
//                            
//                        }  // inside bkg cone   
//                    } //outside leading jet
//                    //     }
//                } //gen particle loop 
//                if(TMath::Abs(genjet_eta)>=etalimit){
//                    if(DoSumPtNorm && gensumpt==0.) continue ;
//                    //    if(gensumpt==0.|| genbkgsumpt>gensumpt) continue ;
//                    //      
//                    my_hists->genSumPt[curr_bin]->Fill(genjet_pt,gensumpt, weight);
//                    
//                    //    if(genbkgsumpt>=gensumpt || genbkgsumsubpt>=gensumsubpt) continue ;
//                    for(int ir = 0 ; ir <nrbin; ir++){  
//                        if(DoSumPtNorm){
//                            genpsi[ir]/=gensumpt;
//                            genrho[ir]/=gensumpt;
//                            genbkgpsi[ir]/=gensumpt;
//                            genbkgrho[ir]/=gensumpt;
//                            gensumchpt[ir]/=gensumpt ;
//                            genbkgsumchpt[ir]/=gensumpt ;
//                        }
//                        else {
//                            genpsi[ir]/=genjet_pt;
//                            genrho[ir]/=genjet_pt;
//                            genbkgpsi[ir]/=genjet_pt;
//                            genbkgrho[ir]/=genjet_pt;
//                            gensumchpt[ir]/=genjet_pt ;
//                            genbkgsumchpt[ir]/=genjet_pt ;
//                            
//                        }
//                        my_hists->genChargeMult[curr_bin][ir]->Fill(jet_pt, gencharge[ir], weight);
//                        my_hists->genSumChPt[curr_bin][ir]->Fill(jet_pt, gensumchpt[ir],weight); 
//                        my_hists->genbkgChargeMult[curr_bin][ir]->Fill(jet_pt, genbkgcharge[ir],weight);
//                        my_hists->genbkgSumChPt[curr_bin][ir]->Fill(jet_pt, genbkgsumchpt[ir],weight); 
//                        
//                        my_hists->genDiffJS[curr_bin][ir]->Fill(jet_pt, genrho[ir],weight);
//                        //       my_hists->genDiffJS[curr_bin][ir]->Fill(jet_pt, genrho[ir],weight);
//                        my_hists->genIntJS[curr_bin][ir]->Fill(jet_pt, genpsi[ir],weight);
//                        // rho/=deltacone ;
//                        my_hists->genbkgDiffJS[curr_bin][ir]->Fill(jet_pt, genbkgrho[ir],weight);
//                        //      my_hists->genbkgDiffJS[curr_bin][ir]->Fill(jet_pt, genbkgrho[ir],weight);
//                        //   my_hists->IntJS[curr_bin]->Fill(jet_pt, 1-psi);
//                        my_hists->genbkgIntJS[curr_bin][ir]->Fill(jet_pt, genbkgpsi[ir],weight);
//                    }  //radius loop in order to fill jet shape histograms   
//                } //eta limit cut for ER bkg
//
//            } //it is MC and gen level 

        } //jet loop
        //  cout << "still working222222\n";
        my_hists->Vertex->Fill(offSel->vz, weight);
        my_hists->CenBin->Fill(offSel->hiBin,weight); 
    }  ///event loop
    
    my_hists->Write();
    //   my_hists->Delete();
    //  delete my_hists;
    std::cout << "working done\n";
}




