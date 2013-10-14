


#include "hiForest.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TNtuple.h"
#include "TMath.h"

#include "TCut.h"
#include <string>

using namespace std;

void analyzeJetTracks(
		 const char* infname = "/d102/yjlee/hiForest2MC/Pythia80_HydjetDrum_mix01_HiForest2_v20.root",
                 const char* outname = "output.root"
		 ){


  TFile* outf = new TFile( outname, "recreate");
  TH2D* hpt2D = new TH2D("hpt2D",";p_{T}^{Lead} (GeV/c);p_{T}^{SubLead} (GeV/c)",100,0,1000,100,0,1000);


  HiForest * t = new HiForest(infname);
  // Here goes more up to date tracking correction instructions
  t->InitTree();

  int maxEvents = t->GetEntries();
  for(int iev = 0; iev < maxEvents; ++iev){
    if(iev%1000==0)cout<<"Processing entry : "<<iev<<" / "<<t->GetEntries()<<endl;
    t->GetEntry(iev);

    int evt = t->hlt.Event;
    int run = t->hlt.Run;

    //Jet Loop
    for(int j = 0; j < t->akPu3PF.nref; ++j){
      double jtpt =  t->akPu3PF.jtpt[j];

    }

    // Track Loop
    for(int i = 0; i < t->track.nTrk; ++i){
      double trkPt =  t->track.trkPt[i];
    }

    // SimTrack Loop
    for(int i = 0; i < t->track.nParticle; ++i){
      double trkPt =  t->track.pPt[i];
    }


    // GenParticle loop
    for(int i = 0; i < t->genparticle.nPar; ++i){
      double trkPt =  t->genparticle.et[i];
    }


  }

}


