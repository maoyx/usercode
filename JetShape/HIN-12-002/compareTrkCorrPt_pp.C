#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "hiForest.h"
#include "TLegend.h"
#include "TSystem.h"
#include "HisMath.C"
#include "commonUtility.h"

void compareTrkCorrPt_pp(
                           TString outdir="fig"
)
{
   TH1::SetDefaultSumw2();
   gSystem->mkdir(outdir,kTRUE);
   float xmin=1,xmax=179.9;
   TString title="pp";
   TString reftitle="PbPb";

    const bool SaveFile=kTRUE ;
   /////////////////////////////////////////////////////////////////////////////////////
   // Load Histograms
   /////////////////////////////////////////////////////////////////////////////////////
 //  HiForest * cpp = new HiForest("/net/hisrv0001/home/zhukova/scratch/HIHighPt/forest/pthat200/mergedFile.root","forest",1);
   HiForest * cpp = new HiForest("/mnt/hadoop/cms/store/user/dgulhan/ppHiIterativeTrack/P01/prod24/Signal_Pythia_pt170/HiForest_v84_merged01/pt170_JEC_ppHiIterativeTrack_P01_prod24_v84_merged_forest_0.root","forest",1);
   cpp->doTrackCorrections = true;
   cpp->doTrackingSeparateLeadingSubleading = false;
   cpp->InitTree();   
   TrackingCorrections * trkCorr = cpp->trackCorrections[0];

   HiForest * cref = new HiForest("/mnt/hadoop/cms/store/user/yenjie/HiForest_v27/Dijet170_HydjetDrum_v27_mergedV1.root","forestref",0);
   cref->doTrackCorrections = true;
   cref->doTrackingSeparateLeadingSubleading = false;
   cref->InitTree();
   TrackingCorrections * trkCorrRef = cref->trackCorrections[0];


   cout << endl << "========= plot =========" << endl;
   Int_t etaPM=5.; // 7+2,-3 for |eta|<1.2, 7+5,-6 for full eta
   Float_t jetPtMin=0;
   Float_t jetPtMax=500;
   Int_t jetBegBin = trkCorr->jetBin_->FindBin(jetPtMin);
 //  Int_t jetEndBin = trkCorr->numJEtBins_;
   Int_t jetEndBin = trkCorr->jetBin_->FindBin(jetPtMax);;
   cout << Form("jet pt %.0f bin: ",jetPtMin) << jetBegBin << " to " << jetEndBin << endl;
   cout << "========================" << endl;

    string infpath=trkCorr->sample_[0]->GetName();
    TString src0=infpath.substr(infpath.find_last_of('/')+1);
    TString src =src0(src0.First("Iter")+4,14);
    cout <<" src =" <<src <<endl ; 
 //   src.ReplaceAll(".root","");
    TString tag = src+"_"+trkCorr->trkCorrModule_+Form("_vs_Pt_%s_%s_jet%.0f_%.0f_ieta%d_wts%d",title.Data(),reftitle.Data(),jetPtMin,jetPtMax, etaPM,trkCorr->weightSamples_);

 //  TString tag = trkCorr->trkCorrModule_+Form("_vs_Pt_%s_%s_jet%.0f_ieta%d_wts%d",title.Data(),reftitle.Data(),jetPtMin,etaPM,trkCorr->weightSamples_);
   
   // Get Eff/fake histograms
   int numCentBin=trkCorr->numCentBins_;
   int numCentBinRef=trkCorrRef->numCentBins_;
    cout <<"ppncen=" << numCentBin <<"  ncen=" << numCentBinRef <<endl ;
	TH1D * vhCorrPtRef[2][5], *vhCorrPt[2][5];
	Int_t colors[10] = {kBlack,kRed,kYellow+2,kGreen+2,kBlue};
   Int_t styles[2] = {kFullCircle,kOpenCircle};

   int icent=0;
 //  int icentRef=numCentBinRef-1;

	for (Int_t lv=0; lv<2; ++lv) {
		for (Int_t c=0; c<numCentBin; ++c) {
			vhCorrPt[lv][c] = (TH1D*) trkCorr->InspectCorr(lv,c,c,jetBegBin,jetEndBin,2,7-etaPM-1,7+etaPM);
			handsomeTH1(vhCorrPt[lv][c],kBlue,1,kOpenCircle);
         vhCorrPt[lv][icent]->SetAxisRange(xmin,xmax,"X");
         vhCorrPt[lv][icent]->SetAxisRange(0,1,"Y");
		}
		for (Int_t c=0; c<numCentBinRef; ++c) {
			vhCorrPtRef[lv][c] = (TH1D*) trkCorrRef->InspectCorr(lv,c,c,jetBegBin,jetEndBin,2,7-etaPM-1,7+etaPM);
			handsomeTH1(vhCorrPtRef[lv][c],colors[c]);
         vhCorrPtRef[lv][icent]->SetAxisRange(xmin,xmax,"X");
         vhCorrPtRef[lv][icent]->SetAxisRange(0,1,"Y");
		}
	}
   
	TCanvas * cEff = new TCanvas("cEff","cEff",500,500);
   cEff->SetLogx();
   vhCorrPt[0][icent]->SetTitle(";Track p_{T} (GeV/c);A #times #epsilon");
   vhCorrPt[0][icent]->SetTitleOffset(1.2);
   vhCorrPt[0][icent]->SetTitleSize(0.055);
	vhCorrPt[0][icent]->Draw("E");
	vhCorrPt[1][icent]->Draw("sameE");
	vhCorrPtRef[0][0]->Draw("sameE");
	vhCorrPtRef[1][0]->Draw("sameE");
                for (Int_t c=0; c<numCentBinRef; ++c) {
	vhCorrPtRef[0][c]->Draw("sameE");
	vhCorrPtRef[1][c]->Draw("sameE");
}
   TLegend *leg0 = new TLegend(0.16,0.84,0.46,0.92);
   leg0->SetFillStyle(0);
   leg0->SetBorderSize(0);
   leg0->SetTextSize(0.04);
   leg0->AddEntry(vhCorrPt[0][0],"PYTHIA+HYDJET","");
   if (etaPM==5)   leg0->AddEntry(vhCorrPt[0][0],Form("Track |#eta|<2.4"),"");
   if (etaPM==2) leg0->AddEntry(vhCorrPt[0][0],Form("Track |#eta|<1.2"),"");
	leg0->Draw();
   TLine * l = new TLine(xmin,1,xmax,1);
   l->SetLineStyle(2);
   l->Draw();
	
   TLegend *leg = new TLegend(0.34,0.32,0.60,0.48);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextSize(0.035);
   leg->AddEntry(vhCorrPt[0][icent],title,"p");
    if (numCentBinRef==2) {
        leg->AddEntry(vhCorrPtRef[0][0],reftitle+ "0-30%","p");
        leg->AddEntry(vhCorrPtRef[0][1],reftitle+ "30-100%","p");
    } else {
        leg->AddEntry(vhCorrPtRef[0][0],reftitle+ "0-10%","p");
        leg->AddEntry(vhCorrPtRef[0][1],reftitle+ "10-30%","p");
        leg->AddEntry(vhCorrPtRef[0][2],reftitle+ "30-50%","p");
        leg->AddEntry(vhCorrPtRef[0][3],reftitle+ "50-70%","p");
        leg->AddEntry(vhCorrPtRef[0][4],reftitle+ "70-100%","p");
    }
   leg->Draw();
   
	drawText("CMS Simulation",0.64,0.89);
	drawText("Fake Rate",0.69,0.26);
   
   cEff->Print(outdir+"/"+tag+".gif");
   cEff->Print(outdir+"/"+tag+".pdf");

	TCanvas * cJet = new TCanvas("cJet","cJet",500,500);
	cJet->SetLogy();

	trkCorr->vhPtHat[1][icent]->Draw("E");
        
       for (Int_t c=0; c<numCentBinRef; ++c) {
	trkCorrRef->vhPtHat[1][c]->SetMarkerColor(kRed+c);
	trkCorrRef->vhPtHat[1][c]->SetMarkerStyle(kOpenCircle+c);
	trkCorrRef->vhPtHat[1][c]->Draw("same E");
        }

	TCanvas * cCent = new TCanvas("cCent","cCent",500,500);
   TH1D * hCent = (TH1D*)trkCorr->sample_[0]->Get("hCent");
   TH1D * hCentRef = (TH1D*)trkCorrRef->sample_[0]->Get("hCent");
	hCentRef->SetMarkerStyle(kOpenCircle);
	hCentRef->Scale(1./hCentRef->GetEntries());
	hCent->Scale(1./hCent->GetEntries());
	hCentRef->Draw("p");
	hCent->Draw("samep");

    if(SaveFile){
        TFile * outf = new TFile(Form("%s/%sTrkEffFake.root", outdir.Data(),tag.Data()), "CREATE");
        for (Int_t lv=0; lv<2; ++lv) {
        vhCorrPt[lv][icent]->Write();
            for (Int_t c=numCentBinRef-1; c>=0; --c) {
                vhCorrPtRef[lv][c]->Write();
            }
        }
        outf->Close();
    }
}
