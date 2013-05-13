#include "TCanvas.h"
#include "TStyle.h"
#include "TLine.h"
#include "/net/hidsk0001/d00/scratch/maoyx/JetShape/CMSSW_4_4_2_patch5/src/HiForest_V2_02_16/hiForest.h"
#include "TLegend.h"
#include "TSystem.h"
#include "HisMath.C"
#include "commonUtility.h"

void drawTrkCorrPtvCent(
                        TString outdir="fig"
                        )
{
    TH1::SetDefaultSumw2();
    gSystem->mkdir(outdir,kTRUE);
    float xmin=1,xmax=179.9;
    TString title="Iterative Tracking";
    const bool SaveFile=kTRUE ;
    /////////////////////////////////////////////////////////////////////////////////////
    // Load Histograms
    /////////////////////////////////////////////////////////////////////////////////////
//    HiForest * c = new HiForest("/net/hidsk0001/d00/scratch/yjlee/merge/v27/pthat170/Dijet170_HydjetDrum_v27_mergedV1.root");
    HiForest * c = new HiForest("/net/hisrv0001/home/zhukova/scratch/HIHighPt/forest/pthat120/mergedFile.root");
    c->doTrackCorrections = true;
    c->InitTree();
    
    TrackingCorrections * trkCorr = c->trackCorrections[0];
    
    cout << endl << "========= plot =========" << endl;
    Int_t etaPM=5; // 7 +2,-3 for |eta|<1.2, 7 =5,-6 for full eta
    Float_t jetPtMin=100;
    Int_t jetBegBin = trkCorr->jetBin_->FindBin(jetPtMin);
    Int_t jetEndBin = trkCorr->numJEtBins_;
    cout << Form("jet pt %.0f bin: ",jetPtMin) << jetBegBin << " to " << jetEndBin << endl;
    cout << "========================" << endl;
    bool doTestCorr = true;
    
    
    string infpath=trkCorr->sample_[0]->GetName();
    TString src=infpath.substr(infpath.find_last_of('/')+1);
    src.ReplaceAll(".root","");
    TString tag = src+"_"+trkCorr->trkCorrModule_+Form("_vs_Pt_jet%.0f_ieta%d_wts%d",jetPtMin,etaPM,trkCorr->weightSamples_);
    
    
    /////////////////////////////////////////////////////////////////////////////////////
    // Inspect Projection
    /////////////////////////////////////////////////////////////////////////////////////
    // Get Eff/fake histograms
    int numCentBin=trkCorr->numCentBins_;
    TH1D * vhCorrPtRef[2][5], *vhCorrPt[2][5];
    Int_t colors[10] = {kBlack,kRed,kYellow+2,kGreen+2,kBlue};
    Int_t styles[2] = {kFullCircle,kOpenCircle};
    for (Int_t lv=0; lv<2; ++lv) {
        for (Int_t c=0; c<numCentBin; ++c) {
            vhCorrPt[lv][c] = (TH1D*) trkCorr->InspectCorr(lv,c,c,jetBegBin,jetEndBin,2,7-etaPM-1,7+etaPM);
            cout << "lv: " << lv << " c: " << c << " " << vhCorrPt[lv][c] << endl;
            vhCorrPt[lv][c]->SetMarkerStyle(styles[lv]);
            handsomeTH1(vhCorrPt[lv][c],colors[c]);
            vhCorrPt[lv][c]->SetAxisRange(xmin,xmax,"X");
        }
    }
    
    // Draw Histograms
    TCanvas * cEff = new TCanvas("cEff","cEff",500,500);
    cEff->SetLogx();
    vhCorrPt[0][0]->SetAxisRange(0,1,"Y");
    vhCorrPt[0][0]->SetTitle(";Track p_{T} (GeV/c);A #times #epsilon");
    vhCorrPt[0][0]->SetTitleOffset(1.2);
    vhCorrPt[0][0]->SetTitleSize(0.055);
    vhCorrPt[0][0]->Draw("E");
    vhCorrPt[1][0]->Draw("sameE");
    for (Int_t lv=0; lv<2; ++lv) {
        for (Int_t c=numCentBin-1; c>=0; --c) {
            vhCorrPt[lv][c]->Draw("sameE");
        }
    }
    TLegend *leg0 = new TLegend(0.16,0.786,0.46,0.92);
    leg0->SetFillStyle(0);
    leg0->SetBorderSize(0);
    leg0->SetTextSize(0.04);
    leg0->AddEntry(vhCorrPt[0][0],"PYTHIA+HYDJET","");
    if (jetPtMin >= 40) leg0->AddEntry(vhCorrPt[0][0],Form("Jet p_{T} #geq %.0f GeV/c",jetPtMin),"");
    leg0->AddEntry(vhCorrPt[0][0],Form("Track %.1f < #eta < %.1f",trkCorr->etaBin_->GetBinLowEdge(7-etaPM-1), trkCorr->etaBin_->GetBinLowEdge(7+etaPM+1)),"");
    leg0->Draw();
    TLine * l = new TLine(xmin,1,xmax,1);
    l->SetLineStyle(2);
    l->Draw();
	
    TLegend *leg = new TLegend(0.34,0.25,0.56,0.55);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(vhCorrPt[0][0],title,"");
    if (numCentBin==2) {
        leg->AddEntry(vhCorrPt[0][0],"0-30%","p");
        leg->AddEntry(vhCorrPt[0][1],"30-100%","p");
    } else if (numCentBin==4) {
        leg->AddEntry(vhCorrPt[0][0],"0-10%","p");
        leg->AddEntry(vhCorrPt[0][1],"10-30%","p");
        leg->AddEntry(vhCorrPt[0][2],"30-50%","p");
        leg->AddEntry(vhCorrPt[0][3],"50-100%","p");
    }
    leg->Draw();
    
	drawText("CMS Simulation",0.64,0.89);
	drawText("Fake Rate",0.69,0.26);
    
    cEff->Print(outdir+"/"+tag+".gif");
    cEff->Print(outdir+"/"+tag+".pdf");
    
    if(SaveFile){
        TFile * outf = new TFile(Form("%s/%sTrkEffFake.root", outdir.Data(),src.Data()), "CREATE");
        for (Int_t lv=0; lv<2; ++lv) {
            for (Int_t c=numCentBin-1; c>=0; --c) {
                vhCorrPt[lv][c]->Write();
            }
        }
        outf->Close();  
    }
    /////////////////////////////////////////////////////////////////////////////////////
    // Inspect Events
    /////////////////////////////////////////////////////////////////////////////////////
	TCanvas * cPtHat = new TCanvas("cPtHat","cPtHat",1000,500);
	cPtHat->Divide(2,1);
	cPtHat->cd(1);
	gPad->SetLogy();
	trkCorr->vhPtHat[0][0]->SetMarkerStyle(kOpenCircle);
	trkCorr->vhPtHat[0][0]->SetMarkerColor(kRed);
	trkCorr->vhPtHat[1][0]->SetTitle(";#hat{p}_{T} (GeV/c);a.u.");
	normHist(trkCorr->vhPtHat[0][0],1,false);//,trkCorr->vhPtHat[0][0]->Integral()/trkCorr->vhPtHat[0][0]->Integral(24,200));
	normHist(trkCorr->vhPtHat[1][0],1,false);//,trkCorr->vhPtHat[1][0]->Integral()/trkCorr->vhPtHat[1][0]->Integral(24,200));
	trkCorr->vhPtHat[1][0]->SetAxisRange(0,500,"X");
	trkCorr->vhPtHat[1][0]->SetAxisRange(1e-7,1,"Y");
	trkCorr->vhPtHat[1][0]->Draw("E");
    // 	trkCorr->vhPtHat[0][0]->Draw("Esame");
    //    TLegend *legev0 = new TLegend(0.53,0.76,0.75,0.90);
    //    legev0->SetFillStyle(0);
    //    legev0->SetBorderSize(0);
    //    legev0->SetTextSize(0.035);
    //    legev0->AddEntry(trkCorr->vhPtHat[0][0],"Raw","p");
    //    legev0->AddEntry(trkCorr->vhPtHat[1][0],"Weighted","p");
    //    legev0->Draw();
	cPtHat->cd(2);
	gPad->SetLogx();
	gPad->SetLogy();
	normHist(trkCorr->hDen1DInsp,1,true);
	normHist(trkCorr->hNum1DInsp,1,true);
	trkCorr->hDen1DInsp->SetMarkerStyle(kOpenCircle);
	trkCorr->hDen1DInsp->SetTitle(";Track p_{T} (GeV/c);a.u.");
	trkCorr->hDen1DInsp->SetAxisRange(1.001,299.9,"X");
	trkCorr->hDen1DInsp->SetAxisRange(1e-13,1e1,"Y");
	trkCorr->hDen1DInsp->Draw("hist");
	trkCorr->hNum1DInsp->Draw("Esame");
	cPtHat->Print(outdir+"/"+tag+"_weighting.gif");
	cPtHat->Print(outdir+"/"+tag+"_weighting.pdf");
    
    /////////////////////////////////////////////////////////////////////////////////////
    // Inspect Each Bin
    /////////////////////////////////////////////////////////////////////////////////////
    TCanvas * cEff2D = new TCanvas("cEff2D","cEff2D",800,800);
	gPad->SetLogy();
    Double_t pt=8,eta=0,jet=101.  ;
    Int_t cBin = 0;
    cEff2D->Divide(2,2);
    cEff2D->cd(1);
    gPad->SetLogy();
    TH2D * hCorr2D = (TH2D*)trkCorr->InspectCorr(0,cBin,cBin,trkCorr->jetBin_->FindBin(jet),trkCorr->jetBin_->FindBin(jet));
    gPad->SetRightMargin(0.15);
    hCorr2D->SetAxisRange(0.5,119.9,"Y");
    hCorr2D->SetAxisRange(0,1,"Z");
    hCorr2D->Draw("colz");
    cEff2D->cd(2);
    gPad->SetLogy();
    TH2D * hCorr2DFak= (TH2D*)trkCorr->InspectCorr(1,cBin,cBin,trkCorr->jetBin_->FindBin(jet),trkCorr->jetBin_->FindBin(jet));
    gPad->SetRightMargin(0.15);
    hCorr2DFak->SetAxisRange(0.5,119.9,"Y");
    hCorr2DFak->SetAxisRange(0,1,"Z");
    hCorr2DFak->Draw("colz");
    cEff2D->cd(3);
    gPad->SetLogy();
    hCorr2D = (TH2D*)trkCorr->InspectCorr(0,cBin,cBin,trkCorr->jetBin_->FindBin(200.),trkCorr->jetBin_->FindBin(jet));
    gPad->SetRightMargin(0.15);
    hCorr2D->SetAxisRange(1,119.9,"Y");
    hCorr2D->SetAxisRange(0,1,"Z");
    hCorr2D->Draw("colz");
    cEff2D->cd(4);
    gPad->SetLogy();
    hCorr2DFak= (TH2D*)trkCorr->InspectCorr(1,cBin,cBin,trkCorr->jetBin_->FindBin(200.),trkCorr->jetBin_->FindBin(jet));
    gPad->SetRightMargin(0.15);
    hCorr2DFak->SetAxisRange(1,119.9,"Y");
    hCorr2DFak->SetAxisRange(0,1,"Z");
    hCorr2DFak->Draw("colz");
    cEff2D->Print(outdir+"/"+tag+"_2D.gif");
    
    /////////////////////////////////////////////////////////////////////////////////////
    // Test corr
    /////////////////////////////////////////////////////////////////////////////////////
    if (doTestCorr) {
        cout << "trk weight: " << trkCorr->GetCorr(pt,eta,jet,cBin) << endl;
        Double_t corr[4];
        for (Int_t i=1; i<=trkCorr->ptBin_->GetNbinsX(); ++i) {
            trkCorr->GetCorr(trkCorr->ptBin_->GetBinCenter(i),eta,jet,cBin,corr);
            cout << "trk pt: " << trkCorr->ptBin_->GetBinLowEdge(i) << " trk eff: " << corr[0] << " trk fak: " << corr[1] << endl;
        }
    }   
}
