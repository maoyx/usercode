

{

  // Make sure you calculate statistical errors correctly
  TH1::SetDefaultSumw2();

  // No histogram titles displayed. 
  TH1D* h1 = new TH1D("h1",";p_{T}^{CaloJet}/p_{T}^{PFJet};Events",50,0,2);
  TH1D* h2 = new TH1D("h2",";p_{T}^{CaloJet}/p_{T}^{PFJet};Events",50,0,2);

  TF1* fGaus = new TF1("f","gaus(0)",0,2);
  fGaus->SetParameter(0,1);
  fGaus->SetParameter(1,1);
  fGaus->SetParameter(2,0.2);
  h1->FillRandom("f",10000);

  fGaus->SetParameter(2,0.22);
  h2->FillRandom("f",10000);

  // Please use RED for MC, keep BLACK for DATA
  h2->SetLineColor(2);
  h2->SetMarkerSize(0);

  // Axis Titles always centered!
  h1->GetXaxis()->CenterTitle();
  h1->GetYaxis()->CenterTitle();
  h2->GetXaxis()->CenterTitle();
  h2->GetYaxis()->CenterTitle();

  // Please have square plots: makes it very convenient manipulating in slides
  TCanvas* c1 = new TCanvas("c1","",600,600);

  h2->SetMaximum(h2->GetMaximum()*1.5);
  h2->Draw("hist"); // MC should be always plotted as line histogram
  h2->Draw("same"); // In case of low MC statistics, show also MC histogram errors
  h1->Draw("same"); // DATA should always be plotted as points

  // No borders, no filling in legends
  // Make sure legend points are far from data points!
  TLegend *leg=new TLegend(0.60,0.80,0.94,0.92);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(43);
  leg->SetTextSize(25);

  leg->AddEntry(h1,"Data","p");
  leg->AddEntry(h2,"MC","l");

  // Add CMS Preliminary for PAS
  // If MC plots only, add "CMS Simulation"
  TLatex *cms = new TLatex(0.06,h2->GetMaximum()*0.92,"CMS Preliminary");
  cms->SetTextFont(63);
  cms->SetTextSize(22);
  cms->Draw();

  leg->Draw();


}

