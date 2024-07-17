void plot()
{
  const char data_path[] = "/exp/dune/app/users/qvuong/data/lownu/contours";
  
  TFile *f = new TFile(Form("%s/contour.root",data_path),"READ");

  TH2D *h0 = (TH2D*)f->Get("hdm1");
  TH2D *h1 = (TH2D*)f->Get("hdm10");
  TH2D *h2 = (TH2D*)f->Get("hdm100");

  double contours[1];
  contours[0] = 9; 

  h0->SetContour(1, contours);
  h1->SetContour(1, contours);
  h2->SetContour(1, contours);

  h0->SetStats(0);
  h1->SetStats(0);
  h2->SetStats(0);

  h0->GetXaxis()->SetRangeUser(1e-4,1.);
  h0->GetYaxis()->SetRangeUser(1e-4,1.);
  h1->GetXaxis()->SetRangeUser(1e-4,1.);
  h1->GetYaxis()->SetRangeUser(1e-4,1.);
  h2->GetXaxis()->SetRangeUser(1e-4,1.);
  h2->GetYaxis()->SetRangeUser(1e-4,1.);
  
  h0->SetLineColor(kBlack);
  h0->SetLineStyle(2);
  h0->SetLineWidth(2);
  h1->SetLineColor(kRed);
  h1->SetLineStyle(1);
  h1->SetLineWidth(2);
  h2->SetLineColor(kBlue);
  h2->SetLineStyle(9);
  h2->SetLineWidth(2);

  h0->SetTitle("3#sigma Contours");
  h0->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  h0->GetYaxis()->SetTitle("U_{e4}^{2}");

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  gStyle->SetNumberContours(999);

  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->SetLogx();
  c->SetLogy();
  c->SetLogz();
  h0->Draw("cont3 same");
  h1->Draw("cont3 same");
  h2->Draw("cont3 same");
  TLegend *lg = new TLegend(0.65,0.75,0.9,0.9);
  lg->AddEntry(h0,"#Deltam^{2} = 1.0");
  lg->AddEntry(h1,"#Deltam^{2} = 10.0");
  lg->AddEntry(h2,"#Deltam^{2} = 100.0");
  lg->Draw();
  c->SaveAs(Form("%s/3sigma.png",data_path));

}
