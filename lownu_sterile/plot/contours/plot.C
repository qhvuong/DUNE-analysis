void plot()
{
  
  TFile *f0 = new TFile("contour100_0.root","READ");
  TFile *f1 = new TFile("contour100_1.root","READ");
  TFile *f2 = new TFile("contour100_2.root","READ");
  TFile *f3 = new TFile("contour100_3.root","READ");

  TH2D *h0 = (TH2D*)f0->Get("h");
  TH2D *h1 = (TH2D*)f1->Get("h");
  TH2D *h2 = (TH2D*)f2->Get("h");
  TH2D *h3 = (TH2D*)f3->Get("h");

  double contours[1];
  contours[0] = 9; 

  h0->SetContour(1, contours);
  h1->SetContour(1, contours);
  h2->SetContour(1, contours);
  h3->SetContour(1, contours);
  
  h0->SetLineColor(kBlack);
  h0->SetLineStyle(2);
  h0->SetLineWidth(2);
  h1->SetLineColor(kGreen);
  h1->SetLineStyle(1);
  h1->SetLineWidth(2);
  h2->SetLineColor(kRed);
  h2->SetLineStyle(9);
  h2->SetLineWidth(2);
  h3->SetLineColor(kBlue);
  h3->SetLineStyle(10);
  h3->SetLineWidth(2);

  h0->SetTitle("3#sigma Contours");

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  gStyle->SetNumberContours(999);

  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->SetLogx();
  c->SetLogy();
  c->SetLogz();
  //h->DrawCopy("colz");
  h0->Draw("cont3 same");
  h1->Draw("cont3 same");
  h2->Draw("cont3 same");
  h3->Draw("cont3 same");
  TLegend *lg = new TLegend(0.55,0.70,0.9,0.9);
  lg->AddEntry(h0,"#Deltam^{2} = 1.0");
  lg->AddEntry(h1,"#Deltam^{2} = 5.0");
  lg->AddEntry(h2,"#Deltam^{2} = 10.0");
  lg->AddEntry(h3,"#Deltam^{2} = 100.0");
  lg->Draw();
  c->SaveAs("3sigma.png");

}
