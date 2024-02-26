void plot1()
{
  double ct[2];
  ct[0] = 9; 
  ct[1] = 25; 
  
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  gStyle->SetNumberContours(999);

  for(int i=0; i<4; i++){ 
 
  TFile *f = new TFile(Form("contour100_%d.root",i),"READ");

  TH2D *h0 = (TH2D*)f->Get("h");


  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->SetLogx();
  c->SetLogy();
  c->SetLogz();
  h0->DrawCopy("colz");
  h0->SetContour(2, ct);
  h0->SetLineColor(kGreen);
  h0->SetLineWidth(2);
  h0->Draw("cont3 same");
  c->SaveAs(Form("contour100_Log_%d.png",i));
  
  //gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  }
}
