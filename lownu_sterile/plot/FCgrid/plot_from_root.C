void plot_from_root()
{
  int i = 1;

  TFile *f = new TFile(Form("NuCut_dm%d.root",i));

  TH2D *h = (TH2D*)f->Get("h");
  h->GetXaxis()->SetTitle("U_{e4}^{2}");
  h->GetYaxis()->SetTitle("U_{m4}^{2}"); 

  double contours[2];
  contours[0] = 3;
  contours[1] = 5;

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  gStyle->SetNumberContours(999);

  TCanvas *c = new TCanvas("c","",800,600);
  gPad->SetGrid();
  c->SetLogx();
  c->SetLogy();
  //c->SetLogz();
  h->DrawCopy("colz");
/*
  h->SetContour(2,contours);
  h->Draw("cont3 same");
  h->SetLineColor(kRed);
*/
  c->SaveAs(Form("contour%d.png",i));
}
