#include <TCanvas.h>
#include <TH1D.h>
#include <TPaveText.h>
#include <DUNEStyle.h>
#include <TLegend.h>


// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}


void plot2D(TCanvas *c, TH2D *h, const char *name, const char *title)
{
  c->Clear();
  c->cd();
  c->SetLogx();
  dunestyle::CenterTitles(h);
  h->Draw("colz");
  dunestyle::Simulation();
  
  TPaveText *pt = new TPaveText(0.65, 0.78, 0.75, 0.88, "NDC");
  pt->SetFillStyle(0);
  pt->SetFillColor(0);  // Transparent background
  pt->SetTextColor(kBlack);
  pt->SetTextSize(0.08);
  pt->SetBorderSize(0); // No border

  pt->AddText(title);  // Add the input text
  pt->Draw();  // Draw the text box on the canvas
  c->SaveAs(Form("%s.png",name));
}


void plot_from_root()
{
  TFile *f = new TFile(Form("LE_1007.root"));

  TH2D *hm = (TH2D*)f->Get("h_m");
  TH2D *he = (TH2D*)f->Get("h_e");

  hm->SetTitle("");
  he->SetTitle("");

  TColor::InvertPalette();
  TCanvas *c = new TCanvas("c","",800,600);
  plot2D(c, hm, "LEm", "#nu_{#mu}");
  plot2D(c, he, "LEe", "#nu_{e}");
  c->Close();
  f->Close();
}
