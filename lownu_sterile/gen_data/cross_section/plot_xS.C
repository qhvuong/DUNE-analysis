#include <iostream>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TTree.h>
#include <TRandom3.h>
#include <THStack.h>
#include <TChain.h>
#include <TLegend.h>
#include <list>
#include <TPaveText.h>
#include "DUNEStyle.h"
#include <TF1.h>


// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}


void plot1D(TCanvas *c, TH1D** hists, const char *cut_name, const char *out_name, const char *sample_name)
{
  c->Clear();
  c->cd();
  double nu_cuts[5] = {50., 10., 1., 0.5, 0.3};
  TLegend * leg = MakeLegend(0.2, 0.6, 0.85, 0.78);
  leg->SetNColumns(2);
  TH1 * hFirst = nullptr;
  for (std::size_t histIdx = 0; histIdx < 5; histIdx++)
  {
    TH1D *h = hists[histIdx];
    auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, histIdx==0 ? 0 : -1);
    h->SetLineColor(color);
    h->SetLineWidth(2.);
    dunestyle::CenterTitles(h);
    //h->SetMaximum(1.);
    auto newh = h->DrawCopy(histIdx == 0 ? "hist" : "hist same");  // need to leak it so it doesn't disappear
    if (!hFirst)
      hFirst = newh;

    // we do this the hard way so the legend has the top-most histogram in the stack first
    leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(newh, Form("%s < %.1f GeV", cut_name, nu_cuts[histIdx]), "l"));
  }
  c->RedrawAxis();  // otherwise the last histogram drawn overlaps with the frame
  leg->Draw();
  hFirst->SetMaximum(100.);
  hFirst->SetMinimum(0.2);
  //hFirst->SetMaximum(hFirst->GetMaximum()*1.4); // make some space for the watermark
  hFirst->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hFirst->GetYaxis()->SetTitle("cross section (#times 10^{-38} cm^{2}/nucleon)");
  hFirst->SetTitleSize(0.05, "X");
  hFirst->SetTitleSize(0.05, "Y");
  dunestyle::Simulation();
 
  TPaveText *pt = new TPaveText(0.65, 0.78, 0.75, 0.88, "NDC");
  pt->SetFillStyle(0);
  pt->SetFillColor(0);  // Transparent background
  pt->SetTextColor(kBlack);
  pt->SetTextSize(0.08);
  pt->SetBorderSize(0); // No border
  pt->AddText(sample_name);  // Add the input text
  pt->Draw();  // Draw the text box on the canvas

  c->SaveAs(Form("%s.png",out_name));
}


void plot_xS()
{

  TFile *f = new TFile("xS.root", "READ");
  TH1D *xSm[5], *xSe[5];

  for (int i = 0; i < 5; i++) {
    xSm[i] = (TH1D*)f->Get(Form("hm%d", i));
    xSe[i] = (TH1D*)f->Get(Form("he%d", i));
   
  }

  //hm[0]->GetYaxis()->SetTitle("Lepton Selection Efficiency");
  //he[0]->GetYaxis()->SetTitle("Lepton Selection Efficiency");
  TCanvas *c = new TCanvas("c","",800,600);
  c->SetLogy();
  plot1D(c, xSm, "true #nu", "xSm", "#nu_{#mu}-CC");
  plot1D(c, xSe, "true #nu", "xSe", "#nu_{e}-CC");
  //c->Close();

}




