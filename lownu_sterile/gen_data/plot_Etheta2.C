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
  TLegend * leg = MakeLegend(0.55, 0.23, 0.85, 0.58);
  TH1 * hFirst = nullptr;
  for (std::size_t histIdx = 0; histIdx < 5; histIdx++)
  {
    TH1D *h = hists[histIdx];
    auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, histIdx==0 ? 0 : -1);
    h->SetLineColor(color);
    h->SetLineWidth(2.);
    dunestyle::CenterTitles(h);
    auto newh = h->DrawCopy(histIdx == 0 ? "" : "same");  // need to leak it so it doesn't disappear
    if (!hFirst)
      hFirst = newh;

    // we do this the hard way so the legend has the top-most histogram in the stack first
    leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(newh, Form("%s < %.1f GeV", cut_name, nu_cuts[histIdx]), "l"));
  }
  c->RedrawAxis();  // otherwise the last histogram drawn overlaps with the frame
  leg->Draw();
  hFirst->SetMaximum(hFirst->GetMaximum()*1.4); // make some space for the watermark
  hFirst->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hFirst->GetYaxis()->SetTitle("Event Selection Efficiency");
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



void plot2D(TCanvas *c, TH2D* hists, const char *name, const char *sample_name)
{
  c->Clear();
  c->cd();
  //hists->SetTitle(Form("#nu < %.1f", nu_cuts[count]));
  dunestyle::CenterTitles(hists);
  hists->Draw("colz");
  hists->SetTitleSize(0.05, "X");
  hists->SetTitleSize(0.06, "Y");

  TF1 *funcPos = new TF1("funcPos", "sqrt(3000./x)", 0.1, 16);  // Avoid x = 0 (undefined)
  funcPos->SetLineColor(kRed);
  funcPos->SetLineWidth(3);
  funcPos->Draw("SAME");

  dunestyle::Simulation();

  TPaveText *pt = new TPaveText(0.65, 0.78, 0.75, 0.88, "NDC");
  pt->SetFillStyle(0);
  pt->SetFillColor(0);  // Transparent background
  pt->SetTextColor(kBlack);
  pt->SetTextSize(0.06);
  pt->SetBorderSize(0); // No border
  pt->AddText(sample_name);  // Add the input text
  pt->Draw();  // Draw the text box on the canvas

  TLegend *legend = MakeLegend(0.5, 0.68, 0.82, 0.78);
  legend->AddEntry(funcPos, "E_{e}#theta_{e}^{2} = 3.0 MeV", "l");
  legend->Draw();

  c->SaveAs(Form("%s.png",name));
}




void plot_Etheta2()
{

  TFile *fCC = new TFile("CC/outFile_forPlots.root", "READ");
  TFile *fnue = new TFile("nuescattering/plot_Etheta2.root", "READ");

  TH2D *hCC, *hCCz;
  TH2D *hnue;

  hCC = (TH2D*)fCC->Get("hThetaVsEe4");
  hCCz = (TH2D*)fCC->Get("hThetaVsEe_z4");

  hnue = (TH2D*)fnue->Get("hThetaVsEe");

  hCC->GetXaxis()->SetTitle("E_{e} (GeV)");
  hCC->GetYaxis()->SetTitle("#theta_{e} (mrad)");

  hCCz->GetXaxis()->SetTitle("E_{e} (GeV)");
  hCCz->GetYaxis()->SetTitle("#theta_{e} (mrad)");

  hnue->GetXaxis()->SetTitle("E_{e} (GeV)");
  hnue->GetYaxis()->SetTitle("#theta_{e} (mrad)");


  TColor::InvertPalette();
  TCanvas *c = new TCanvas("c","",800,600);

  plot2D(c, hCCz, "ThetaVsEe_CC", "#nu_{e}-CC");
  plot2D(c, hnue, "ThetaVsEe_nue", "#nu + e");

  
}




