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
#include <TParameter.h>
#include "DUNEStyle.h"

// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}

void plot_test()
{
  const int N_nucut = 5;
  double nucut[N_nucut] = {10., 3., 1., 0.5, 0.3};

  TFile *in = new TFile("out.root","READ");
  std::vector<TH1D> hin, hout;

  for(int j=0; j<N_nucut; j++){
    TH1D* h_in = (TH1D*)in->Get(Form("hin%d",j));
    if (h_in) hin.emplace_back(*h_in); 

    TH1D* h_out = (TH1D*)in->Get(Form("hout%d",j));
    if (h_out) hout.emplace_back(*h_out);
  }
  in->Close();
  delete in;

  auto getColor = [](int idx) {
    return dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, idx);
  };

  
  for (int j = 0; j < N_nucut; ++j) {
    // Create a canvas for each pair of histograms
    TCanvas *c = new TCanvas(Form("c%d", j), Form("Comparison %d", j), 800, 600);
    //gPad->Update();

    // Set histogram styles for better distinction
    hin[j].SetLineColor(getColor(0));
    hin[j].SetLineWidth(2);
    hout[j].SetLineColor(getColor(1));
    hout[j].SetLineWidth(2);

    hin[j].SetTitle(Form("#nu < %.1f GeV", nucut[j]));
    hin[j].GetXaxis()->SetTitle("Elep_reco (GeV)");
    hin[j].GetYaxis()->SetTitleOffset(1.25);
    hin[j].GetYaxis()->SetTitle("counts");
    dunestyle::CenterTitles(&hin[j]);

    hin[j].SetMaximum(hin[j].GetMaximum()*1.25);

    // Draw the histograms
    hin[j].DrawCopy("hist");       // Draw the first histogram
    hout[j].DrawCopy("hist same"); // Overlay the second histogram

    TLegend * leg = MakeLegend(0.4, 0.65, 0.9, 0.78);
    leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(&hin[j], Form("LepPDG + reco_numu"), "l"));
    leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(&hout[j], Form("Only reco_numu"), "l"));
    leg->SetTextSize(0.04);
    leg->Draw();
    
    dunestyle::Simulation();

    c->RedrawAxis();
    c->SaveAs(Form("io%d.png", j));
  }




/*
  auto getColor = [](int idx) {
    return dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, idx);
  };

  TCanvas *c = new TCanvas("c","",800,600);
  // stacked histogram
  
  TLegend * leg = MakeLegend(0.68, 0.45, 0.9, 0.87);
  for (std::size_t histIdx = 0; histIdx < 1; histIdx++){
    c->Clear();
    c->cd();
    TH1D& hi = hin[histIdx];
    TH1D& ho = hout[histIdx];
    hi.SetLineColor(getColor(0));
    ho.SetLineColor(getColor(0));
    dunestyle::CenterTitles(&hi);
    auto newhi = hi.DrawCopy("same");
    auto newho = ho.DrawCopy("same");
    //leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(newhi, Form("LepPDG + reco_numu"), "l"));
    //leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(newho, Form("Only reco_numu"), "l"));
    //c->RedrawAxis();
    //leg->Draw();
    dunestyle::WIP();
  }
  
*/
  /*
  TH1 * hFirst = nullptr;
  for (std::size_t histIdx = 0; histIdx < hin.size(); histIdx++)
  {
    TH1D& h = hin[histIdx];
    auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, histIdx==0 ? 0 : -1);
    h.SetLineColor(color);
    h.SetFillStyle(0);
    dunestyle::CenterTitles(&h);
    auto newh = h.DrawCopy(histIdx == 0 ? "" : "same");  // need to leak it so it doesn't disappear
    if (!hFirst)
      hFirst = newh;

    // we do this the hard way so the legend has the top-most histogram in the stack first
    leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(newh, Form("Hist #%zu", histIdx+1), "l"));
  }
  c->RedrawAxis();  // otherwise the last histogram drawn overlaps with the frame
  leg->Draw();
  hFirst->SetMaximum(hFirst->GetMaximum()*1.25); // make some space for the watermark
  dunestyle::WIP();
  */

}



