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
  TLegend * leg = MakeLegend(0.55, 0.4, 0.85, 0.75);
  TH1 * hFirst = nullptr;
  for (std::size_t histIdx = 0; histIdx < 5; histIdx++)
  {
    TH1D *h = hists[histIdx];
    auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, histIdx==0 ? 0 : -1);
    h->SetLineColor(color);
    h->SetLineWidth(2.);
    dunestyle::CenterTitles(h);
    auto newh = h->DrawCopy(histIdx == 0 ? "hist" : "hist same");  // need to leak it so it doesn't disappear
    if (!hFirst)
      hFirst = newh;

    // we do this the hard way so the legend has the top-most histogram in the stack first
    leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(newh, Form("%s < %.1f GeV", cut_name, nu_cuts[histIdx]), "l"));
  }
  c->RedrawAxis();  // otherwise the last histogram drawn overlaps with the frame
  leg->Draw();
  hFirst->SetMaximum(hFirst->GetMaximum()*1.4); // make some space for the watermark
  //hFirst->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  //hFirst->GetYaxis()->SetTitle("Lepton Selection Efficiency");
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



void plot2D(TCanvas *c, TH2D* hists, const char *name, int count, const char *sample_name)
{
  double nu_cuts[5] = {50., 10., 1., 0.5, 0.3};
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

  TPaveText *pt = new TPaveText(0.65, 0.68, 0.75, 0.88, "NDC");
  pt->SetFillStyle(0);
  pt->SetFillColor(0);  // Transparent background
  pt->SetTextColor(kBlack);
  pt->SetTextSize(0.06);
  pt->SetBorderSize(0); // No border
  pt->AddText(Form("#nu < %.1f GeV", nu_cuts[count]));  // Add the input text
  pt->AddText(sample_name);
  pt->Draw();  // Draw the text box on the canvas


  c->SaveAs(Form("%s.png",name));
}




void plot_eff()
{

  TFile *f = new TFile("outFile_forPlots.root", "READ");
  //f->Print();

  TH1D *hmEv[5], *heEv[5];
  TH1D *hm[5], *he[5];
  TH1D *hmEff[5], *heEff[5];
  TH1D *hmResE[5], *heResE[5];


  TH2D *hmEvRecoVsEv[5], *heEvRecoVsEv[5];
  TH2D *hThetaVsEe[5], *hThetaVsEeReco[5];
  TH2D *hThetaVsEe_z[5], *hThetaVsEeReco_z[5];


  for (int histIdx = 0; histIdx < 5; histIdx++) {
    hm[histIdx] = (TH1D*)f->Get(Form("hm%d", histIdx));
    he[histIdx] = (TH1D*)f->Get(Form("he%d", histIdx));

    hmEff[histIdx] = (TH1D*)hm[histIdx]->Clone();
    heEff[histIdx] = (TH1D*)he[histIdx]->Clone();

    hmEff[histIdx]->SetName(Form("hmEff%d", histIdx));
    heEff[histIdx]->SetName(Form("heEff%d", histIdx));

    hmEv[histIdx] = (TH1D*)f->Get(Form("hmEv%d", histIdx));
    heEv[histIdx] = (TH1D*)f->Get(Form("heEv%d", histIdx));

    hmEff[histIdx]->Divide(hmEv[histIdx]);
    heEff[histIdx]->Divide(heEv[histIdx]);

    hmResE[histIdx] = (TH1D*)f->Get(Form("hmResE%d", histIdx));
    heResE[histIdx] = (TH1D*)f->Get(Form("heResE_wc%d", histIdx));

    hmResE[histIdx]->GetXaxis()->SetTitle("(E_{#nu}^{reco} - E_{#nu}^{true})/E_{#nu}^{true}");
    hmResE[histIdx]->GetYaxis()->SetTitle("Events / 10^{21} POT / 0.32 GeV");

    heResE[histIdx]->GetXaxis()->SetTitle("(E_{#nu}^{reco} - E_{#nu}^{true})/E_{#nu}^{true}");
    heResE[histIdx]->GetYaxis()->SetTitle("Events / 10^{21} POT / 0.32 GeV");

    hmEvRecoVsEv[histIdx] = (TH2D*)f->Get(Form("hmEvRecoVsEv%d", histIdx));
    heEvRecoVsEv[histIdx] = (TH2D*)f->Get(Form("heEvRecoVsEv%d", histIdx));

    hThetaVsEe[histIdx] = (TH2D*)f->Get(Form("hThetaVsEe%d", histIdx));
    hThetaVsEeReco[histIdx] = (TH2D*)f->Get(Form("hThetaVsEeReco%d", histIdx));
    hThetaVsEe_z[histIdx] = (TH2D*)f->Get(Form("hThetaVsEe_z%d", histIdx));
    hThetaVsEeReco_z[histIdx] = (TH2D*)f->Get(Form("hThetaVsEeReco_z%d", histIdx)); 


    hmEvRecoVsEv[histIdx]->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    hmEvRecoVsEv[histIdx]->GetYaxis()->SetTitle("E_{#nu}^{reco} (GeV)");

    heEvRecoVsEv[histIdx]->GetXaxis()->SetTitle("E_{#nu} (GeV)");
    heEvRecoVsEv[histIdx]->GetYaxis()->SetTitle("E_{#nu}^{reco} (GeV)");

    hThetaVsEe[histIdx]->GetXaxis()->SetTitle("E_{e} (GeV)");
    hThetaVsEe[histIdx]->GetYaxis()->SetTitle("#theta_{e} (mrad)");

    hThetaVsEeReco[histIdx]->GetXaxis()->SetTitle("E_{e} (GeV)");
    hThetaVsEeReco[histIdx]->GetYaxis()->SetTitle("#theta_{e}^{reco} (mrad)");

    hThetaVsEe_z[histIdx]->GetXaxis()->SetTitle("E_{e} (GeV)");
    hThetaVsEe_z[histIdx]->GetYaxis()->SetTitle("#theta_{e} (mrad)");

    hThetaVsEeReco_z[histIdx]->GetXaxis()->SetTitle("E_{e} (GeV)");
    hThetaVsEeReco_z[histIdx]->GetYaxis()->SetTitle("#theta_{e}^{reco} (mrad)");

  }


  TColor::InvertPalette();
  TCanvas *c = new TCanvas("c","",800,600);


  plot1D(c, hmResE, "true #nu", "mResE", "#nu_{#mu}-CC");
  plot1D(c, heResE, "true #nu", "eResE", "#nu_{e}-CC");


/*
  for(int i=0; i<5; i++){
    plot2D(c, hmEvRecoVsEv[i], Form("mEvRecoVsEv%d",i), i, "#nu_{#mu}-CC");
    plot2D(c, heEvRecoVsEv[i], Form("eEvRecoVsEv%d",i), i, "#nu_{e}-CC");
    plot2D(c, hThetaVsEe[i], Form("ThetaVsEe%d",i), i, "#nu_{e}-CC");
    plot2D(c, hThetaVsEeReco[i], Form("ThetaVsEeReco%d",i), i, "#nu_{e}-CC");
    plot2D(c, hThetaVsEe_z[i], Form("ThetaVsEe_z%d",i), i, "#nu_{e}-CC");
    plot2D(c, hThetaVsEeReco_z[i], Form("ThetaVsEeReco_z%d",i), i, "#nu_{e}-CC");
  }


  hm[0]->GetYaxis()->SetTitle("Lepton Selection Efficiency");
  he[0]->GetYaxis()->SetTitle("Lepton Selection Efficiency");
  plot1D(c, hm, "true #nu", "m", "#nu_{#mu}-CC");
  plot1D(c, he, "true #nu", "e", "#nu_{e}-CC");
  c->Close();



  TFile *out = new TFile("eff.root","RECREATE");
  for(int i=0; i<5; i++){
    hm[i]->Write();
    he[i]->Write();
    hmEff[i]->Write();
    heEff[i]->Write();
  }
  out->Close();
*/
}




