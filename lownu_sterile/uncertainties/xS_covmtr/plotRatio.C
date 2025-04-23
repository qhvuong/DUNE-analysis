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

void plot1D_singleDial(TCanvas *c, TH1D *hn, TH1D *ho, TH1D *hp, const char *cut_name, const char *out_name, const char *sample_name) {
  c->Clear();
  auto color1 = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 1);
  auto color2 = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 2);
  auto color3 = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 3);
  hn->SetLineColor(color1);
  ho->SetLineColor(kBlack);
  hp->SetLineColor(color2);

  hn->SetLineWidth(2.);
  ho->SetLineWidth(2.);
  hp->SetLineWidth(2.);

  ho->SetTitle(Form("%s #pm1#sigma", cut_name));

  // Create two pads
  TPad *pad1 = new TPad("pad1", "Upper Pad", 0, 0.25, 1, 1.0);  // Upper pad for histograms
  TPad *pad2 = new TPad("pad2", "Lower Pad", 0, 0, 1, 0.3);    // Lower pad for ratio
  pad1->SetBottomMargin(0.02);  // Reduce bottom margin for upper pad
  pad2->SetTopMargin(0.02);     // Reduce top margin for lower pad
  pad2->SetBottomMargin(0.3);   // Increase bottom margin for axis labels in lower pad
  pad1->Draw();
  pad2->Draw();

  // Draw histograms on the upper pad
  pad1->cd();
  //pad1->SetLogy();
  TLegend *leg = MakeLegend(0.6, 0.45, 0.85, 0.75);
  ho->Draw("hist");
  hn->Draw("hist same");
  hp->Draw("hist same");

  leg->AddEntry(hp, "+1 #sigma");
  leg->AddEntry(ho, "nominal");
  leg->AddEntry(hn, "-1 #sigma");
  c->RedrawAxis();
  leg->Draw();
  ho->SetMaximum(ho->GetMaximum() * 1.5);
  ho->GetXaxis()->SetTitle("E_{lep} (GeV)");
  ho->GetYaxis()->SetTitle("Entries / 10^{21} POT");
  ho->SetTitleSize(0.05, "X");
  ho->SetTitleSize(0.05, "Y");
  
  dunestyle::Simulation()->SetTextSize(0.08);
 
  TPaveText *pt = new TPaveText(0.65, 0.78, 0.75, 0.88, "NDC");
  pt->SetFillStyle(0);
  pt->SetFillColor(0);  // Transparent background
  pt->SetTextColor(kBlack);
  pt->SetTextSize(0.08);
  pt->SetBorderSize(0); // No border
  pt->AddText(sample_name);  // Add the input text
  pt->Draw();  // Draw the text box on the canvas

  // Draw ratio plot on the lower pad
  pad2->cd();
  TH1D *hnom = (TH1D*)ho->Clone();
  TH1D *hpos = (TH1D*)hp->Clone();
  TH1D *hneg = (TH1D*)hn->Clone();

  hpos->Divide(ho);
  hneg->Divide(ho);
  hnom->Divide(ho);

  hpos->SetLineWidth(2);
  hneg->SetLineWidth(2);
  hnom->SetLineWidth(2);


  hnom->SetMaximum(1.5);
  hnom->SetMinimum(0.5);
  hnom->SetTitle(0);
  hnom->GetYaxis()->SetTitle("Shifted/Nom");
  hnom->GetYaxis()->SetNdivisions(505);
  hnom->GetYaxis()->SetTitleSize(0.13);
  hnom->GetYaxis()->SetLabelSize(0.08);
  hnom->GetYaxis()->SetTitleOffset(0.4);
  hnom->GetXaxis()->SetTitle("E_{lep} (GeV)");
  hnom->GetXaxis()->SetTitleSize(0.14);
  hnom->GetXaxis()->SetLabelSize(0.08);
  hnom->GetXaxis()->SetTitleOffset(1.0);

  hnom->Draw("hist");
  hneg->Draw("hist same");
  hpos->Draw("hist same");
  c->RedrawAxis();
  c->SaveAs(Form("ratios/%s.png",out_name));
}


void plot1D_multipleDials(TCanvas *c, TH1D** hists, TH1D *hnom, const char *out_name, const char *sample_name) {
  c->Clear();
  //auto colorNOM = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 6);
  hnom->SetLineColor(kBlack);
  hnom->SetLineWidth(2.);

  // Create two pads
  TPad *pad1 = new TPad("pad1", "Upper Pad", 0, 0.25, 1, 1.0);  // Upper pad for histograms
  TPad *pad2 = new TPad("pad2", "Lower Pad", 0, 0, 1, 0.3);    // Lower pad for ratio
  pad1->SetBottomMargin(0.02);  // Reduce bottom margin for upper pad
  pad2->SetTopMargin(0.02);     // Reduce top margin for lower pad
  pad2->SetBottomMargin(0.3);   // Increase bottom margin for axis labels in lower pad
  pad1->Draw();
  pad2->Draw();

  // Draw histograms on the upper pad
  pad1->cd();
  std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"};
  const char *name[] = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"}; 
  int N_wgt = namelist.size();

  TLegend *leg = MakeLegend(0.45, 0.4, 0.85, 0.75);
  TH1 *hFirst = nullptr;
  for (std::size_t histIdx = 0; histIdx < 5; histIdx++)
  {
    TH1D *h = hists[histIdx];
    auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, histIdx == 0 ? 1 : -1);
    h->SetLineColor(color);
    h->SetLineWidth(2.);
    dunestyle::CenterTitles(h);
    auto newh = h->DrawCopy(histIdx == 0 ? "hist" : "hist same");
    if (!hFirst)
      hFirst = newh;

    leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(newh, Form("%s", name[histIdx]), "l"));
  }
  hnom->Draw("hist same");
  pad1->RedrawAxis();
  leg->AddEntry(hnom, "nominal");
  leg->Draw();
  hFirst->SetMaximum(hFirst->GetMaximum() * 1.5);
  hFirst->GetXaxis()->SetTitle("E_{lep} (GeV)");
  hFirst->GetYaxis()->SetTitle("Entries / 10^{21} POT");
  hFirst->SetTitleSize(0.05, "X");
  hFirst->SetTitleSize(0.05, "Y");
  
  dunestyle::Simulation();

  TPaveText *pt = new TPaveText(0.65, 0.78, 0.75, 0.88, "NDC");
  pt->SetFillStyle(0);
  pt->SetFillColor(0);
  pt->SetTextColor(kBlack);
  pt->SetTextSize(0.08);
  pt->SetBorderSize(0);
  pt->AddText(sample_name);
  pt->Draw();

  // Draw ratio plot on the lower pad
  pad2->cd();
  TH1D *hNominal = (TH1D*)hnom->Clone();
  for (int histIdx = 0; histIdx < 5; histIdx++)
  {
    TH1D *hShifted = hists[histIdx];
    TH1D *hRatio = (TH1D*)hShifted->Clone(Form("hRatio%d", histIdx));
    hRatio->Divide(hnom);
    hRatio->SetLineColor(hShifted->GetLineColor());
    hRatio->SetLineWidth(2);
    hRatio->SetMaximum(1.4);
    hRatio->SetMinimum(0.7);
    hRatio->GetYaxis()->SetTitle("Shifted/Nom");
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.13);
    hRatio->GetYaxis()->SetLabelSize(0.08);
    hRatio->GetYaxis()->SetTitleOffset(0.4);
    hRatio->GetXaxis()->SetTitle("E_{lep} (GeV)");
    hRatio->GetXaxis()->SetTitleSize(0.14);
    hRatio->GetXaxis()->SetLabelSize(0.08);
    hRatio->GetXaxis()->SetTitleOffset(1.0);

    hRatio->Draw(histIdx == 0 ? "hist" : "hist same");
  }
  hNominal->Divide(hnom);
  hNominal->SetLineColor(hnom->GetLineColor());
  hNominal->SetLineWidth(2.);
  hNominal->Draw("hist same");

  c->SaveAs(Form("ratios/%s.png", out_name));
}


void plotRatio()
{

  TFile *f     = new TFile("/exp/dune/app/users/qvuong/data/lownu/input_dfiles/CCtest_output.root");
  int cutNu = 4;

  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();

  std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};
  const char *name[] = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"}; 

  //std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"};
  //const char *name[] = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"}; 

  int N_wgt = namelist.size();

  std::cout << N_wgt << "\n";

  TH1D *ratio[5];
  TH1D *hAVGn[5], *hAVGp[5];
  int Nbins = 56;
  //TH1D *hAVG = new TH1D("hAVG","",Nbins,0,Nbins);

  for(int iw=0; iw<N_wgt; iw++){

  TH1D *mn       = (TH1D*)f->Get(Form("%s_mElep%d_n",name[iw],cutNu));
  TH1D *mo       = (TH1D*)f->Get(Form("%s_mElep%d_o",name[iw],cutNu));
  TH1D *mp       = (TH1D*)f->Get(Form("%s_mElep%d_p",name[iw],cutNu));

  TH1D *en       = (TH1D*)f->Get(Form("%s_eElep%d_n",name[iw],cutNu));
  TH1D *eo       = (TH1D*)f->Get(Form("%s_eElep%d_o",name[iw],cutNu));
  TH1D *ep       = (TH1D*)f->Get(Form("%s_eElep%d_p",name[iw],cutNu));

  TH1D *hmn = (TH1D*)mn->Clone();
  TH1D *hmo = (TH1D*)mo->Clone();
  TH1D *hmp = (TH1D*)mp->Clone();

  TH1D *hen = (TH1D*)en->Clone();
  TH1D *heo = (TH1D*)eo->Clone();
  TH1D *hep = (TH1D*)ep->Clone();

  plot1D_singleDial(c, hmn, hmo, hmp, name[iw], Form("m_%s",name[iw]), "#nu_{#mu}-CC");
  plot1D_singleDial(c, hen, heo, hep, name[iw], Form("e_%s",name[iw]), "#nu_{e}-CC");

  //hAVGp[iw] = (TH1D*)hp->Clone();
  //hAVGn[iw] = (TH1D*)hn->Clone();

  }
/*
  TH1D *hnom = (TH1D*)f->Get(Form("%s_eElep%d_o",name[0],cutNu));
  plot1D_multipleDials(c, hAVGp, hnom, "elRatiop", "#nu_{e}-CC");
  plot1D_multipleDials(c, hAVGn, hnom, "elRation", "#nu_{e}-CC");
*/
}
