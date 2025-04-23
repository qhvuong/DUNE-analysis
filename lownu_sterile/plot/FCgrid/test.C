#include <iostream>
#include <TMath.h>
#include "DUNEStyle.h"

// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}


void test() {
    /*
    // Degrees of freedom
    const int dof = 4;

    // Define the cumulative probabilities for 1, 2, and 3 sigma
    const double prob1sigma = 0.6827;  // 1 sigma (68.27%)
    const double prob2sigma = 0.9000;  // 2 sigma (95.45%)
    const double prob3sigma = 0.9973;  // 3 sigma (99.73%)

    // Calculate the corresponding chi-squared values for 1, 2, 3 sigma
    double chi2_1sigma = TMath::ChisquareQuantile(prob1sigma, dof);
    double chi2_2sigma = TMath::ChisquareQuantile(prob2sigma, dof);
    double chi2_3sigma = TMath::ChisquareQuantile(prob3sigma, dof);

    // Print the results
    std::cout << "For " << dof << " degrees of freedom:" << std::endl;
    std::cout << "1 sigma (68.27%) corresponds to chi^2 = " << chi2_1sigma << std::endl;
    std::cout << "2 sigma (95.45%) corresponds to chi^2 = " << chi2_2sigma << std::endl;
    std::cout << "3 sigma (99.73%) corresponds to chi^2 = " << chi2_3sigma << std::endl;

    return 0;
    */



    const char out_path[] = "/exp/dune/app/users/qvuong/data/lownu/Feldman_Cousins";

    TFile *f = new TFile(Form("FCtot.root"), "READ");
    TFile *f1 = new TFile(Form("FCtot_nodet.root"), "READ");


    TH1D *h = (TH1D*)f->Get("h");
    TH1D *h1 = (TH1D*)f1->Get("h");

    h->Scale(1.0 / h->Integral()); 
    h1->Scale(1.0 / h1->Integral()); 
    h->SetMaximum(h->GetMaximum()*10.);

    auto color0 = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 1);
    auto color1 = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 2);

    h->SetLineColor(color0);
    h1->SetLineColor(color1);

    h->SetLineStyle(kSolid);
    h1->SetLineStyle(kSolid);

    h->SetLineWidth(2);
    h1->SetLineWidth(2);

    double mean = h->GetMean();
    double sigma = h->GetStdDev();

    double mean1 = h1->GetMean();
    double sigma1 = h1->GetStdDev();

    double level = mean + 3 * sigma; // 99.7% (3σ)
    double level1 = mean1 + 3 * sigma1; // 99.7% (3σ)


    TLine *line = new TLine(level, 0, level, h->GetMaximum());
    TLine *line1 = new TLine(level1, 0, level1, h->GetMaximum());

    line->SetLineColor(color0);
    line1->SetLineColor(color1);

    line->SetLineStyle(kDashed);
    line1->SetLineStyle(kDashed);

    line->SetLineWidth(2);
    line1->SetLineWidth(2);


    TCanvas *c = new TCanvas("c","",800,600);
    c->SetGrid();
    c->SetLogy();
    h->Draw();
    h1->Draw("same");
    line->Draw("same");
    line1->Draw("same");
    TLegend *leg = MakeLegend(0.6, 0.65, 0.85, 0.85);
    leg->AddEntry(h, "All uncertainties");
    leg->AddEntry(h1, "No detector uncertainty");
    leg->Draw();
    dunestyle::CenterTitles(h);
    dunestyle::WIP();
    c->SaveAs("FCchi2.png");
}
