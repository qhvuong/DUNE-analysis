#include "DUNEStyle.h"
#include <TLegend.h>
#include <TPaveText.h>

// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}






void plotContours(TCanvas *c, TH2D *h0, TH2D *h1, const char *name)
{
    c->Clear();
    c->cd();
    c->SetLogx();
    c->SetLogy();
    c->SetLogz();

    // Create a legend
    TLegend *leg = MakeLegend(0.18, 0.18, 0.55, 0.35);

    // Set up the drawing area manually
    double xMin = h0->GetXaxis()->GetXmin();
    double xMax = h0->GetXaxis()->GetXmax();
    double yMin = h0->GetYaxis()->GetXmin();
    double yMax = h0->GetYaxis()->GetXmax();
    gPad->DrawFrame(xMin, yMin, xMax, yMax, ";U_{#mu4}^{2};U_{e4}^{2}");

    // Confidence levels (sigma^2 values for contours)
    std::vector<double> sigmaLevels = {6.13333, 12.53333, 31.73333};  // 1σ, 2σ, 3σ levels
    std::vector<int> linestyles = {kSolid, kDashed, kDotted};

    for (size_t sigmaIdx = 0; sigmaIdx < sigmaLevels.size(); sigmaIdx++) 
    {
        double sigmaSquaredLevel = sigmaLevels[sigmaIdx];

        auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, sigmaIdx);

        // Process the first histogram (h0)
        std::vector<TGraph*> graphs0 = dunestyle::GetContourGraphs(h0, sigmaSquaredLevel);
        if (!graphs0.empty()) {
            for (TGraph *g : graphs0) {
                g->SetLineColor(color);
                //g->SetLineStyle(linestyles[sigmaIdx]);
                g->SetLineWidth(2);
                g->Draw("same");
            }
            leg->AddEntry(graphs0.back(), Form("%zu#sigma", sigmaIdx+1), "l");
        }
        /*
        // Process the second histogram (h1)
        std::vector<TGraph*> graphs1 = dunestyle::GetContourGraphs(h1, sigmaSquaredLevel);
        if (!graphs1.empty()) {
            for (TGraph *g : graphs1) {
                g->SetLineColor(color1);
                g->SetLineStyle(linestyles[sigmaIdx]);
                g->SetLineWidth(2);
                g->Draw("same");
            }
            leg->AddEntry(graphs1.back(), Form("No Detector Uncertainty: %zu#sigma", sigmaIdx+1), "l");
        }
        */
    }

    // Draw the legend
    leg->Draw();
    dunestyle::WIP()->SetTextSize(0.06);

    TPaveText *pt = new TPaveText(0.65, 0.7, 0.75, 0.8, "NDC");
    pt->SetFillStyle(0);
    pt->SetFillColor(0);  // Transparent background
    pt->SetTextColor(kBlack);
    pt->SetTextSize(0.05);
    pt->SetBorderSize(0); // No border
    pt->AddText("#Deltam^{2} = 5.0 eV^{2}");  // Add the input text
    pt->Draw();  // Draw the text box on the canvas

    c->SaveAs(Form("%s.png", name));
}





void plot2()
{

    TFile *f0 = new TFile(Form("fixed23_02.root"),"READ");
    TFile *f1 = new TFile(Form("fixed23_02_nodet.root"),"READ");
    TH2D* h0 = (TH2D*)f0->Get(Form("h_2"));
    TH2D* h1 = (TH2D*)f1->Get(Form("h_2"));

    h0->SetTitle(Form("3#sigma Contours (#Deltam^{2} = 5.0 eV^{2})"));
    h0->GetXaxis()->SetTitle("U_{#mu4}^{2}");
    h0->GetYaxis()->SetTitle("U_{e4}^{2}");

    h1->SetTitle(Form("3#sigma Contours (#Deltam^{2} = 5.0 eV^{2})"));
    h1->GetXaxis()->SetTitle("U_{#mu4}^{2}");
    h1->GetYaxis()->SetTitle("U_{e4}^{2}");



  

  //gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  //gStyle->SetNumberContours(999);




  TCanvas *c = new TCanvas("c","",800,800);

  //plotContours(c, hL0, "ut0");
  //plotContours(c, hL1, "ut1");
  plotContours(c, h0, h1, "contours");

  //c->Close();



  
  /*
  //h0->Draw("cont3 same");
  h2->Draw("cont3 same");
  h1->Draw("cont3 same");
  //h3->Draw("cont3 same");
  TLegend *lg = new TLegend(0.65,0.75,0.9,0.9);
  //lg->AddEntry(h0,"#Deltam^{2} = 1.0");
  lg->AddEntry(h1,"OLD #Deltam^{2} = 10.0 eV^{2}");
  lg->AddEntry(h2,"NEW #Deltam^{2} = 10.0 eV^{2}");
  //lg->AddEntry(h3,"#Deltam^{2} = 100.0");
  lg->Draw();
  c->SaveAs(Form("3sigma_new.png"));
  */

}


