#include "DUNEStyle.h"
#include <TLegend.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TPolyLine.h>

// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(
  return leg;
}

void draw90ContoursPads_mue(TCanvas *c, TH2D* hists[6], const char *out_name) {
    c->Clear();

    const int N_HISTS = 6;  // Number of histograms per canvas
    double dm[N_HISTS] = {0.1, 0.5, 1.0, 5.0, 9.9, 90.2};
    const char* expNAMES[4] = {"MICROBOONE", "STEREO", "PROSPECT", "MINOS"};
    const char* sensNAMES[4] = {"MICROBOONE 95% CL", "STEREO 95% CL", "PROSPECT 95% CL", "MINOS+ 90% CL"};
    double sens[4][N_HISTS]= {
                {0.06683198783377048, 0.012719085976607295, 0.004651022504075903, 0.0019230246120471976, 0.001954493050767267, 0.0023150843640797674}, //4Ue42Um42 95%
                {0.8823915399615171, 0.2994825101247634, 0.0701635896779922, 0.09406476191238712, 0.7923889449166047, 10.},  //Ue42 95%
                {10., 0.3122166697660653, 0.0763768816631674, 0.06905369387506037, 0.14742195161523622, 10.}, //Ue42 95%
                {0.0070997237900409845, 0.0062363746782495766, 0.006184341858927206, 0.015343187563158006, 0.009001630974621208, 0.019043389042955775} //Um42 90%
    };

    // Set up the drawing area manually
    double xMin = hists[0]->GetXaxis()->GetXmin();
    double xMax = hists[0]->GetXaxis()->GetXmax();
    double yMin = hists[0]->GetYaxis()->GetXmin();
    double yMax = hists[0]->GetYaxis()->GetXmax();

    const char* xtitle = hists[0]->GetXaxis()->GetTitle();
    const char* ytitle = hists[0]->GetYaxis()->GetTitle();

    const int N_POINTS = 100;

    // Confidence levels (sigma^2 values for contours)
    double sigmaLevels = 9.33333;  // 1σ, 2σ, 3σ levels

    c->Divide(3, 2, 0.001, 0.001); 

    for (int i = 0; i < N_HISTS; ++i) {
        c->cd(i+1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gPad->SetTopMargin(0.12);     // Shrink top margin
        gPad->SetBottomMargin(0.15);  // Shrink bottom margin
        gPad->SetLeftMargin(0.15);    // Shrink left margin
        gPad->SetRightMargin(0.15);   // Shrink right margin
        gPad->SetTicks();
        TLegend *leg = MakeLegend(0.18, 0.18, 0.55, 0.35);

        gPad->DrawFrame(xMin, yMin, xMax, yMax, Form("%s;%s;%s", hists[i]->GetTitle(), xtitle, ytitle));

        std::vector<TGraph*> graphDUNE = dunestyle::GetContourGraphs(hists[i], sigmaLevels);
        for (TGraph* g : graphDUNE) {
            g->SetLineColor(kBlack);
            //g->SetLineStyle(linestyles[setIdx]);
            g->SetLineWidth(2);
            g->Draw("L same");
        }
        leg->AddEntry(graphDUNE.back(), Form("DUNE ND 90%% CL"), "l");

        Color_t color[4];
        for(int icolor=0; icolor<4; icolor++){
            color[icolor] = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, icolor==0 ? 1 : -1);
        }

        for (int expIdx = 0; expIdx < 4; expIdx++) 
        {
            auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, expIdx==0 ? 1 : -1);
            if(sens[expIdx][i] == 10.) continue;
            TGraph* g = new TGraph();
            g->SetLineColor(color);
            g->SetLineWidth(2);
            for (int j = 0; j < N_POINTS; ++j) {
                double x = xMin * pow(xMax / xMin, double(j) / (N_POINTS - 1));
                double y = yMin * pow(yMax / yMin, double(j) / (N_POINTS - 1));

                if (strcmp(expNAMES[expIdx], "MICROBOONE") == 0) y = sens[expIdx][i] / (4*x);
                else if (strcmp(expNAMES[expIdx], "MINOS") == 0) x = sens[expIdx][i];
                else y = ( 1-sqrt(1-sens[expIdx][i]) ) / 2.;
                g->SetPoint(j, x, y);
            }
            //funcGraphs.push_back(g);
            g->Draw("L same");
            leg->AddEntry(g, Form("%s", sensNAMES[expIdx]), "l");
        }
        leg->Draw();
        dunestyle::WIP()->SetTextSize(0.06);
    }

    c->SaveAs(Form("%s.png",out_name));
}

void drawSix2DHistos(TCanvas *c, TH2D* hists[6], const char *out_name) {
    c->Clear();

    double globalMin = hists[0]->GetMinimum();
    double globalMax = hists[0]->GetMaximum();
    for (int i = 1; i < 6; ++i) {
        globalMin = std::min(globalMin, hists[i]->GetMinimum());
        globalMax = std::max(globalMax, hists[i]->GetMaximum());
    }
    std::cout << globalMin << "\t" << globalMax << "\n";

    c->Divide(3, 2, 0.001, 0.001);

    for (int i = 0; i < 6; ++i) {
        c->cd(i+1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gPad->SetTopMargin(0.12);     // Shrink top margin
        gPad->SetBottomMargin(0.15);  // Shrink bottom margin
        gPad->SetLeftMargin(0.15);    // Shrink left margin
        gPad->SetRightMargin(0.15);   // Shrink right margin
        gPad->SetTicks();
        hists[i]->SetMinimum(0.001);
        hists[i]->SetMaximum(2E3);
        hists[i]->Draw("colz");  // no 'z' to hide individual color bars
    }

    c->SaveAs(Form("%s.png",out_name));
}

void drawContoursPads(TCanvas *c, TH2D* hists[6], const char *out_name) {
    c->Clear();

    // Set up the drawing area manually
    double xMin = hists[0]->GetXaxis()->GetXmin();
    double xMax = hists[0]->GetXaxis()->GetXmax();
    double yMin = hists[0]->GetYaxis()->GetXmin();
    double yMax = hists[0]->GetYaxis()->GetXmax();

    const char* xtitle = hists[0]->GetXaxis()->GetTitle();
    const char* ytitle = hists[0]->GetYaxis()->GetTitle();

    // Confidence levels (sigma^2 values for contours)
    std::vector<double> sigmaLevels = {5.6, 12.5333, 38.6667};  // 1σ, 2σ, 3σ levels

    c->Divide(3, 2, 0.001, 0.001); 

    for (int i = 0; i < 6; ++i) {
        c->cd(i+1);
        gPad->SetLogx();
        gPad->SetLogy();
        gPad->SetLogz();
        gPad->SetTopMargin(0.12);     // Shrink top margin
        gPad->SetBottomMargin(0.15);  // Shrink bottom margin
        gPad->SetLeftMargin(0.15);    // Shrink left margin
        gPad->SetRightMargin(0.15);   // Shrink right margin
        gPad->SetTicks();
        TLegend *leg = MakeLegend(0.18, 0.18, 0.55, 0.35);

        gPad->DrawFrame(xMin, yMin, xMax, yMax, Form("%s;%s;%s", hists[i]->GetTitle(), xtitle, ytitle));

        for (size_t sigmaIdx = 0; sigmaIdx < sigmaLevels.size(); sigmaIdx++) 
        {
            double sigmaSquaredLevel = sigmaLevels[sigmaIdx];

            auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, sigmaIdx==0 ? 0 : -1);

            // Process the first histogram (h0)
            std::vector<TGraph*> graphs = dunestyle::GetContourGraphs(hists[i], sigmaSquaredLevel);
            if (!graphs.empty()) {
                for (TGraph *g : graphs) {
                    g->SetLineColor(color);
                    //g->SetLineStyle(linestyles[sigmaIdx]);
                    g->SetLineWidth(2);
                    g->Draw("same");
                }
                leg->AddEntry(graphs.back(), Form("%zu#sigma", sigmaIdx+1), "l");
            }
        }
        leg->Draw();
        dunestyle::WIP()->SetTextSize(0.06);
    }

    c->SaveAs(Form("%s.png",out_name));
}



void plot_from_root()
{
  const int N_samples = 6;
  double U_vals[N_samples] = {0.00104563, 0.0104497, 0.0194204, 0.0514304, 0.104431, 0.194082};
  double M_vals[N_samples] = {0.103576, 0.507293, 1.01218, 4.95746, 9.89143, 90.2109};

  TH2D** hMe = new TH2D*[N_samples];
  TH2D** hMm = new TH2D*[N_samples];
  TH2D** hme = new TH2D*[N_samples];

  TFile *fMe = new TFile(Form("contours_um42.root"),"READ");
  TFile *fMm = new TFile(Form("contours_ue42.root"),"READ");
  TFile *fme = new TFile(Form("contours_dm2.root"),"READ");

  for(int i=0; i<N_samples; i++){
    hMe[i] = (TH2D*)fMe->Get(Form("h_Me_%d",i));
    hMm[i] = (TH2D*)fMm->Get(Form("h_Mm_%d",i));
    hme[i] = (TH2D*)fme->Get(Form("h_me_%d",i));

    hme[i]->SetTitle(Form("#Deltam^{2} = %.1f eV^{2}", M_vals[i]));
    hme[i]->GetXaxis()->SetTitle("U_{#mu4}^{2}");
    hme[i]->GetYaxis()->SetTitle("U_{e4}^{2}");

    hMe[i]->SetTitle(Form("U_{#mu4}^{2} = %.3f", U_vals[i]));
    hMe[i]->GetXaxis()->SetTitle("U_{e4}^{2}");
    hMe[i]->GetYaxis()->SetTitle("#Deltam^{2}");

    hMm[i]->SetTitle(Form("U_{e4}^{2} = %.3f", U_vals[i]));
    hMm[i]->GetXaxis()->SetTitle("U_{#mu4}^{2}");
    hMm[i]->GetYaxis()->SetTitle("#Deltam^{2}");
  }

  TColor::InvertPalette();
  TCanvas *c = new TCanvas("c","",1400,800);

/*
  drawSix2DHistos(c, hMe, "contour_Me");
  drawSix2DHistos(c, hme, "contour_mue");
  drawSix2DHistos(c, hMm, "contour_Mmu");

  drawContoursPads(c, hMe, "sigmaContours_Me");
  drawContoursPads(c, hme, "sigmaContours_mue");
  drawContoursPads(c, hMm, "sigmaContours_Mmu");
*/

  draw90ContoursPads(c, hme, "90CL_mue");


}


