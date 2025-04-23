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



void overlayFunctionOnCanvas(double A) {
    // Assume the canvas with histograms is already drawn
    TCanvas *c = gPad->GetCanvas();  // Get current active canvas

    if (!c) {
        std::cerr << "No active canvas found!" << std::endl;
        return;
    }

    c->cd();  // Select the current canvas

    // Get X and Y axis limits from the existing histograms
    double xMin = gPad->GetUxmin();
    double xMax = gPad->GetUxmax();
    double yMin = gPad->GetUymin();
    double yMax = gPad->GetUymax();

    std::cout << "X Range: [" << xMin << ", " << xMax << "], Y Range: [" << yMin << ", " << yMax << "]\n";

    // Define number of points for smooth curve
    const int N = 100;  
    TGraph *graph = new TGraph();

    for (int i = 0; i < N; i++) {
        // Generate log-spaced x values
        double x = xMin * pow(xMax / xMin, double(i) / (N - 1));
        double y = 4 / (A * x);  // Compute y

        if (y >= yMin && y <= yMax) {  // Ensure y is within the visible range
            graph->SetPoint(graph->GetN(), x, y);
        }
    }

    // Set graph style
    graph->SetLineColor(kRed);
    graph->SetLineWidth(2);
    graph->SetLineStyle(2);  // Dashed line

    // Draw the function on top of existing histograms
    graph->Draw("L SAME");  // "L" for line, "SAME" to overlay

    // Refresh canvas
    //c->Update();
}




void plotContours(TCanvas *c, TH2D **h0, TH2D **h1, const char *name, double contourLevel)
{
    c->Clear();
    c->cd();
    c->SetLogx();
    c->SetLogy();
    c->SetLogz();

    const int N_POINTS = 100;

    TLegend *leg = MakeLegend(0.16, 0.18, 0.5, 0.65);

    const int N_HISTS = 5;  // Number of histograms per set
    double dm[N_HISTS] = {0.5, 1.0, 5.0, 10.0, 100.0};
    //double MICROBOONE_green[N_HISTS] = {0.0027972328383502575, 0.0009949978637502671, 0.000804496955881946, 0.0006832494538781662, 0.0007367531158285329};
    //double MICROBOONE_magenta[N_HISTS] = {0.0031797923608898947, 0.0012834157557737058, 0.0010032923066174695, 0.0011964803957119961, 0.0012558832114312577};
    double MICROBOONE[N_HISTS] = {0.01255542751665414, 0.004664660387808786, 0.001858961240723379, 0.00206226787648898, 0.0023441512391716925};

    double STEREO[N_HISTS] = {0.2713712113607749, 0.0653907064421432, 0.0826638834911441, 0.7342370605769079, 1.};
    
    double MINOS[N_HISTS] = {0.00616, 0.00609, 0.01445, 0.00888, 0.01915};
    double Ut42[2] = {0.0, 0.5};
    std::vector<TH2D*> histSet1, histSet2;

    for (int i = 0; i < N_HISTS; i++) {
        histSet1.push_back(h0[i]);
        histSet2.push_back(h1[i]);
    }

    // Set up the drawing area manually
    double xMin = histSet1[0]->GetXaxis()->GetXmin();
    double xMax = histSet1[0]->GetXaxis()->GetXmax();
    double yMin = histSet1[0]->GetYaxis()->GetXmin();
    double yMax = histSet1[0]->GetYaxis()->GetXmax();
    std::cout << "X Range: [" << xMin << ", " << xMax << "], Y Range: [" << yMin << ", " << yMax << "]\n";

    gPad->DrawFrame(xMin, yMin, xMax, yMax, ";U_{#mu4}^{2};U_{e4}^{2}");

    std::vector<int> linestyles = {kSolid, kDashed};

    for (int i = 0; i < N_HISTS; i++) {
        auto color = dunestyle::colors::NextColor();

        for (int setIdx = 0; setIdx < 1; setIdx++) {
            TH2D* h2d = (setIdx == 0) ? histSet1[i] : histSet2[i];

            // Directly use sigma^2 level
            double sigmaSquaredLevel = contourLevel;  // e.g., 3 * 3 = 9 for 3σ
            std::cout << "Processing Histogram " << i << ", Set " << setIdx << std::endl;
            std::cout << "Contour Level (sigma^2): " << sigmaSquaredLevel << std::endl;

            // Get contour graphs at the specified sigma^2 level
            std::vector<TGraph*> graphs = dunestyle::GetContourGraphs(h2d, sigmaSquaredLevel);
            if (graphs.empty()) {
                std::cerr << "No contours found at level " << sigmaSquaredLevel << " for histogram " << i 
                          << ", Set " << setIdx << std::endl;
                continue;
            }

            for (TGraph* g : graphs) {
                g->SetLineColor(color);
                g->SetLineStyle(linestyles[setIdx]);
                g->SetLineWidth(2);
                g->Draw("same");
            }

            leg->AddEntry(graphs.back(), Form("DUNE ND 95%% CL #Deltam^{2} = %.1f", dm[i]), "l");
        }

        TGraph *graphMICROBOONE = new TGraph();
        TGraph *graphSTEREO = new TGraph();
        TGraph *graphMINOS = new TGraph();
        for (int j = 0; j < N_POINTS; j++) {
            double x = xMin * pow(xMax / xMin, double(j) / (N_POINTS - 1));
            double ySTEREO = (1-sqrt(1-STEREO[i]))/2;
            double yMICROBOONE = MICROBOONE[i] / (4*x);
            double y = yMin * pow(yMax / yMin, double(j) / (N_POINTS - 1));

            if (ySTEREO >= yMin && ySTEREO <= yMax) graphSTEREO->SetPoint(graphSTEREO->GetN(), x, ySTEREO);
            if (yMICROBOONE >= yMin && yMICROBOONE <= yMax) graphMICROBOONE->SetPoint(graphMICROBOONE->GetN(), x, yMICROBOONE);
            graphMINOS->SetPoint(graphMINOS->GetN(), MINOS[i], y);
        }
        
        /*
        // Set graph style
        graphSTEREO->SetLineColor(color);
        graphSTEREO->SetLineWidth(2);
        graphSTEREO->SetLineStyle(10);  // Dashed line
        graphSTEREO->Draw("L SAME");
        leg->AddEntry(graphSTEREO, Form("STEREO 95%% CL #Deltam^{2} = %.1f", dm[i]), "l");

        graphMICROBOONE->SetLineColor(color);
        graphMICROBOONE->SetLineWidth(2);
        graphMICROBOONE->SetLineStyle(3);  // Dashed line
        leg->AddEntry(graphMICROBOONE, Form("MicroBooNE 95%% CL #Deltam^{2} = %.1f", dm[i]), "l");
        graphMICROBOONE->Draw("L SAME");
        */
        
        graphMINOS->SetLineColor(color);
        graphMINOS->SetLineWidth(2);
        graphMINOS->SetLineStyle(8);  // Dashed line
        graphMINOS->Draw("L SAME");
        leg->AddEntry(graphMINOS, Form("MINOS+ 90%% CL #Deltam^{2} = %.1f", dm[i]), "l");
        
    }
    leg->Draw();
    //dunestyle::CenterTitles(histSet1[0]);
    dunestyle::WIP()->SetTextSize(0.06);
    c->SaveAs(Form("%s.png", name));

    
}




void plot1()
{
  double dm[5] = {0.5, 1, 5, 10, 100};
  int N_samples = 5;

  TH2D** h0 = new TH2D*[N_samples];
  TH2D** h1 = new TH2D*[N_samples];
  TH2D** hL0 = new TH2D*[N_samples];
  TH2D** hL1 = new TH2D*[N_samples];

  for(int i=0; i<N_samples; i++){
    TFile *f0 = new TFile(Form("fixed23_0%d.root", i),"READ");
    TFile *f1 = new TFile(Form("fixed23_1%d.root", i),"READ");
    h0[i] = (TH2D*)f0->Get(Form("h_%d",i));
    //hL0[i] = (TH2D*)f0->Get(Form("h_%d",i));
    h1[i] = (TH2D*)f1->Get(Form("h_%d",i));
    //hL1[i] = (TH2D*)f1->Get(Form("h_%d",i));

    h0[i]->SetTitle(Form("3#sigma Contours (#Deltam^{2} = %.1f eV^{2})",dm[i]));
    h0[i]->GetXaxis()->SetTitle("U_{#mu4}^{2}");
    h0[i]->GetYaxis()->SetTitle("U_{e4}^{2}");

    h1[i]->SetTitle(Form("3#sigma Contours (#Deltam^{2} = %.1f eV^{2})",dm[i]));
    h1[i]->GetXaxis()->SetTitle("U_{#mu4}^{2}");
    h1[i]->GetYaxis()->SetTitle("U_{e4}^{2}");
  }

  //std::cout << h0[0]->GetNbinsX();

  TCanvas *c = new TCanvas("c","",1000,800);

  //plotContours(c, h0, h1, "95CL", 12.5333);
  plotContours(c, h0, h1, "90CL", 9.86667);


}


