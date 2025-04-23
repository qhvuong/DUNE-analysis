#include "DUNEStyle.h"



void plotContours(TCanvas *c, TH2D **h0, TH2D **h1, const char *name, double contourLevel)
{
    c->Clear();
    c->cd();
    c->SetLogx();
    c->SetLogy();
    c->SetLogz();

    TLegend *leg = MakeLegend(0.7, 0.45, 0.85, 0.8);

    const int N_HISTS = 5;  // Number of histograms per set
    double dm[N_HISTS] = {0.5, 1.0, 5.0, 10.0, 100.0};
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
    gPad->DrawFrame(xMin, yMin, xMax, yMax, "3#sigma contours;U_{#mu4}^{2};U_{e4}^{2}");

    std::vector<int> linestyles = {kSolid, kDashed};

    for (int i = 0; i < N_HISTS; i++) {
        auto color = dunestyle::colors::NextColor();

        for (int setIdx = 0; setIdx < 2; setIdx++) {
            TH2D* h2d = (setIdx == 0) ? histSet1[i] : histSet2[i];

            // Directly use sigma^2 level
            double sigmaSquaredLevel = contourLevel * contourLevel;  // e.g., 3 * 3 = 9 for 3σ
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

            leg->AddEntry(graphs.back(), Form("U_{#tau4}^{2} = %.1f, #Deltam^{2} = %.1f", Ut42[setIdx], dm[i]), "l");
        }
    }

    leg->Draw();
    //dunestyle::CenterTitles(histSet1[0]);
    dunestyle::WIP();
    c->SaveAs(Form("%s.png", name));
}





void plot_0()
{
  int N_tot = 10000;
  int N_samples = 5;

  double ue42, um42, ut42, dm2, bf_ue42, bf_um42, bf_ut42, bf_dm2, chi2, nochi2, dchi2;
  const char out_path[] = "/pnfs/dune/scratch/users/qvuong/output/contour_grid";
  double dm[5] = {0.5, 1, 5, 10, 100};
  

  double U_max = 0.7;
  double U_min = 0.0001;

  int N = 100;
  double log_min = log10(U_min);
  double log_max = log10(U_max);
  double binWidth = (log_max - log_min) / N;

  double bin_edges[N+1], bin_center[N];
  for(int i=0; i<N+1; i++){
    bin_edges[i] = pow(10, (log_min + i * binWidth));
    //bin_edges[i+1] = U_min + pow(10, (log_min + (i+1) * binWidth));
    //std::cout << bin_edges[i] << "\t";
  }

  TH2D** h = new TH2D*[N_samples];
  for(int i=0; i<N_samples; i++){
    h[i] = new TH2D(Form("h_%d",i),"Sensitivity Contour",N,bin_edges,N,bin_edges);
    h[i]->SetStats(0);
  }

  //gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  //gStyle->SetNumberContours(999);

  //TCanvas *c = new TCanvas("c","",800,600);
  //c->SetLogx();
  //c->SetLogy();


  for(int j=0; j<N_samples; j++){
  std::cout << "dm2 = " << dm[j] << "\n";

  for(int i=0; i<N_tot; i++){
    string name = Form("%s/fixed23_0%d/output_%d.txt",out_path,j,i);
    //std::cout << name << "\n";
    ifstream f(name);
    if(i%1000==0) std::cout << i*100./N_tot << " percent...\n";
    if(!f) {
      std::cout << "failed: " << i << "\n";
      continue;
    }
    else{
      f >> ue42 >> um42 >> ut42 >> dm2 >> nochi2;
      dchi2 = sqrt(nochi2 - chi2);
      h[j]->Fill(um42, ue42, dchi2);
      f.close();
    }
  }

  //h[j]->Draw("colz");

  TFile *out = new TFile(Form("fixed23_0%d.root", j), "RECREATE");
  h[j]->Write();
  out->Close();


  }
}



/*  
  TFile *f0 = new TFile(Form("contour_dm1.root"),"READ");
  TFile *f1 = new TFile(Form("contour_dm5.root"),"READ");
  TFile *f2 = new TFile(Form("contour_dm10.root"),"READ");
  TFile *f3 = new TFile(Form("contour_dm100.root"),"READ");

  TH2D *h0 = (TH2D*)f0->Get("h");
  TH2D *h1 = (TH2D*)f1->Get("h");
  TH2D *h2 = (TH2D*)f2->Get("h");
  TH2D *h3 = (TH2D*)f3->Get("h");

  double contours[1];
  contours[0] = 3; 

  h0->SetContour(1, contours);
  h1->SetContour(1, contours);
  h2->SetContour(1, contours);
  h3->SetContour(1, contours);

  h0->SetStats(0);
  h1->SetStats(0);
  h2->SetStats(0);
  h3->SetStats(0);

  h0->GetXaxis()->SetRangeUser(1e-4,0.7);
  h0->GetYaxis()->SetRangeUser(1e-4,0.7);
  h1->GetXaxis()->SetRangeUser(1e-4,0.7);
  h1->GetYaxis()->SetRangeUser(1e-4,0.7);
  h2->GetXaxis()->SetRangeUser(1e-4,0.7);
  h2->GetYaxis()->SetRangeUser(1e-4,0.7);
  h3->GetXaxis()->SetRangeUser(1e-4,0.7);
  h3->GetYaxis()->SetRangeUser(1e-4,0.7);
 
  h0->SetLineColor(kBlack);
  h0->SetLineStyle(2);
  h0->SetLineWidth(2);
  h1->SetLineColor(kRed);
  h1->SetLineStyle(1);
  h1->SetLineWidth(2);
  h2->SetLineColor(kBlue);
  h2->SetLineStyle(9);
  h2->SetLineWidth(2);
  h3->SetLineColor(kGreen);
  //h3->SetLineStyle(9);
  h3->SetLineWidth(2);

  h0->SetTitle("3#sigma Contours");
  h0->GetXaxis()->SetTitle("U_{#mu4}^{2}");
  h0->GetYaxis()->SetTitle("U_{e4}^{2}");

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  gStyle->SetNumberContours(999);

  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  c->SetLogx();
  c->SetLogy();
  c->SetLogz();
  h0->Draw("cont3 same");
  h1->Draw("cont3 same");
  h2->Draw("cont3 same");
  h3->Draw("cont3 same");
  TLegend *lg = new TLegend(0.65,0.75,0.9,0.9);
  lg->AddEntry(h0,"#Deltam^{2} = 1.0");
  lg->AddEntry(h1,"#Deltam^{2} = 5.0");
  lg->AddEntry(h2,"#Deltam^{2} = 10.0");
  lg->AddEntry(h3,"#Deltam^{2} = 100.0");
  lg->Draw();
  c->SaveAs(Form("3sigma.png"));

}

*/
