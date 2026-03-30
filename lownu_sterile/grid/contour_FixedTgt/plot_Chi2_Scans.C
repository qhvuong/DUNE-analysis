#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include "TH2D.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TColor.h"
#include "TMarker.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TFile.h"
#include "DUNEStyle.h"

#define UE42 0.04
#define UM42 0.01
#define UT42 0.2

void plot_Chi2_Scans_grid() {
  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);
  TColor::InvertPalette();
  gStyle->SetNumberContours(999);

  const std::string prefix = "";
  const std::string suffix = "_noNuE";
  // const std::string suffix = "";

  // Plotting order (xidx, yidx):
  // 1) dm2 vs ue42, 2) dm2 vs um42, 3) dm2 vs ut42, 4) um42 vs ue42
  struct Mode { std::string name; int xidx, yidx; };
  const std::vector<Mode> modes = {
    {"ue42_dm2", 0, 3}, // top-left:  y=dm2  vs x=|Ue4|^2
    {"um42_dm2", 1, 3}, // top-right: y=dm2  vs x=|Umu4|^2
    {"ut42_dm2", 2, 3}, // bottom-left: y=dm2 vs x=|Utau4|^2
    {"ue42_um42", 0, 1} // bottom-right: y=|Umu4|^2 vs x=|Ue4|^2
  };

  // Binning
  const int nBins = 100;
  const double minVals[4] = {1e-4, 1e-4, 1e-4, 0.1};
  const double maxVals[4] = {0.7,   0.7,   0.7,   100.0};
  auto logspace = [](double min, double max, int nbins, double* edges) {
    const double log_min = std::log10(min);
    const double log_max = std::log10(max);
    for (int i = 0; i <= nbins; ++i)
      edges[i] = std::pow(10.0, log_min + i * (log_max - log_min) / nbins);
  };

  auto labelOf = [](int id)->std::string {
    if (id == 0) return "|U_{e4}|^{2}";
    if (id == 1) return "|U_{#mu4}|^{2}";
    if (id == 2) return "|U_{#tau4}|^{2}";
    return "#Deltam^{2} (eV^{2})";
  };

  // --- Auto-detect true dm2 from the ue42_um42 plane file ---
  double dm2_true = std::nan("");
  {
    std::ifstream f((prefix + "ue42_um42_nucut3" + suffix + ".txt").c_str());
    if (!f.is_open()) {
      std::cerr << "ERROR: cannot open " << (prefix + "ue42_um42_nucut3" + suffix + ".txt") << "\n";
      return;
    }
    // read the first non-empty, non-comment line
    std::string line;
    while (std::getline(f, line)) {
      if (line.empty() || line[0]=='#') continue;
      std::istringstream iss(line);
      double pars0, pars1, pars2, pars3, chi2, tgtchi2;
      if (iss >> pars0 >> pars1 >> pars2 >> pars3 >> chi2 >> tgtchi2) {
        dm2_true = pars3; // 4th number is dm2
        break;
      }
    }
    f.close();
  }
  if (!std::isfinite(dm2_true)) {
    std::cerr << "ERROR: failed to read a valid dm2 from ue42_um42 file\n";
    return;
  }

  // Map dm2 to label
  const double tol = 1e-6;
  std::string dm2_tag = "largeDm2"; // default
  if (std::fabs(dm2_true - 0.1) < tol)      dm2_tag = "smallDm2";
  else if (std::fabs(dm2_true - 6.0) < tol) dm2_tag = "medDm2";
  else if (std::fabs(dm2_true - 80.0) < tol)dm2_tag = "largeDm2";
  gSystem->mkdir(dm2_tag.c_str(), kTRUE);

  TFile* fout = TFile::Open(Form("%s/chi2_scans%s.root", dm2_tag.c_str(), suffix.c_str()), "RECREATE");
  if (!fout || fout->IsZombie()) {
  std::cerr << "ERROR: cannot create output ROOT file\n";
  return;
  }


  // True values (with auto-detected dm2)
  double fixedVals[4] = {UE42, UM42, UT42, dm2_true};

  // === Canvas with title band + 2×2 grid ===
  TCanvas* c_all = new TCanvas(Form("c_all_%s", dm2_tag.c_str()), "Chi2 scans", 1200, 1000);

  // Title pad (top)
  TPad* pTitle = new TPad("pTitle", "title", 0.00, 0.92, 1.00, 1.00);
  pTitle->SetFillStyle(0);
  pTitle->SetBorderSize(0);
  pTitle->Draw();

  // Grid pad (bottom)
  TPad* pGrid = new TPad("pGrid", "grid", 0.00, 0.00, 1.00, 0.92);
  pGrid->SetFillStyle(0);
  pGrid->SetBorderSize(0);
  pGrid->Draw();
  pGrid->cd();
  pGrid->Divide(2, 2, 0.003, 0.003);

  // Global title
  std::string globalTitle = Form(
    "True values: |U_{e4}|^{2} = %.2f, |U_{#mu4}|^{2} = %.2f, "
    "|U_{#tau4}|^{2} = %.1f, #Deltam^{2} = %.1f eV^{2}",
    fixedVals[0], fixedVals[1], fixedVals[2], fixedVals[3]
  );

  // If we’re using a suffix (e.g. "_noNuE"), annotate the title
  if (!suffix.empty()) {
    globalTitle += " (no #nu-e elastic)";
  }

  pTitle->cd();
  {
    TLatex t; t.SetNDC(); t.SetTextAlign(22); t.SetTextSize(0.4);
    t.DrawLatex(0.5, 0.50, globalTitle.c_str());
  }

  // Draw 4 subplots
  for (int i = 0; i < 4; ++i) {
    const auto& mode = modes[i];
    pGrid->cd(i + 1);

    // margins per subpad
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.15);

    std::ifstream file((prefix + mode.name + "_nucut3" + suffix + ".txt").c_str());
    if (!file.is_open()) {
      std::cerr << "Cannot open file: " << prefix + mode.name + "_nucut3" + suffix + ".txt" << std::endl;
      continue;
    }

    double xEdges[nBins+1], yEdges[nBins+1];
    logspace(minVals[mode.xidx], maxVals[mode.xidx], nBins, xEdges);
    logspace(minVals[mode.yidx], maxVals[mode.yidx], nBins, yEdges);

    TH2D* h2 = new TH2D(Form("h2_%s_%s", mode.name.c_str(), dm2_tag.c_str()), "",
                        nBins, xEdges, nBins, yEdges);
    h2->GetXaxis()->SetTitle(labelOf(mode.xidx).c_str());
    h2->GetYaxis()->SetTitle(labelOf(mode.yidx).c_str());
    h2->GetXaxis()->SetTitleOffset(1.10);
    h2->GetYaxis()->SetTitleOffset(1.10);

    double pars[4], chi2, tgtchi2;
    while (file >> pars[0] >> pars[1] >> pars[2] >> pars[3] >> chi2 >> tgtchi2) {
      const double dchi2 = chi2 - tgtchi2;
      h2->Fill(pars[mode.xidx], pars[mode.yidx], dchi2);
    }
    file.close();

    gPad->SetLogx();
    gPad->SetLogy();
    h2->SetMinimum(0.0);
    h2->SetMaximum(200.0);
    dunestyle::CenterTitles(h2);
    h2->Draw("COLZ");

    fout->cd();
    h2->Write();

    // DUNE tag only on top-left
    if (i == 0) {
      TLatex* sim = dunestyle::Simulation();
      if (sim) sim->SetTextSize(0.1);
    }

    // Red star at the true values in this plane
    TMarker* m = new TMarker(fixedVals[mode.xidx], fixedVals[mode.yidx], 29);
    m->SetMarkerColor(kRed);
    m->SetMarkerSize(2.0);
    m->Draw("same");
  }

  c_all->SaveAs(Form("%s/dchi2_scans_grid%s.png", dm2_tag.c_str(), suffix.c_str()));

  c_all->Close();
  fout->Write();
  fout->Close();
}
