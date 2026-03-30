#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>
#include <string>
#include <array>
#include <limits>
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
  const std::string suffix = "";  // leave "" if none
  // const std::string suffix = "_noNuE";  // leave "" if none

  struct Mode { std::string name; int xidx, yidx; };
  const std::vector<Mode> modes = {
    {"ue42_dm2", 0, 3}, // top-left:  y=dm2  vs x=|Ue4|^2
    {"um42_dm2", 1, 3}, // top-right: y=dm2  vs x=|Umu4|^2
    {"ut42_dm2", 2, 3}, // bot-left:  y=dm2  vs x=|Utau4|^2
    {"ue42_um42", 0, 1} // bot-right: y=|Umu4|^2 vs x=|Ue4|^2
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
  auto file_exists = [](const std::string& p) {
    std::ifstream f(p.c_str());
    return f.good();
  };

  // --- Determine dm2_tag from filenames like ue42_um42_medDm2[...].txt ---
  const std::string mode_ref = "ue42_um42";
  // const std::array<std::string,3> tags = {"smallDm2","medDm2","largeDm2"};
  // const std::array<std::string,1> tags{{"smallDm2"}};
  const std::array<std::string,1> tags{{"medDm2"}};
  // const std::array<std::string,1> tags{{"largeDm2"}};

  std::string dm2_tag, chosen_path;
  for (const auto& t : tags) {
    const std::string p1 = prefix + mode_ref + "_" + t + suffix + ".txt";
    const std::string p2 = prefix + mode_ref + "_" + t + ".txt";
    if (file_exists(p1)) { dm2_tag = t; chosen_path = p1; break; }
    if (file_exists(p2)) { dm2_tag = t; chosen_path = p2; break; }
  }
  if (dm2_tag.empty()) {
    std::cerr << "ERROR: cannot infer dm2_tag; none of "
              << mode_ref << "_{small,med,large}Dm2[<suffix>].txt exist.\n";
    return;
  }

  // Optional: validate tag vs dm2 value in file
  double dm2_true = std::numeric_limits<double>::quiet_NaN();
  {
    std::ifstream f(chosen_path.c_str());
    std::string line;
    while (std::getline(f, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream iss(line);
      double pars0, pars1, pars2, pars3, chi2, tgtchi2;
      if (iss >> pars0 >> pars1 >> pars2 >> pars3 >> chi2 >> tgtchi2) {
        dm2_true = pars3; // 4th number is dm2
        break;
      }
    }
  }
  const double dm2_expected =
      (dm2_tag == "smallDm2") ? 0.1 :
      (dm2_tag == "medDm2")   ? 6.0 :
      (dm2_tag == "largeDm2") ? 80.0 : std::numeric_limits<double>::quiet_NaN();
  const double tol = 1e-6;
  if (!std::isfinite(dm2_true) || std::fabs(dm2_true - dm2_expected) > tol) {
    std::cerr << "WARNING: tag/file mismatch: tag=" << dm2_tag
              << " implies dm2=" << dm2_expected
              << " but file shows dm2=" << dm2_true << "\n";
  }

  gSystem->mkdir(dm2_tag.c_str(), kTRUE);
  TFile* fout = TFile::Open(Form("dchi2_scans_%s%s.root", dm2_tag.c_str(), suffix.c_str()), "RECREATE");
  if (!fout || fout->IsZombie()) {
    std::cerr << "ERROR: cannot create output ROOT file\n";
    return;
  }

  // True values (with validated/parsed dm2)
  double fixedVals[4] = {UE42, UM42, UT42, dm2_true};

  // === Canvas with title band + 2×2 grid ===
  TCanvas* c_all = new TCanvas(Form("c_all_%s", dm2_tag.c_str()), "Chi2 scans", 1200, 1000);

  // Title pad (top)
  TPad* pTitle = new TPad("pTitle", "title", 0.00, 0.92, 1.00, 1.00);
  pTitle->SetFillStyle(0); pTitle->SetBorderSize(0); pTitle->Draw();

  // Grid pad (bottom)
  TPad* pGrid = new TPad("pGrid", "grid", 0.00, 0.00, 1.00, 0.92);
  pGrid->SetFillStyle(0); pGrid->SetBorderSize(0); pGrid->Draw();
  pGrid->cd(); pGrid->Divide(2, 2, 0.003, 0.003);

  // Global title
  std::string globalTitle = Form(
    "True values: |U_{e4}|^{2} = %.2f, |U_{#mu4}|^{2} = %.2f, "
    "|U_{#tau4}|^{2} = %.1f, #Deltam^{2} = %.1f eV^{2}",
    fixedVals[0], fixedVals[1], fixedVals[2], fixedVals[3]
  );
  if (!suffix.empty()) globalTitle += " (no #nu-e elastic)";

  pTitle->cd();
  { TLatex t; t.SetNDC(); t.SetTextAlign(22); t.SetTextSize(0.5); t.DrawLatex(0.5, 0.50, globalTitle.c_str()); }

  // Helper: find the best-existing path for a given mode
  auto mode_path = [&](const std::string& mode) {
    const std::string a = prefix + mode + "_" + dm2_tag + suffix + ".txt"; if (file_exists(a)) return a;
    const std::string b = prefix + mode + "_" + dm2_tag + ".txt";           if (file_exists(b)) return b;
    const std::string c = prefix + mode + "_nucut3" + suffix + ".txt";      if (file_exists(c)) return c; // legacy
    const std::string d = prefix + mode + "_nucut3.txt";                     return d;                      // last try
  };

  // Draw 4 subplots
  for (int i = 0; i < 4; ++i) {
    const auto& mode = modes[i];
    pGrid->cd(i + 1);

    // margins per subpad
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.14);
    gPad->SetLeftMargin(0.16);
    gPad->SetRightMargin(0.15);

    const std::string path = mode_path(mode.name);
    std::ifstream file(path.c_str());
    if (!file.is_open()) {
      std::cerr << "Cannot open file: " << path << std::endl;
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
    // gPad->SetLogz();
    h2->SetMinimum(0);
    h2->SetMaximum(180);
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

  c_all->SaveAs(Form("dchi2_scans_grid_%s%s.png", dm2_tag.c_str(), suffix.c_str()));

  c_all->Close();
  fout->Write();
  fout->Close();
}
