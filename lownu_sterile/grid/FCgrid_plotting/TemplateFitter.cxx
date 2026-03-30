#include "TemplateFitter.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TFile.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TDecompSVD.h"
#include "TDecompChol.h"
#include "TGraph.h"
#include "TLegend.h"
#include <TRandom3.h>
#include <TMatrixDEigen.h>
#include <TLine.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <DUNEStyle.h>
// using namespace std;
#include "TGaxis.h"
#include <vector>
#include <cassert>
#include <cmath>
#include <algorithm>
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include <iostream>
#include <iomanip>
#include <TString.h>
#include <TAxis.h>

namespace {  // file‑local helpers only (no header change)

  TH1D* MakeEmptyHist(const char* name, const std::vector<double>& edges) {
    auto* h = new TH1D(name, "", static_cast<int>(edges.size()) - 1, edges.data());
    h->SetDirectory(nullptr);
    return h;
  }

  static std::vector<double>
  BuildCompositeEdgesAuto(const std::vector<double>& EA,
                          const std::vector<double>& EB,
                          const std::vector<double>& EC,
                          double gap = 0.0)
  {
    auto span = [](const std::vector<double>& E){ return E.back() - E.front(); };
    const double offA = 0.0;
    const double offB = offA + span(EA) + gap;
    const double offC = offB + span(EB) + gap;

    std::vector<double> edges;
    edges.reserve( (EA.size()-1) + (EB.size()-1) + (EC.size()-1) + 1 );

    edges.insert(edges.end(), EA.begin(), EA.end());
    for (size_t i = 1; i < EB.size(); ++i) edges.push_back(EB[i] + offB);
    for (size_t i = 1; i < EC.size(); ++i) edges.push_back(EC[i] + offC);
    return edges;
  }

  // NEW: generic emptiness check (TH1/TH2/TH3). Treats nullptr as empty.
  inline bool IsEmpty(const TH1* h, double eps = 0.0, bool include_under_over = true) {
    if (!h) return true;

    const int nx = h->GetNbinsX();
    const int ny = h->GetNbinsY();
    const int nz = h->GetNbinsZ();

    const int x0 = include_under_over ? 0 : 1;
    const int x1 = include_under_over ? nx + 1 : nx;
    const int y0 = (ny > 0 && include_under_over) ? 0 : 1;
    const int y1 = (ny > 0 && include_under_over) ? ny + 1 : (ny > 0 ? ny : 1);
    const int z0 = (nz > 0 && include_under_over) ? 0 : 1;
    const int z1 = (nz > 0 && include_under_over) ? nz + 1 : (nz > 0 ? nz : 1);

    for (int ix = x0; ix <= x1; ++ix)
      for (int iy = y0; iy <= y1; ++iy)
        for (int iz = z0; iz <= z1; ++iz)
          if (std::abs(h->GetBinContent(ix, iy, iz)) > eps)
            return false;

    return true;
  }

} // namespace


// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}

TemplateFitter::TemplateFitter(TemplateGroup_CC (&ccTemplates)[nbins_Ev], TemplateGroup_nue (&nueTemplates)[nbins_Ev], const std::vector<int>& dropBins) 
: dropBins_(dropBins)
{
  for (int i = 0; i < nbins_Ev; ++i) {
    CC[i] = ccTemplates[i];
    nue[i] = nueTemplates[i];
  }
}


void TemplateFitter::setEnergyBins(double bins[nbins_Ev + 1]) {
  for (int i = 0; i <= nbins_Ev; ++i)
    binEdges[i] = bins[i];
}

void TemplateFitter::setCovmtr(const CovarianceMatrix& covInput) {
  cov = covInput;
  for (int i = 0; i < nbins; ++i) {
    for (int j = 0; j < nbins; ++j) {
      flmx(i, j) = cov.FluxMx[i][j];
      sigmx(i, j) = cov.CrossSectionMx[i][j];
      detmx(i, j) = cov.DetectorMx[i][j];
    }
  }
  sysmx = flmx + detmx + sigmx;
}

void TemplateFitter::setLEParameters(LE_FitParameters LEPars) {
  LEfit = LEPars;
}


double TemplateFitter::getAvgPme( double energy, double ue42, double um42, double ut42, double dm2, double ft[7] ) const 
{
  if(dm2<par_min[3]) return 0.;
  else{
  double k = 1.27*dm2/energy;
  double L0=0.34, L1=0.35, L2=0.55, L3=0.6, L;
  double a=ft[0], b=ft[1], c=ft[2], d=ft[3], avg1=ft[4], avg2=ft[5], norm=ft[6];
  double A = 4 * ue42 * um42;

  L=L1;
  double prob_u = - avg1 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);
  L=L0;
  double prob_l = - avg1 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);

  double prob1 = prob_u - prob_l;


  L=L2;
  double term1_u = 2*pow(k,4) * L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_u = 6*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_u = 3*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );
  L=L1;
  double term1_l = 2*pow(k,4) * L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_l = 6*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_l = 3*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );

  prob_u = A/(48*pow(k,4)) * (term1_u - term2_u - term3_u);
  prob_l = A/(48*pow(k,4)) * (term1_l - term2_l - term3_l);

  double prob2 = prob_u - prob_l;

  
  L=L3;
  prob_u = - avg2 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);
  L=L2;
  prob_l = - avg2 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);

  double prob3 = prob_u - prob_l;

  double prob = norm*(prob1 + prob2 + prob3);

  return prob;}
}


double TemplateFitter::getAvgPmt( double energy, double ue42, double um42, double ut42, double dm2, double ft[7] ) const
{
  if(dm2<par_min[3]) return 0.;
  else{
  double k = 1.27*dm2/energy;
  double L0=0.34, L1=0.35, L2=0.55, L3=0.6, L;
  double a=ft[0], b=ft[1], c=ft[2], d=ft[3], avg1=ft[4], avg2=ft[5], norm=ft[6];
  double A = 4 * ut42 * um42;

  L=L1;
  double prob_u = - avg1 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);
  L=L0;
  double prob_l = - avg1 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);

  double prob1 = prob_u - prob_l;


  L=L2;
  double term1_u = 2*pow(k,4) * L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_u = 6*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_u = 3*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );
  L=L1;
  double term1_l = 2*pow(k,4) * L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_l = 6*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_l = 3*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );

  prob_u = A/(48*pow(k,4)) * (term1_u - term2_u - term3_u);
  prob_l = A/(48*pow(k,4)) * (term1_l - term2_l - term3_l);

  double prob2 = prob_u - prob_l;

  
  L=L3;
  prob_u = - avg2 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);
  L=L2;
  prob_l = - avg2 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);

  double prob3 = prob_u - prob_l;

  double prob = norm*(prob1 + prob2 + prob3);

  return prob;}
}

double TemplateFitter::getAvgPet( double energy, double ue42, double um42, double ut42, double dm2, double ft[7] ) const
{
  if(dm2<=par_min[3]) return 0.;
  else{
  double k = 1.27*dm2/energy;
  double L0=0.34, L1=0.35, L2=0.55, L3=0.6, L;
  double a=ft[0], b=ft[1], c=ft[2], d=ft[3], avg1=ft[4], avg2=ft[5], norm=ft[6];
  double A = 4 * ut42 * ue42;

  L=L1;
  double prob_u = - avg1 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);
  L=L0;
  double prob_l = - avg1 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);

  double prob1 = prob_u - prob_l;


  L=L2;
  double term1_u = 2*pow(k,4) * L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_u = 6*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_u = 3*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );
  L=L1;
  double term1_l = 2*pow(k,4) * L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_l = 6*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_l = 3*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );

  prob_u = A/(48*pow(k,4)) * (term1_u - term2_u - term3_u);
  prob_l = A/(48*pow(k,4)) * (term1_l - term2_l - term3_l);

  double prob2 = prob_u - prob_l;

  
  L=L3;
  prob_u = - avg2 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);
  L=L2;
  prob_l = - avg2 * A * ( sin(2*k*L) - 2*k*L ) / (4*k);

  double prob3 = prob_u - prob_l;

  double prob = norm*(prob1 + prob2 + prob3);

  return prob;}
}

double TemplateFitter::getAvgPee( double energy, double ue42, double um42, double ut42, double dm2, double ft[7] ) const
{
  if(dm2<par_min[3]) return 1.;
  else{
  double k = 1.27*dm2/energy;
  double L0=0.34, L1=0.35, L2=0.55, L3=0.6, L;
  double a=ft[0], b=ft[1], c=ft[2], d=ft[3], avg1=ft[4], avg2=ft[5], norm=ft[6];
  double A = 4 * ue42 * (1.-ue42);

  L=L1;
  double prob_u = avg1 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L0;
  double prob_l = avg1 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);

  double prob1 = prob_u - prob_l;


  L=L2;
  double term1_u = -2*(A-2)*pow(k,4)*L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_u = 6*A*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_u = 3*A*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );
  L=L1;
  double term1_l = -2*(A-2)*pow(k,4)*L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_l = 6*A*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_l = 3*A*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );

  prob_u = 1/(48*pow(k,4)) * (term1_u + term2_u + term3_u);
  prob_l = 1/(48*pow(k,4)) * (term1_l + term2_l + term3_l);

  double prob2 = prob_u - prob_l;


  L=L3;
  prob_u = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L2;
  prob_l = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);

  double prob3 = prob_u - prob_l;

  double prob = norm*(prob1 + prob2 + prob3);

  return prob;}
}


double TemplateFitter::getAvgPmm( double energy, double ue42, double um42, double ut42, double dm2, double ft[7] ) const
{
  if(dm2<par_min[3]) return 1.;
  else{
  double k = 1.27*dm2/energy;
  double L0=0.34, L1=0.35, L2=0.55, L3=0.6, L;
  double a=ft[0], b=ft[1], c=ft[2], d=ft[3], avg1=ft[4], avg2=ft[5], norm=ft[6];
  double A = 4 * um42 * (1.-um42);

  L=L1;
  double prob_u = avg1 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L0;
  double prob_l = avg1 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);

  double prob1 = prob_u - prob_l;


  L=L2;
  double term1_u = -2*(A-2)*pow(k,4)*L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_u = 6*A*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_u = 3*A*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );
  L=L1;
  double term1_l = -2*(A-2)*pow(k,4)*L * ( 12*a + 6*b*L + 4*c*L*L + 3*d*L*L*L );
  double term2_l = 6*A*k*sin(2*k*L) * ( 2*a*k*k + 2*b*k*k*L + 2*c*k*k*L*L - c + 2*d*k*k*L*L*L - 3*d*L );
  double term3_l = 3*A*cos(2*k*L) * ( 2*k*k*(b+2*c*L) + d*(6*k*k*L*L-3) );

  prob_u = 1/(48*pow(k,4)) * (term1_u + term2_u + term3_u);
  prob_l = 1/(48*pow(k,4)) * (term1_l + term2_l + term3_l);

  double prob2 = prob_u - prob_l;


  L=L3;
  prob_u = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);
  L=L2;
  prob_l = avg2 * ( A*sin(2*k*L) - 2*(A-2)*k*L ) / (4*k);

  double prob3 = prob_u - prob_l;

  double prob = norm*(prob1 + prob2 + prob3);

  return prob;}
}



// void TemplateFitter::plotOscillationComponents(TCanvas* c, const double par[4],
//                                                TH1D* CC_tp_mm, TH1D* CC_tp_me, TH1D* CC_tp_ee, TH1D* CC_tp_em,
//                                                TH1D* nue_tp_mm, TH1D* nue_tp_me, TH1D* nue_tp_ee, TH1D* nue_tp_em,
//                                                bool plotBestfit) const
void TemplateFitter::plotOscillationComponents(
    TCanvas* c, const double par[4],
    // CCνμ (μ-like sample)
    TH1D* CCm_sig_em, TH1D* CCm_sig_mm,
    // CCνμ background components (from NC)
    TH1D* CCm_bkg_nc,
    // CCνe signal (e-like sample)
    TH1D* CCe_sig_me, TH1D* CCe_sig_ee,
    // CCνe background components (from NC and from ν–e leakage)
    TH1D *CCe_bkg_nc, TH1D* CCe_bkg_me, TH1D* CCe_bkg_ee, TH1D* CCe_bkg_em, TH1D* CCe_bkg_mm, 
    // ν–e signal (inverse Eθ² selection)
    TH1D* nue_sig_mm, TH1D* nue_sig_me, TH1D* nue_sig_ee, TH1D* nue_sig_em,
    // ν–e background (from CCνe leakage into ν–e)
    TH1D* nue_bkg_ee, TH1D* nue_bkg_em,
    bool plotBestfit) const
{
  // Build CC lepton-energy edges (54 bins → 55 edges)
  assert(nbinsCC == 54 && "This binning logic assumes 54 CC bins");
  assert(nbins_nue == 7  && "Expecting 7 bins for nu+e segment");

  std::vector<double> CCEdges(nbinsCC + 1);
  CCEdges[0] = 0.0;
  for (int i = 0; i < nbinsCC; ++i) {
    double step = 0.0;
    if      (i < 1)       step = 0.45;
    else if (i < 35)      step = 0.10;
    else if (i < 36)      step = 0.15;
    else if (i < 41)      step = 0.20;
    else if (i < 46)      step = 0.40;
    else if (i < 51)      step = 0.80;
    else if (i < 53)      step = 1.50;
    else                  step = 2.00;
    CCEdges[i + 1] = CCEdges[i] + step;
    // std::cout << i << "\t" << CCEdges[i] << "\t" << CCEdges[i + 1] << "\n";
  }
  // sanity
  assert(std::abs(CCEdges.front()) < 1e-12);
  assert(std::abs(CCEdges.back() - 16.0) < 1e-9);

  // EA and EB: identical CC edges
  // std::vector<double> EA = CCEdges;
  const double cutE = 5.0;
  int i0 = int(std::lower_bound(CCEdges.begin(), CCEdges.end(), cutE) - CCEdges.begin()); // -> 41
  // A segment edges: 0..5.0 (inclusive), so take edges[0..i0]
  std::vector<double> EA_trunc(CCEdges.begin(), CCEdges.begin() + i0 + 1);

  std::vector<double> EB = CCEdges;

  // EC (ν+e) edges — use your actual list here
  const double edgesC_arr[8] = { 0.3, 0.6, 0.92, 1.3, 1.75, 2.45, 3.9, 16.0 };
  std::vector<double> EC(edgesC_arr, edgesC_arr + 8);

  // more sanity
  // assert(static_cast<int>(EA.size()) == nbinsCC + 1);
  assert(int(EA_trunc.size()) == i0 + 1);
  assert(std::abs(EA_trunc.front() - 0.0) < 1e-12 && std::abs(EA_trunc.back() - 5.0) < 1e-9);
  assert(int(EB.size()) == nbinsCC + 1);
  assert(int(EC.size()) == nbins_nue + 1);

  // Final stitched axis (A |gap| B |gap| C)
  const double Emax = 16.0;
  const double gap  = 0.0;
  // std::vector<double> edges = BuildCompositeEdges(EA, EB, EC, Emax, gap);
  std::vector<double> edges = BuildCompositeEdgesAuto(EA_trunc, EB, EC, gap);
  // assert(static_cast<int>(edges.size()) == nbins + 1);


  // // ---- 2) Targets on energy axis ----
  // TH1D* tgtE   = MakeEmptyHist("tgtE",   edges);
  // TH1D* FCtgtE = MakeEmptyHist("FCtgtE", edges);
  // for (int i = 0; i < nbins; ++i) {
  //   tgtE  ->SetBinContent(i + 1, this->target(i, 0));
  //   FCtgtE->SetBinContent(i + 1, this->FCtarget(i, 0));
  // }

  // ---- 3) Components on energy axis (copy contents by index) ----
  // TH1D* h_mmE = MakeEmptyHist("h_mmE", edges);
  // TH1D* h_meE = MakeEmptyHist("h_meE", edges);
  // TH1D* h_eeE = MakeEmptyHist("h_eeE", edges);
  // TH1D* h_emE = MakeEmptyHist("h_emE", edges);
  // TH1D* h_ncE = MakeEmptyHist("h_ncE", edges);


  TH1D* hCCm_sig    = MakeEmptyHist("hCCm_sig", edges);
  TH1D* hCCm_bkg_NC = MakeEmptyHist("hCCm_bkg", edges);
  TH1D* hCCe_sig    = MakeEmptyHist("hCCe_sig", edges);
  TH1D* hCCe_bkg_NC = MakeEmptyHist("hCCe_bkg", edges); 
  TH1D* hCCe_bkg    = MakeEmptyHist("hCCe_bkg_nue", edges); 
  TH1D* hnue_sig    = MakeEmptyHist("hnue_sig", edges);
  TH1D* hnue_bkg    = MakeEmptyHist("hnue_bkg", edges); 

  TH1D* hOsc_me     = MakeEmptyHist("hOsc_me", edges);
  TH1D* hOsc_em     = MakeEmptyHist("hOsc_em", edges);

  TH1D* hRatio      = MakeEmptyHist("hRatio", edges);

  // --- compute sizes/offsets after truncation ---
  const int nA = static_cast<int>(EA_trunc.size()) - 1;  // bins in first segment (e.g. 41)
  const int nB = nbinsCC;                                // 54
  const int nC = nbins_nue;                              // 8
  const int nbins_total = nA + nB + nC;
  auto map_to_original = [&](int i_new)->int {
    // original stitching was [A:54][B:54][C:nbins_nue]
    // we removed (54 - nA) bins from A, so bins after A shift forward by that amount
    return (i_new < nA) ? i_new : i_new + (nbinsCC - nA);
  };

  // sanity with the histogram you created from `edges`
  assert(hCCm_sig->GetNbinsX() == nbins_total);

  // --- fill components by segment with correct indices ---
  for (int i = 0; i < nbins_total; ++i) {
    // fill ratio with correct original bin index
    const int j = map_to_original(i); // 0-based index in original stitched hist
    hRatio->SetBinContent(i + 1, hRatioToNull->GetBinContent(j + 1));
    // if (i == 0 || i == nA || i == nA+nB) {
    //   std::cout << "i_new=" << i << " -> j_orig=" << j << "\n";
    // }

    if (i < nA) {
      // segment A: 0 .. nA-1  (use CC bins 1..nA)
      const int idx = i;                  // 0-based into source
      double sig_mm = CCm_sig_mm->GetBinContent(idx + 1);
      double sig_em = CCm_sig_em->GetBinContent(idx + 1);
      hCCm_sig->SetBinContent(i + 1, sig_mm);
      hOsc_em->SetBinContent(i + 1, sig_em);
      
      hCCm_bkg_NC->SetBinContent(i + 1, CCm_bkg_nc->GetBinContent(idx+1));

      // hRatio->SetBinContent(i + 1, hRatioToNull->GetBinContent(idx+1));
    }
    else if (i < nA + nB) {
      // segment B: nA .. nA+nB-1  (use full CC B bins 1..54)
      const int idx = i - nA;             // 0..nB-1
      // signal in CC νe
      double sig_ee = CCe_sig_ee->GetBinContent(idx+1); // νe→νe
      double sig_me = CCe_sig_me->GetBinContent(idx+1); // νμ→νe

      // background leaking into CC νe (from ν–e selection)
      double bkg_mm = CCe_bkg_mm->GetBinContent(idx+1); // νμ-comp
      double bkg_em = CCe_bkg_em->GetBinContent(idx+1); // νe→νμ
      double bkg_ee = CCe_bkg_ee->GetBinContent(idx+1); // νe→νe
      double bkg_me = CCe_bkg_me->GetBinContent(idx+1); // νμ→νe

      // background from NC interaction
      double bkg_nc = CCe_bkg_nc->GetBinContent(idx+1); // νμ→νe

      // fill by oscillation channel
      hCCe_sig->SetBinContent(i+1, sig_ee);
      hOsc_me->SetBinContent(i+1, sig_me);
      hCCe_bkg->SetBinContent(i+1, bkg_mm + bkg_me + bkg_em + bkg_ee);
      hCCe_bkg_NC->SetBinContent(i+1, bkg_nc); // background contributes to em channel in CC νe region

      // hRatio->SetBinContent(i + 1, hRatioToNull->GetBinContent(idx+1));

    }
    else {
      // segment C: nA+nB .. end  (use nu+e bins 1..8)
      const int idx = i - nA - nB;        // 0..nC-1
      // signal in ν–e
      double s_mm = nue_sig_mm->GetBinContent(idx+1);
      double s_me = nue_sig_me->GetBinContent(idx+1);
      double s_ee = nue_sig_ee->GetBinContent(idx+1);
      double s_em = nue_sig_em->GetBinContent(idx+1);

      // background into ν–e (from CC νe failing inverse Eθ² cut)
      double b_ee = nue_bkg_ee->GetBinContent(idx+1);
      double b_em = nue_bkg_em->GetBinContent(idx+1);

      // fill by oscillation channel
      hnue_sig->SetBinContent(i+1, s_ee + s_mm);
      hOsc_me->SetBinContent(i+1, s_me);
      hOsc_em->SetBinContent(i+1, s_em);
      hnue_bkg->SetBinContent(i+1, b_em + b_ee);

      // hRatio->SetBinContent(i + 1, hRatioToNull->GetBinContent(idx+1));

    }
  }

  // std::cout << "Original ratio bin 55 = "
  //           << hRatioToNull->GetBinContent(55) << "\n";
  // std::cout << "Original ratio bin 56 = "
  //           << hRatioToNull->GetBinContent(56) << "\n";

  // ---- 2) Targets on energy axis ----
  TH1D* tgtE   = MakeEmptyHist("tgtE",   edges);
  TH1D* FCtgtE = MakeEmptyHist("FCtgtE", edges);
  // --- OPTIONAL: remap targets from original 116-bin stitching ---
  // original order was: [A:54][B:54][C:8]; we removed (54 - nA) bins from A
  // so new bin i corresponds to original index j = i (in A) or i + (54 - nA) (in B or C)
  for (int i = 0; i < nbins_total; ++i) {
    int j = (i < nA) ? i : i + (nbinsCC - nA);  // skip removed A bins
    tgtE  ->SetBinContent(i + 1, this->target(j, 0));
    FCtgtE->SetBinContent(i + 1, this->FCtarget(j, 0));
  }

  auto is_empty = [](const TH1* h) {
    if (!h) return true;
    const int n = h->GetNbinsX();
    for (int b = 0; b <= n + 1; ++b) {          // include under/overflow
      if (h->GetBinContent(b) != 0.0) return false;
    }
    return true;
  };

  // ---- 4) Styling ----
  auto color_CCm = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 1);
  auto color_CCe = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 2);
  auto color_nue = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 3);
  auto color_nc  = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 4);
  auto color_me  = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 5);
  auto color_em  = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 6);

  auto SetColor = [](TH1* h, Color_t color, double alpha = 1.0, double lw = 0.5) {
    h->Scale(1.0, "width");
    h->SetFillColorAlpha(color, alpha);
    h->SetLineColor(color);
    h->SetLineWidth(lw);
  };

  SetColor(hCCm_sig,      color_CCm);
  SetColor(hCCm_bkg_NC,   color_nc);
  SetColor(hCCe_sig,      color_CCe);
  SetColor(hCCe_bkg_NC,   color_nc);
  SetColor(hCCe_bkg,      color_nue);
  SetColor(hnue_sig,      color_nue);
  SetColor(hnue_bkg,      color_CCe);
  SetColor(hOsc_me,       color_me);
  SetColor(hOsc_em,       color_em);

  // ---- 5) Stack & draw ----
  THStack* h_stack = new THStack(
      "h_stack",
      Form("Oscillation Hypothesis: U_{e4}^{2} = %.3f, U_{#mu4}^{2} = %.3f, U_{#tau4}^{2} = %.1f, #Deltam^{2} = %.1f eV^{2} (#nu < %.1f GeV)",
          par[0], par[1], par[2], par[3], kNuCuts[kNuCut]));

  // gStyle->SetOptTitle(1);            // make sure title is on
  // gStyle->SetTitleFontSize(0.1);   // adjust size globally for titles

  // This order provides plots with different dm2 etc
  h_stack->Add(hOsc_em);
  h_stack->Add(hOsc_me);
  h_stack->Add(hCCe_bkg);
  h_stack->Add(hCCm_bkg_NC);
  h_stack->Add(hCCe_bkg_NC);
  h_stack->Add(hnue_bkg);
  h_stack->Add(hCCm_sig);
  h_stack->Add(hCCe_sig);
  h_stack->Add(hnue_sig);
  
  // // This order provides plots with small contributions at the bottom
  // h_stack->Add(hCCe_bkg);
  // h_stack->Add(hOsc_em);
  // h_stack->Add(hCCm_bkg_NC);
  // h_stack->Add(hCCe_bkg_NC);
  // h_stack->Add(hOsc_me);
  // h_stack->Add(hCCm_sig);
  // h_stack->Add(hCCe_sig);
  // h_stack->Add(hnue_bkg);
  // h_stack->Add(hnue_sig);


  // NC components
  double nc_int =
      hCCm_bkg_NC->Integral(0, hCCm_bkg_NC->GetNbinsX()+1) +
      hCCe_bkg_NC->Integral(0, hCCe_bkg_NC->GetNbinsX()+1);

  // Total = sum of all components in the stack
  double tot_int =
      hCCe_bkg->Integral(0, hCCe_bkg->GetNbinsX()+1) +
      hOsc_em->Integral(0, hOsc_em->GetNbinsX()+1) +
      hCCm_bkg_NC->Integral(0, hCCm_bkg_NC->GetNbinsX()+1) +
      hCCe_bkg_NC->Integral(0, hCCe_bkg_NC->GetNbinsX()+1) +
      hOsc_me->Integral(0, hOsc_me->GetNbinsX()+1) +
      hCCm_sig->Integral(0, hCCm_sig->GetNbinsX()+1) +
      hCCe_sig->Integral(0, hCCe_sig->GetNbinsX()+1) +
      hnue_bkg->Integral(0, hnue_bkg->GetNbinsX()+1) +
      hnue_sig->Integral(0, hnue_sig->GetNbinsX()+1);

  double ratio = (tot_int > 0.0) ? nc_int / tot_int : 0.0;

  std::cout << "NC / Total = " << ratio << std::endl;

  // Nice log-range (adjust as you like)
  h_stack->SetMinimum(20.);
  h_stack->SetMaximum(h_stack->GetMaximum() * 8.);


  // ---------------------- Canvas & pads (contiguous) ----------------------
  c->cd();
  c->Clear();
  c->SetMargin(0,0,0,0);
  c->SetBorderMode(0);

  TPad* pTop = new TPad("pTop","", 0.00, 0.26, 1.00, 1.00);  // top 70%
  TPad* pBot = new TPad("pBot","", 0.00, 0.00, 1.00, 0.26);  // bottom 30%

  for (auto* p : {pTop, pBot}) {
    p->SetBorderMode(0);
    p->SetFrameBorderMode(0);
    p->SetFillStyle(4000);   // transparent, avoids seam
    p->SetLeftMargin(0.12);
    p->SetRightMargin(0.05);
  }
  pTop->SetTopMargin(0.08);
  pTop->SetBottomMargin(0.00);  // touch bottom pad
  pTop->SetLogy();

  pBot->SetTopMargin(0.00);     // touch top pad
  pBot->SetBottomMargin(0.36);  // room for labels + title
  pBot->SetGridy(true);

  pTop->Draw();
  pBot->Draw();

  // ---------------------- TOP pad: stack only (no x-title) ----------------------
  pTop->cd();

  // Reasonable log range
  h_stack->SetMinimum(20.);
  h_stack->SetMaximum(h_stack->GetMaximum() * 8.);
  h_stack->Draw("hist");

  // Remove native x axis text and title in the TOP pad
  if (TH1* frameTop = h_stack->GetHistogram()) {
    frameTop->GetXaxis()->SetLabelSize(0);
    frameTop->GetXaxis()->SetTitleSize(0);   // <-- hide x-axis title on top
    frameTop->GetXaxis()->SetTickLength(0);
    frameTop->GetYaxis()->SetTitle("entries / GeV / yr.POT");
  }

  // Title box polish (optional)
  gPad->Update();
  if (auto* pt = (TPaveText*)gPad->GetPrimitive("title")) {
    pt->SetTextSize(0.044);
    pt->SetBorderSize(0);
    pt->SetFillStyle(0);
    pt->SetX1NDC(0.12); pt->SetX2NDC(0.88);
    pt->SetY1NDC(0.93); pt->SetY2NDC(0.99);
    gPad->Modified(); gPad->Update();
  }

  // --- segment geometry and separators (A|B and B|C) ---
  auto span = [](const std::vector<double>& E){ return E.back() - E.front(); };
  // const double gap  = 0.0;      // your gap
  const double spanA = span(EA_trunc), spanB = span(EB), spanC = span(EC);
  const double offA  = 0.0;
  const double offB  = offA + spanA + gap;
  const double offC  = offB + spanB + gap;

  // vertical separators in TOP pad (from bottom to top of pad)
  const double yBottomTop = std::pow(10.0, gPad->GetUymin());
  const double yTopTop    = std::pow(10.0, gPad->GetUymax());

  auto* L1_top = new TLine(offA + spanA, yBottomTop, offA + spanA, yTopTop);
  auto* L2_top = new TLine(offB + spanB, yBottomTop, offB + spanB, yTopTop);
  for (TLine* L : {L1_top, L2_top}) { L->SetLineColor(kBlack); L->SetLineWidth(2); L->SetLineStyle(2); L->Draw(); }

  // Segment labels (optional)
  const double yMidTop = std::sqrt(yBottomTop * yTopTop);
  TLatex latex; latex.SetNDC(false); latex.SetTextAlign(22); latex.SetTextSize(0.045);
  latex.DrawLatex(offA + 0.5*spanA, yMidTop, "#nu_{#mu}-CC");
  latex.DrawLatex(offB + 0.5*spanB, yMidTop, "#nu_{e}-CC");
  latex.DrawLatex(offC + 0.5*spanC, yMidTop, "#nu+e");

  // Legend (unchanged)
  TLegend* leg = MakeLegend(0.60, 0.60, 0.85, 0.88);
  leg->AddEntry(hCCm_sig, "Inclusive CC #nu_{#mu}");
  leg->AddEntry(hCCe_sig, "Inclusive CC #nu_{e}");
  leg->AddEntry(hCCm_bkg_NC, "NC background");
  leg->AddEntry(hnue_sig, "#nu + e elastic");
  if (hOsc_em && !IsEmpty(hOsc_em)) leg->AddEntry(hOsc_em, "Oscillated #nu_{e} #rightarrow #nu_{#mu}");
  if (hOsc_me && !IsEmpty(hOsc_me)) leg->AddEntry(hOsc_me, "Oscillated #nu_{#mu} #rightarrow #nu_{e}");
  leg->Draw();

  dunestyle::CenterTitles(h_stack->GetHistogram());
  auto* sim = dunestyle::Simulation();
  sim->SetTextSize(0.08);
  gPad->RedrawAxis();

  // ---------------------- BOTTOM pad: ratio + the ONLY x-title ----------------------
  pBot->cd();

  // Make a drawable copy of stitched ratio hRatio
  TH1D* r = (TH1D*)hRatio->Clone("hRatio_draw");
  r->SetDirectory(nullptr);
  r->SetTitle("");
  r->SetLineWidth(2);
  r->SetMarkerStyle(20);
  r->SetMarkerSize(0.7);

  // Y around 1
  r->SetMinimum(0.8);
  r->SetMaximum(1.2);

  // This pad owns the ONLY x-axis title
  r->GetXaxis()->SetTitle("Reconstructed lepton energy [GeV]");
  r->GetXaxis()->SetTitleSize(0.15);
  r->GetXaxis()->SetTitleOffset(1.1);
  // We'll draw custom composite labels, so hide native labels/ticks
  r->GetXaxis()->SetLabelSize(0);
  r->GetXaxis()->SetTickLength(0);

  r->GetYaxis()->SetTitle("Ratio To Null");
  r->GetYaxis()->SetNdivisions(505);
  r->GetYaxis()->SetTitleSize(0.12);
  r->GetYaxis()->SetTitleOffset(0.4);
  r->GetYaxis()->SetLabelSize(0.10);

  dunestyle::CenterTitles(r);
  r->Draw("hist");

  // Unity line
  auto* unity = new TLine(edges.front(), 1.0, edges.back(), 1.0);
  unity->SetLineStyle(2);
  unity->SetLineWidth(2);
  unity->SetLineColor(kBlack);
  unity->Draw("same");

  // Composite x-axis labels in bottom pad only
  const double y0_bot  = r->GetMinimum();    // linear coords here
  const int    nDivA   = 203;
  const int    nDivBC  = 505;
  // const double Emax    = 16.0;
  const double eps0    = 1e-3;

  auto* axA2 = new TGaxis(offA, y0_bot, offA + spanA, y0_bot, 0.0, spanA, nDivA, "+");
  axA2->SetLabelSize(0.10); axA2->SetLabelOffset(0.02); axA2->SetTickSize(0.06); axA2->Draw();

  auto* axB2 = new TGaxis(offB, y0_bot, offB + spanB, y0_bot, 0.0, Emax, nDivBC, "+");
  axB2->SetLabelSize(0.10); axB2->SetLabelOffset(0.02); axB2->SetTickSize(0.06); axB2->Draw();

  auto* axC2 = new TGaxis(offC, y0_bot, offC + spanC, y0_bot, eps0, Emax, nDivBC, "+");
  axC2->SetLabelSize(0.10); axC2->SetLabelOffset(0.02); axC2->SetTickSize(0.06); axC2->Draw();

  // Draw vertical separators in BOTTOM pad too (so they span both pads)
  const double yBottomBot = r->GetMinimum();
  const double yTopBot    = r->GetMaximum();

  auto* L1_bot = new TLine(offA + spanA, yBottomBot, offA + spanA, yTopBot);
  auto* L2_bot = new TLine(offB + spanB, yBottomBot, offB + spanB, yTopBot);
  for (TLine* L : {L1_bot, L2_bot}) { L->SetLineColor(kBlack); L->SetLineWidth(2); L->SetLineStyle(2); L->Draw(); }

  pBot->Update();
  c->Update();

  
}



// Build keep indices from drop list
inline std::vector<int> BuildKeepIndices(int N, std::vector<int> drop) {
    std::sort(drop.begin(), drop.end());
    drop.erase(std::unique(drop.begin(), drop.end()), drop.end());
    for (int d : drop) {
        if (d < 0 || d >= N)
            throw std::out_of_range("drop index out of range");
    }
    std::vector<int> keep;
    keep.reserve(N - (int)drop.size());
    int j = 0;
    for (int i = 0; i < N; ++i) {
        if (j < (int)drop.size() && drop[j] == i) { ++j; continue; }
        keep.push_back(i);
    }
    return keep;
}

// Reduce N×N or N×1
inline TMatrixD ReduceMatrix(const TMatrixD& A,
                            const std::vector<int>& drop_idx)
{
    const int M = A.GetNrows();
    const int N = A.GetNcols();

    // Column vector case: N×1
    if (N == 1) {
        auto keep = BuildKeepIndices(M, drop_idx);
        TMatrixD R((int)keep.size(), 1);
        for (int i = 0; i < (int)keep.size(); ++i)
            R(i,0) = A(keep[i], 0);
        return R;
    }

    // Square matrix case: N×N
    if (M == N) {
        auto keep = BuildKeepIndices(N, drop_idx);
        TMatrixD R((int)keep.size(), (int)keep.size());
        for (int i = 0; i < (int)keep.size(); ++i)
            for (int j = 0; j < (int)keep.size(); ++j)
                R(i,j) = A(keep[i], keep[j]);
        return R;
    }

    throw std::invalid_argument("ReduceMatrix: input must be NxN or Nx1");
}



void DrawMatrix(const TMatrixD& M,
                const char* title = nullptr,  // 1D: "Vector;Index;Value"  |  2D: "Matrix;Column j;Row i"
                const char* out   = nullptr,  // e.g. "plot.png" (nullptr = don't save)
                bool autoLog      = true)
{
  const int nRows = M.GetNrows();
  const int nCols = M.GetNcols();

  static int uid = 0;
  const TString cname = Form("c_matrix_%d", uid);
  const TString hname = Form("h_matrix_%d", uid);
  ++uid;

  gStyle->SetOptStat(0);

  // --- 1D case: vector (Nx1 or 1xN) ---
  if (nCols == 1 || nRows == 1) {
    const int N = (nCols == 1) ? nRows : nCols;

    auto* c = new TCanvas(cname, cname, 900, 500);
    c->SetRightMargin(0.06);

    auto* h = new TH1D(hname,
                       title ? title : "Vector;Index;Value",
                       N, -0.5, N - 0.5);

    double minVal = 0.0, maxVal = 0.0;
    bool first = true;
    for (int i = 0; i < N; ++i) {
      const double v = (nCols == 1) ? M(i,0) : M(0,i);
      h->SetBinContent(i + 1, v);
      std::cout << i+1 << "\t" << v << "\n";
      if (first) { minVal = maxVal = v; first = false; }
      else { if (v < minVal) minVal = v; if (v > maxVal) maxVal = v; }
    }

    // Cosmetics
    h->GetXaxis()->CenterTitle(true);
    h->GetYaxis()->CenterTitle(true);
    h->SetLineWidth(2);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.7);

    if (autoLog && minVal > 0.0) c->SetLogy();
    h->Draw("HIST");  // or "E1" if you prefer points

    if (out && out[0]) c->SaveAs(out);
    return;
  }

  // --- 2D case: matrix (nRows x nCols) ---
  auto* c = new TCanvas(cname, cname, 950, 800);
  c->SetRightMargin(0.15);

  auto* h = new TH2D(hname,
                     title ? title : "Matrix;Column j;Row i",
                     nCols, -0.5, nCols - 0.5,
                     nRows, -0.5, nRows - 0.5);

  double minVal = 0.0, maxVal = 0.0;
  bool first = true;
  for (int i = 0; i < nRows; ++i) {
    for (int j = 0; j < nCols; ++j) {
      const double v = M(i, j);
      h->SetBinContent(j + 1, i + 1, v);  // X=j, Y=i
      if (first) { minVal = maxVal = v; first = false; }
      else { if (v < minVal) minVal = v; if (v > maxVal) maxVal = v; }
    }
  }


  gStyle->SetNumberContours(999);
  // h->SetContour(999);
  h->SetMinimum(1.);
  h->SetMaximum(1e10);
  c->SetLogz();
  // if (autoLog && minVal > 0.0) c->SetLogz();

  h->GetXaxis()->CenterTitle(true);
  h->GetYaxis()->CenterTitle(true);
  h->GetZaxis()->SetTitleOffset(1.1);

  h->Draw("COLZ");
  if (out && out[0]) c->SaveAs(out);
}



static void CompareInverses(const TMatrixD& cov, const char* tag = "cov")
{
  const int N = cov.GetNrows();
  if (N != cov.GetNcols()) {
    std::cerr << "[cmp] matrix not square\n";
    return;
  }

  // Symmetrize (numerical hygiene): covSym = 0.5*(cov + cov^T)
  TMatrixD covT(TMatrixD::kTransposed, cov);
  TMatrixD covSym = cov; covSym += covT; covSym *= 0.5;

  // avg diag (manual "trace")
  double sumdiag = 0.0;
  for (int i = 0; i < N; ++i) sumdiag += covSym(i,i);
  const double avg_var = (N > 0) ? sumdiag / N : 0.0;

  // SVD spectrum for conditioning and tolerance
  TDecompSVD svd_spec(covSym);
  if (!svd_spec.Decompose()) {
    std::cerr << "[cmp] SVD spectrum decomposition failed\n";
    return;
  }
  TVectorD s = svd_spec.GetSig();
  double smax = 0.0, sminpos = 1e300;
  for (int i = 0; i < s.GetNrows(); ++i) {
    smax = std::max(smax, s[i]);
    if (s[i] > 0.0) sminpos = std::min(sminpos, s[i]);
  }
  const double cond = (sminpos > 0.0) ? (smax / sminpos) : 1e300;
  const double tol  = smax * 1e-12;  // your current pseudo-inverse cutoff

  auto print_diag10 = [&](const char* name, const TMatrixD& A) {
    std::cout << name << " diag (first 10):\n";
    std::cout.setf(std::ios::scientific);
    std::cout << std::setprecision(6);
    for (int i = 0; i < std::min(10, N); ++i)
      std::cout << "  [" << i << "] = " << A(i,i) << "\n";
  };

  std::cout << "[cmp] tag=" << tag
            << " N=" << N
            << " avg_diag=" << avg_var
            << " cond~" << cond
            << " tol=" << tol << "\n";

  // --- Cholesky inverse on symmetric copy ---
  TMatrixDSym covs(N);
  for (int i = 0; i < N; ++i)
    for (int j = 0; j <= i; ++j)
      covs(i,j) = covSym(i,j);

  bool chol_ok = false;
  TMatrixD inv_chol(N,N); inv_chol.Zero();
  {
    TDecompChol chol(covs);
    if (chol.Decompose()) {
      TMatrixDSym invs = chol.Invert();   // symmetric inverse
      inv_chol = TMatrixD(invs);          // convert to TMatrixD
      chol_ok = true;
    } else {
      std::cout << "[cmp] Cholesky failed (not SPD)\n";
    }
  }

  // --- SVD inverse with NO truncation (should match Cholesky for SPD) ---
  TMatrixD inv_svd_full(N,N); inv_svd_full.Zero();
  {
    TDecompSVD svd(covSym); svd.Decompose();
    TMatrixD U = svd.GetU(), V = svd.GetV();
    TVectorD sf = svd.GetSig();
    TMatrixD Sinv(N,N); Sinv.Zero();
    for (int i = 0; i < N; ++i)
      Sinv(i,i) = (sf[i] != 0.0) ? 1.0 / sf[i] : 0.0;
    inv_svd_full = V * Sinv * TMatrixD(TMatrixD::kTransposed, U);
  }

  // --- SVD pseudo-inverse with tolerance (your current behavior) ---
  TMatrixD inv_svd_tol(N,N); inv_svd_tol.Zero();
  {
    TDecompSVD svd(covSym); svd.Decompose();
    TMatrixD U = svd.GetU(), V = svd.GetV();
    TVectorD sf = svd.GetSig();
    TMatrixD Sinv(N,N); Sinv.Zero();
    for (int i = 0; i < N; ++i)
      if (sf[i] > tol) Sinv(i,i) = 1.0 / sf[i];
    inv_svd_tol = V * Sinv * TMatrixD(TMatrixD::kTransposed, U);
  }

  if (chol_ok) print_diag10("[cmp] inv_chol",     inv_chol);
  print_diag10("[cmp] inv_svd_full", inv_svd_full);
  print_diag10("[cmp] inv_svd_tol ", inv_svd_tol);
}

// Manual 1-norm (max column sum of |.|)
static double OneNorm(const TMatrixD& A) {
  const int m = A.GetNrows(), n = A.GetNcols();
  double best = 0.0;
  for (int j = 0; j < n; ++j) {
    double col = 0.0;
    for (int i = 0; i < m; ++i) col += std::abs(A(i,j));
    best = std::max(best, col);
  }
  return best;
}

// Frobenius norm
static double FroNorm(const TMatrixD& A) {
  const int m = A.GetNrows(), n = A.GetNcols();
  long double s = 0.0L;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j) {
      const double v = A(i,j);
      s += (long double)v * (long double)v;
    }
  return std::sqrt((double)s);
}

// Max absolute element
static double MaxAbs(const TMatrixD& A) {
  const int m = A.GetNrows(), n = A.GetNcols();
  double best = 0.0;
  for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
      best = std::max(best, std::abs(A(i,j)));
  return best;
}

// Check how well `inv` inverts `cov`
void CheckInverseResidual(const TMatrixD& cov, const TMatrixD& inv, const char* tag = "inv") {
  const int N = cov.GetNrows();
  if (N != cov.GetNcols() || inv.GetNrows() != N || inv.GetNcols() != N) {
    std::cerr << "[resid] dimension mismatch\n";
    return;
  }

  // Symmetrize: Csym = 0.5*(C + C^T)
  TMatrixD covT(TMatrixD::kTransposed, cov);
  TMatrixD covSym = cov; covSym += covT; covSym *= 0.5;

  // Build identity
  TMatrixD I(TMatrixD::kUnit, covSym);

  // R1 = Csym * inv - I
  TMatrixD R1 = covSym; R1 *= inv; R1 -= I;

  // R2 = inv * Csym - I
  TMatrixD R2 = inv; R2 *= covSym; R2 -= I;

  // Some scale info (avg diag of Csym)
  double sumdiag = 0.0; for (int i = 0; i < N; ++i) sumdiag += covSym(i,i);
  const double avg_diag = (N > 0) ? sumdiag / N : 0.0;

  std::cout.setf(std::ios::scientific);
  std::cout << std::setprecision(6);
  std::cout << "[resid] tag=" << tag
            << "  N=" << N
            << "  avg_diag(C)=" << avg_diag << "\n";

  auto printNorms = [&](const char* name, const TMatrixD& R){
    std::cout << "  " << name
              << "  ||.||_F=" << FroNorm(R)
              << "  ||.||_1=" << OneNorm(R)
              << "  max|.|="  << MaxAbs(R) << "\n";
  };

  printNorms("R1 = C*inv - I", R1);
  printNorms("R2 = inv*C - I", R2);
}






TMatrixD TemplateFitter::BuildPredictionVector(const double *par, bool plotBestfit, const char* plotName) const {
  // 1. Clone base histograms
  TH1D* baseCC  = (TH1D*)CC[0].mCC_cc->Clone("baseCC");  baseCC->Reset();
  TH1D* baseNue = (TH1D*)nue[0].m_nue->Clone("baseNue"); baseNue->Reset();

  // ===== CC νμ (unchanged) =====
  TH1D* CCm_tp     = (TH1D*)baseCC->Clone("CCm_tp");
  TH1D* CCm_sig    = (TH1D*)baseCC->Clone("CCm_sig");
  TH1D* CCm_sig_mm = (TH1D*)baseCC->Clone("CCm_sig_mm");
  TH1D* CCm_sig_em = (TH1D*)baseCC->Clone("CCm_sig_em");
  // ===== CC νμ background (unchanged) =====
  TH1D* CCm_bkg    = (TH1D*)baseCC->Clone("CCm_bkg");
  TH1D* CCm_bkg_nc = (TH1D*)baseCC->Clone("CCm_bkg_mm");

  // ===== CC νe SIGNAL (from CC eCC, non-inv) =====
  TH1D* CCe_tp     = (TH1D*)baseCC->Clone("CCe_tp");
  TH1D* CCe_sig    = (TH1D*)baseCC->Clone("CCe_sig");
  TH1D* CCe_sig_ee = (TH1D*)baseCC->Clone("CCe_sig_ee");
  TH1D* CCe_sig_me = (TH1D*)baseCC->Clone("CCe_sig_me"); // off-diagonal (appearance)

  // ===== CC νe BACKGROUND (from ν–e that fail ν–e cut → non-inv) =====
  TH1D* CCe_bkg    = (TH1D*)baseCC->Clone("CCe_bkg");
  TH1D* CCe_bkg_nc = (TH1D*)baseCC->Clone("CCe_bkg_nc");
  TH1D* CCe_bkg_ee = (TH1D*)baseCC->Clone("CCe_bkg_ee");
  TH1D* CCe_bkg_em = (TH1D*)baseCC->Clone("CCe_bkg_em");
  TH1D* CCe_bkg_me = (TH1D*)baseCC->Clone("CCe_bkg_me");
  TH1D* CCe_bkg_mm = (TH1D*)baseCC->Clone("CCe_bkg_mm");

  // ===== ν–e SIGNAL (uses inv histos) =====
  TH1D* nue_tp     = (TH1D*)baseNue->Clone("nue_tp");
  TH1D* nue_sig    = (TH1D*)baseNue->Clone("nue_sig");
  TH1D* nue_sig_me = (TH1D*)baseNue->Clone("nue_sig_me");
  TH1D* nue_sig_mm = (TH1D*)baseNue->Clone("nue_sig_mm");
  TH1D* nue_sig_em = (TH1D*)baseNue->Clone("nue_sig_em");
  TH1D* nue_sig_ee = (TH1D*)baseNue->Clone("nue_sig_ee");

  // ===== ν–e BACKGROUND (CC eCC leaking into ν–e selection → inv) =====
  TH1D* nue_bkg    = (TH1D*)baseNue->Clone("nue_bkg");
  TH1D* nue_bkg_em = (TH1D*)baseNue->Clone("nue_bkg_em");
  TH1D* nue_bkg_ee = (TH1D*)baseNue->Clone("nue_bkg_ee");
  TH1D* nue_bkg_nc = (TH1D*)baseNue->Clone("nue_bkg_NC");
  // (Typically no ν–e background into me/mm from CC μ; leave those at 0 unless you intend some)

  // 2. Loop over energy bins and compute probabilities
  for (int i = 10; i < nbins_Ev; ++i) {
    int iL;
    if(nbins_Ev_LE == 29) {
      // THIS IS OLD MAPPING FOR 29 ENERGY BANDS
      iL = (i < 240) ? i / 10 : (i < 400 ? 24 + (i - 240) / 40 : 28); 
    }
    else if(nbins_Ev_LE == 12){
      // upper edges (exclusive) for iL = 0..11
      static const int iBound[12] = {20, 40, 60, 80, 120, 160, 202, 220, 240, 250, 280, 430};
      iL = std::upper_bound(iBound, iBound + 12, i) - iBound;
    }

    double ft_m[7], ft_e[7];
    for (int k = 0; k < 7; ++k) {
      ft_m[k] = LEfit.LEfit_m[iL][k];
      ft_e[k] = LEfit.LEfit_e[iL][k];
    }

    double me = 0, mt = 0, em = 0, et = 0, ee = 0, mm = 0;
    for (int j = 0; j <= 1000; ++j) {
      double e = binEdges[i] + j * (binEdges[i+1] - binEdges[i]) / 1000.0;
      me += getAvgPme(e, par[0], par[1], par[2], par[3], ft_m);
      mt += getAvgPmt(e, par[0], par[1], par[2], par[3], ft_m);
      em += getAvgPme(e, par[0], par[1], par[2], par[3], ft_e);
      et += getAvgPet(e, par[0], par[1], par[2], par[3], ft_e);
      ee += getAvgPee(e, par[0], par[1], par[2], par[3], ft_e);
      mm += getAvgPmm(e, par[0], par[1], par[2], par[3], ft_m);
    }

    double Pme = me / 1001.0;
    double Pem = em / 1001.0;
    double Pee = ee / 1001.0;
    double Pet = et / 1001.0;
    double Pmm = mm / 1001.0;
    double Pmt = mt / 1001.0;

    // double Pm_act = Pem + Pmm + Pmt;

    // Building vmu-CC
    CCm_sig_mm->Add(CC[i].mCC_cc, Pmm);                 // νμ CC survival
    CCm_bkg_nc->Add(CC[i].mCC_nc, Pme+Pmm+Pmt);         // NC survival (mainly νμ -> weight by mixing between νμ with active flavors)
    CCm_sig_em->Add(CC[i].eCC_cc, Pem);                 // νe→νμ appearance

    // Building ve-CC
    // CCe_sig_me->Add(CC[i].mCC_NoCut, Pme);
    // CCe_sig_ee->Add(CC[i].eCC, Pee);
    // CCνe SIGNAL (from CC eCC, non-inv)
    CCe_sig_ee->Add(CC[i].eCC_cc, Pee);                 // νe CC survival
    CCe_bkg_nc->Add(CC[i].eCC_nc, Pme+Pmm+Pmt);         // NC survival (mainly νμ -> weight by mixing between νμ with active flavors)
    CCe_sig_me->Add(CC[i].mCC_w, Pme);                  // νμ→νe appearance (using no reco cut true CCνμ)

    // CCνe BACKGROUND (from ν–e failing the ν–e cut → non-inv)
    // Tag by oscillation channel (α→β for the neutrino), but ALL of these stay in the CCνe sample.
    CCe_bkg_mm->Add(CC[i].eCC_BKGm,    Pmm + Pmt);  // νμ(±ντ) survival    → CCνe bkg
    CCe_bkg_me->Add(CC[i].eCC_BKGm_w,  Pme);        // νμ→νe (appearance) → CCνe bkg
    CCe_bkg_ee->Add(CC[i].eCC_BKGe,    Pee);        // νe→νe (survival)   → CCνe bkg
    CCe_bkg_em->Add(CC[i].eCC_BKGe_w,  Pem + Pet);        // νe→νμ (disappearance) → CCνe bkg 


    // Building v-e
    // v-e SIGNAL
    nue_sig_mm->Add(nue[i].m_nue,   Pmm + Pmt);
    nue_sig_me->Add(nue[i].m_nue_w, Pme);
    nue_sig_ee->Add(nue[i].e_nue,   Pee);
    nue_sig_em->Add(nue[i].e_nue_w, Pem + Pet);

    // v-e BACKGROUND (from ve-CC failing the inverse Etheta2 cut)
    nue_bkg_ee->Add(nue[i].nue_BKGe_cc, Pee);
    nue_bkg_em->Add(nue[i].nue_BKGe_cc, Pem);
    nue_bkg_nc->Add(nue[i].nue_BKGe_nc, Pme+Pmm+Pmt); // NC survival (mainly νμ -> weight by mixing between νμ with active flavors)
  }

  CCm_sig->Add(CCm_sig_em); CCm_sig->Add(CCm_sig_mm);
  CCm_bkg->Add(CCm_bkg_nc);
  CCm_tp->Add(CCm_sig); CCm_tp->Add(CCm_bkg);

  CCe_sig->Add(CCe_sig_me); CCe_sig->Add(CCe_sig_ee);
  CCe_bkg->Add(CCe_bkg_nc); CCe_bkg->Add(CCe_bkg_me); CCe_bkg->Add(CCe_bkg_ee); CCe_bkg->Add(CCe_bkg_em); CCe_bkg->Add(CCe_bkg_mm);
  CCe_tp->Add(CCe_sig); CCe_tp->Add(CCe_bkg); 

  nue_sig->Add(nue_sig_me); nue_sig->Add(nue_sig_ee); nue_sig->Add(nue_sig_em); nue_sig->Add(nue_sig_mm);
  nue_bkg->Add(nue_bkg_nc); nue_bkg->Add(nue_bkg_ee); nue_bkg->Add(nue_bkg_em);
  nue_tp->Add(nue_sig); nue_tp->Add(nue_bkg); 


  // // Sum the parent pool feeding the first νeCC reco bin through mCC_w
  // double mCCw_bin1_sum = 0.0;
  // for (int ii = 10; ii < nbins_Ev; ++ii) {
  //   mCCw_bin1_sum += CC[ii].mCC_w->GetBinContent(1);
  // }

  // std::cout << plotName
  //           << " mCC_w parent sum(bin1) = " << mCCw_bin1_sum
  //           << " ; CCe_sig_me(bin1) = " << CCe_sig_me->GetBinContent(1)
  //           << " ; implied <Pme> = " << (CCe_sig_me->GetBinContent(1) / mCCw_bin1_sum)
  //           << " ; theory max = " << (4.0 * par[0] * par[1])
  //           << "\n";

  // hRatioToNull->Reset();
  // hRatioToNull->Reset();
  // 3. Fill TMatrixD
  TMatrixD prediction(nbins, 1);
  for (int bx = 0; bx < nbins; ++bx) {
    if (bx < nbinsCC)
      prediction[bx][0] = CCm_tp->GetBinContent(bx + 1);
    else if (bx < 2 * nbinsCC)
      prediction[bx][0] = CCe_tp->GetBinContent(bx - nbinsCC + 1);
    else
      prediction[bx][0] = nue_tp->GetBinContent(bx - 2 * nbinsCC + 1);

    hRatioToNull->SetBinContent(bx+1, prediction[bx][0]);
  }

  if(par[0]==0. && par[1]==0. && par[2]==0. && par[3]==0.){
    for(int i=0; i<nbins; i++){
      hNull->SetBinContent(i+1, prediction[i][0]);
    }
  }

  if(hNull && !IsEmpty(hNull)) hRatioToNull->Divide(hNull);

  // std::cout << "Original ratio bin 55 = "
  //           << hRatioToNull->GetBinContent(55) << "\n";
  // std::cout << "Original ratio bin 56 = "
  //           << hRatioToNull->GetBinContent(56) << "\n";

  if ( plotBestfit ) {
    //std::cout << plotName << std::endl;
    TCanvas* c = new TCanvas("c", "", 800, 700);
    plotOscillationComponents( c, par, CCm_sig_em,CCm_sig_mm,CCm_bkg_nc, CCe_sig_me,CCe_sig_ee, CCe_bkg_nc,CCe_bkg_me,CCe_bkg_ee,CCe_bkg_em,CCe_bkg_mm, nue_sig_mm,nue_sig_me,nue_sig_ee,nue_sig_em, nue_bkg_ee,nue_bkg_em, plotBestfit );
    c->SaveAs(Form("%s.png",plotName));
    c->Close();
  }
  // else if(plotBestfit==false){
  //   if(plotTarget==true){
  //     TCanvas* c = new TCanvas("c", "", 900, 600);
  //     plotOscillationComponents( c, par, CC_tp_mm, CC_tp_me, CC_tp_ee, CC_tp_em, nue_tp_mm, nue_tp_me, nue_tp_ee, nue_tp_em );
  //     c->SaveAs(Form("%s.png",plotName));
  //     c->Close();
  //   }
  // }

  return prediction;
}


double TemplateFitter::getTarget(double par[4], int nucut, bool originalTgt) {
  // 1) Build targets (full size)
  target = BuildPredictionVector(par);    // N×1
  kNuCut = nucut;
  // FCtarget = target;
  double area = 0., throw_area = 0.;
  for (int i = 0; i < nbins; ++i) {
    FCtarget(i,0) = target(i,0) * (1.0 + cov.FCWeights[i]);
    area       += target(i, 0);
    throw_area += FCtarget(i, 0);
  }
  const TMatrixD& tgtFull = originalTgt ? target : FCtarget;

  // 2) Stat (full)
  statmx.Zero();
  for (int i = 0; i < nbins; ++i) statmx(i,i) = tgtFull(i,0);

  // 3) Absolute systematics (build locally to avoid double-scaling): Use MC-predicted statistics to scale fractional systematic covar matrix
  TMatrixD sys_abs(nbins, nbins);
  for (int i = 0; i < nbins; ++i) {
    for (int j = 0; j < nbins; ++j)
      sys_abs(i,j) = /* fractional */ sysmx(i,j) * target(i,0) * target(j,0);
    // if(i<10) std::cout << i << "\t" << statmx(i,i) << "\t" << sys_abs(i,i) << "\n";
  }

  // 4) Total absolute covariance matrix
  covmx = statmx + sys_abs;
  // DrawMatrix(covmx, "Full covariance;Column j;Row i", "cov_full.png");

  // 5) Apply drop only to the cached members we need later
  const TMatrixD tgt_final = dropBins_.empty()
      ? tgtFull
      : ReduceMatrix(tgtFull, dropBins_);
  // std::cout << tgt_final.GetNrows() << "\t" << tgt_final.GetNcols() << "\n";
  // DrawMatrix(tgt_final, "Target;Bin index;Counts", "target.png");


  TMatrixD cov = dropBins_.empty()
      ? covmx
      : ReduceMatrix(covmx, dropBins_);
  // DrawMatrix(cov, "Reduced covariance;Column j;Row i", "cov_reduced.png");

  // CompareInverses(cov, "after_drop");

  myTarget.ResizeTo(tgt_final);
  myTarget = tgt_final;
  // // Print myTarget (assumed column vector)
  // std::cout << "myTarget (" << myTarget.GetNrows() << "x" << myTarget.GetNcols() << ")\n";
  // for (int i = 0; i < myTarget.GetNrows(); ++i) {
  //   std::cout << "  [" << i << "] = " << myTarget(i,0) << "\n";
  // }

  const int Nk = myTarget.GetNrows();

  // --- Symmetrize numerically: covSymD = 0.5*(cov + cov^T) ---
  TMatrixD covT(TMatrixD::kTransposed, cov);
  TMatrixD covSymD = cov; covSymD += covT; covSymD *= 0.5;

  // Average diag for ridge scale
  double avg_diag = 0.0;
  for (int i = 0; i < Nk; ++i) avg_diag += covSymD(i,i);
  avg_diag = (Nk > 0) ? avg_diag / Nk : 0.0;

  // Helper: ||C*inv - I||_1
  auto resid1 = [](const TMatrixD& C, const TMatrixD& inv) {
    TMatrixD I(TMatrixD::kUnit, C);
    TMatrixD R = C; R *= inv; R -= I;
    return R.Norm1();
  };

  bool chol_ok = false;
  double lambda = 0.0;                 // ridge
  const double target_resid = 1e-6;    // accept threshold (tune if needed)

  for (int tries = 0; tries < 6 && !chol_ok; ++tries) {
    if (tries > 0) lambda = (lambda == 0.0) ? (1e-12 * avg_diag) : (lambda * 10.0);

    // Build symmetric matrix with ridge: covSymD_reg = covSymD + λI
    TMatrixD covSymD_reg = covSymD;
    for (int i = 0; i < Nk; ++i) covSymD_reg(i,i) += lambda;

    // Pack into TMatrixDSym
    TMatrixDSym covs(Nk);
    for (int i = 0; i < Nk; ++i)
      for (int j = 0; j <= i; ++j)
        covs(i,j) = covSymD_reg(i,j);

    TDecompChol chol(covs);
    if (!chol.Decompose()) continue;

    // Invert via Cholesky
    TMatrixDSym invs = chol.Invert();
    TMatrixD inv = TMatrixD(invs);

    // Accept only if it actually inverts the ridged matrix
    if (resid1(covSymD_reg, inv) < target_resid) {
      invmx.ResizeTo(Nk, Nk);
      invmx = inv;
      // Optional: check residual & print
      CheckInverseResidual(covSymD_reg, invmx, "chol+ridge");
      chol_ok = true;
    }
  }

  // Fallback: SVD pseudo-inverse (same symmetrized input)
  if (!chol_ok) {
    TDecompSVD svd(covSymD);
    if (!svd.Decompose())
      throw std::runtime_error("getTarget: SVD decomposition failed");

    TVectorD s = svd.GetSig();
    double smax = 0.0; for (int i = 0; i < s.GetNrows(); ++i) smax = std::max(smax, s[i]);
    const double tol = smax * 1e-12;   // scale-aware cutoff

    TMatrixD U = svd.GetU(), V = svd.GetV();
    TMatrixD Sinv(Nk, Nk); Sinv.Zero();
    for (int i = 0; i < Nk; ++i) if (s[i] > tol) Sinv(i,i) = 1.0 / s[i];

    invmx.ResizeTo(Nk, Nk);
    invmx = V * Sinv * TMatrixD(TMatrixD::kTransposed, U);

    Optional: CheckInverseResidual(covSymD, invmx, "svd");
  }

  // // SVD pseudo-inverse (scale-aware tol)
  // TDecompSVD svd(cov);
  // if (!svd.Decompose()) throw std::runtime_error("getTarget: SVD decomposition failed");
  // TVectorD s = svd.GetSig();
  // double smax = 0.0; for (int i = 0; i < s.GetNrows(); ++i) smax = std::max(smax, s[i]);
  // const double tol = smax * 1e-12;

  // TMatrixD U = svd.GetU();
  // TMatrixD V = svd.GetV();
  // TMatrixD Sinv(Nk, Nk); Sinv.Zero();
  // for (int i = 0; i < Nk; ++i) if (s[i] > tol) Sinv(i,i) = 1.0 / s[i];

  // invmx.ResizeTo(Nk, Nk);
  // invmx = V * Sinv * TMatrixD(TMatrixD::kTransposed, U);
  // CheckInverseResidual(cov, invmx, "svd");    // for the SVD-built inverse
  return (area == 0.) ? 0. : (throw_area/area);
  // return 0;
}

double TemplateFitter::CalcChi2Core(const TMatrixD& prediction) const {
  const TMatrixD pred = dropBins_.empty()
      ? prediction
      : ReduceMatrix(prediction, dropBins_);
  // std::cout << pred.GetNrows() << "\t" << pred.GetNcols() << "\n";
  // DrawMatrix(pred, "Prediction;Bin index;Counts", "pred.png");

  TMatrixD diff = pred - myTarget;
  TMatrixD diffT(TMatrixD::kTransposed, diff);
  TMatrixD chi2 = diffT * invmx * diff;
  return chi2(0,0);
}

double TemplateFitter::CalculateChi2(double par[4], bool plotBestfit, const char* outName) {
  TMatrixD prediction = BuildPredictionVector(par, plotBestfit, outName);
  // return 0.;
  return CalcChi2Core(prediction);
}

double TemplateFitter::getChi2(const double* par) {
  TMatrixD prediction = BuildPredictionVector(par);
  double chi2_val = CalcChi2Core(prediction);

  double penalty = 0.0;

  // // Penalty on extreme dm2
  // double epsilon = 1E-8;
  double threshold = 0.05;
  // if (par[3] < threshold) penalty = (threshold / (par[3] + epsilon)) * (threshold / (par[3] + epsilon));
  if (par[3] < threshold) {
    double log_dm2 = std::log10(par[3] + 1e-8);
    double log_center = std::log10(threshold);
    double sigma = 0.3;
    penalty = std::pow((log_dm2 - log_center) / sigma, 2);
  }
  else if (par[3] > 500.) penalty = (par[3] - 500.) * (par[3] - 500.);

  // Penalty on unphysical parameter regions
  double lambda = 1E6;
  double sumU = par[0] + par[1] + par[2];
  if (sumU > 1.0) {
    double excess = sumU - 1.0;
    penalty += lambda * excess * excess;  // strong quadratic penalty
  }

  // std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << par[3] << "\t" << chi2_val << "\t" << chi2_val + penalty << "\n";

  return chi2_val + penalty;
}

bool TemplateFitter::doFit(int nFunctionCalls, int nIterations, double Tolerance, const char* varName[4], double seed[4], double stepSize[4], std::vector<double>& parBf, std::vector<double>& parError, double& chi2, bool fixedUt42) {
  ROOT::Math::Minimizer* fitter = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  fitter->SetMaxFunctionCalls(nFunctionCalls);
  fitter->SetMaxIterations(nIterations);
  fitter->SetTolerance(Tolerance);

  for (int i = 0; i < 4; ++i) {
    fitter->SetVariable(i, varName[i], seed[i], stepSize[i]);
    fitter->SetVariableLimits(i, par_min[i], par_max[i]);
  }
  if (fixedUt42) {
    fitter->SetVariableValue(2, 0.0);  // set value explicitly
    fitter->FixVariable(2);            // freeze at this value
  } 

  ROOT::Math::Functor functor(this, &TemplateFitter::getChi2, 4);
  fitter->SetFunction(functor);
  fitter->Minimize();

  const double* bestfit = fitter->X();
  const double* bestfitError = fitter->Errors();

  parBf.resize(4);
  parError.resize(4);

  for (int i = 0; i < 4; ++i) {
    parBf[i] = bestfit[i];
    parError[i] = bestfitError[i];
  }

  chi2 = fitter->MinValue();
  return true;
}



// TMatrixD TemplateFitter::BuildPredictionVector(const double *par, bool plotBestfit, const char* plotName) const {
//   // 1. Clone base histograms
//   TH1D* baseCC = (TH1D*)CC[0].mCC->Clone("baseCC"); baseCC->Reset();
//   TH1D* CC_tp_m  = (TH1D*)baseCC->Clone("CC_tp_m");
//   TH1D* CC_tp_e  = (TH1D*)baseCC->Clone("CC_tp_e");
//   TH1D* CC_tp_me = (TH1D*)baseCC->Clone("CC_tp_me");
//   TH1D* CC_tp_ee = (TH1D*)baseCC->Clone("CC_tp_ee");
//   TH1D* CC_tp_em = (TH1D*)baseCC->Clone("CC_tp_em");
//   TH1D* CC_tp_mm = (TH1D*)baseCC->Clone("CC_tp_mm");

//   TH1D* baseNue = (TH1D*)nue[0].m_nue->Clone("baseNue"); baseNue->Reset();
//   TH1D* nue_tp       = (TH1D*)baseNue->Clone("nue_tp");
//   TH1D* nue_tp_me    = (TH1D*)baseNue->Clone("nue_tp_me");
//   TH1D* nue_tp_em    = (TH1D*)baseNue->Clone("nue_tp_em");
//   TH1D* nue_tp_ee    = (TH1D*)baseNue->Clone("nue_tp_ee");
//   TH1D* nue_tp_mm    = (TH1D*)baseNue->Clone("nue_tp_mm");
//   TH1D* nue_tp_os    = (TH1D*)baseNue->Clone("nue_tp_os");
//   TH1D* nue_tp_unos  = (TH1D*)baseNue->Clone("nue_tp_unos");

//   // 2. Loop over energy bins and compute probabilities
//   for (int i = 10; i < nbins_Ev; ++i) {
//     int iL = (i < 240) ? i / 10 : (i < 400 ? 24 + (i - 240) / 40 : 28);
//     double ft_m[7], ft_e[7];
//     for (int k = 0; k < 7; ++k) {
//       ft_m[k] = LEfit.LEfit_m[iL][k];
//       ft_e[k] = LEfit.LEfit_e[iL][k];
//     }

//     double me = 0, mt = 0, em = 0, ee = 0, mm = 0;
//     for (int j = 0; j <= 1000; ++j) {
//       double e = binEdges[i] + j * (binEdges[i+1] - binEdges[i]) / 1000.0;
//       me += getAvgPme(e, par[0], par[1], par[2], par[3], ft_m);
//       mt += getAvgPmt(e, par[0], par[1], par[2], par[3], ft_m);
//       em += getAvgPme(e, par[0], par[1], par[2], par[3], ft_e);
//       ee += getAvgPee(e, par[0], par[1], par[2], par[3], ft_e);
//       mm += getAvgPmm(e, par[0], par[1], par[2], par[3], ft_m);
//     }

//     double Pme = me / 1001.0;
//     double Pem = em / 1001.0;
//     double Pee = ee / 1001.0;
//     double Pmm = mm / 1001.0;
//     double Pmt = mt / 1001.0;

//     CC_tp_me->Add(CC[i].mCC_NoCut, Pme);
//     CC_tp_ee->Add(CC[i].eCC, Pee);
//     CC_tp_em->Add(CC[i].eCC, Pem);
//     CC_tp_mm->Add(CC[i].mCC, Pmm);

//     nue_tp_me->Add(nue[i].m_nue_w, Pme);
//     nue_tp_mm->Add(nue[i].m_nue, Pmm + Pmt);
//     nue_tp_em->Add(nue[i].e_nue_w, Pem);
//     nue_tp_ee->Add(nue[i].e_nue, Pee);
//   }

//   for (int i = 41; i < nbinsCC; ++i) {
//     CC_tp_mm->SetBinContent(i + 1, 0.0);
//     CC_tp_em->SetBinContent(i + 1, 0.0);
//   }

//   CC_tp_e->Add(CC_tp_me); CC_tp_e->Add(CC_tp_ee);
//   CC_tp_m->Add(CC_tp_em); CC_tp_m->Add(CC_tp_mm);
//   nue_tp_os->Add(nue_tp_me); nue_tp_os->Add(nue_tp_em);
//   nue_tp_unos->Add(nue_tp_ee); nue_tp_unos->Add(nue_tp_mm);
//   nue_tp->Add(nue_tp_os); nue_tp->Add(nue_tp_unos);

//   // 3. Fill TMatrixD
//   TMatrixD prediction(nbins, 1);
//   for (int bx = 0; bx < nbins; ++bx) {
//     if (bx < nbinsCC)
//       prediction[bx][0] = CC_tp_m->GetBinContent(bx + 1);
//     else if (bx < 2 * nbinsCC)
//       prediction[bx][0] = CC_tp_e->GetBinContent(bx - nbinsCC + 1);
//     else
//       prediction[bx][0] = nue_tp->GetBinContent(bx - 2 * nbinsCC + 1);
//   }

//   if ( plotBestfit ) {
//     //std::cout << plotName << std::endl;
//     TCanvas* c = new TCanvas("c", "", 900, 600);
//     plotOscillationComponents( c, par, CC_tp_mm, CC_tp_me, CC_tp_ee, CC_tp_em, nue_tp_mm, nue_tp_me, nue_tp_ee, nue_tp_em, plotBestfit );
//     c->SaveAs(Form("%s.png",plotName));
//     c->Close();
//   }
//   // else if(plotBestfit==false){
//   //   if(plotTarget==true){
//   //     TCanvas* c = new TCanvas("c", "", 900, 600);
//   //     plotOscillationComponents( c, par, CC_tp_mm, CC_tp_me, CC_tp_ee, CC_tp_em, nue_tp_mm, nue_tp_me, nue_tp_ee, nue_tp_em );
//   //     c->SaveAs(Form("%s.png",plotName));
//   //     c->Close();
//   //   }
//   // }

//   //std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << par[3] << "\n";

//   return prediction;
// }



// double TemplateFitter::getTarget(double par[4], bool originalTgt) {
//   target = BuildPredictionVector(par);


//   for(int i=0;i<nbins;i++){
//     for(int j=0;j<nbins;j++){
//       sysmx(i,j) = sysmx(i,j) * target(i,0) * target(j,0);
//     }
//     // statmx(i, i) = target(i, 0);
//   }

//   FCtarget = target;
//   double area = 0., throw_area = 0.;

//   for (int i = 0; i < nbins; ++i) {
//     FCtarget(i, 0) *= (1.0 + cov.FCWeights[i]);
//     statmx(i, i) = FCtarget(i, 0);
//     area       += target(i, 0);
//     throw_area += FCtarget(i, 0);
//   }

//   if (originalTgt) myTarget = target;
//   else myTarget = FCtarget;

//   covmx = statmx + sysmx;
//   TDecompSVD svd(covmx);
//   if (!svd.Decompose()) {
//     std::cerr << "SVD Decomposition failed!" << std::endl;
//   }

//   TVectorD singularValues = svd.GetSig();
//   TMatrixD U = svd.GetU();
//   TMatrixD V = svd.GetV();
//   TMatrixD Ut(TMatrixD::kTransposed, U);

//   TMatrixD Sinv(nbins, nbins);
//   Sinv.Zero();
//   for (int i = 0; i < nbins; ++i) {
//     if (singularValues[i] > 1e-12)
//       Sinv(i, i) = 1.0 / singularValues[i];
//   }

//   invmx = V * Sinv * Ut;
//   return throw_area/area;
// }

// double TemplateFitter::CalcChi2Core(const TMatrixD& prediction) const {
//   TMatrixD diff = prediction - myTarget;
//   TMatrixD diffT(TMatrixD::kTransposed, diff);
//   TMatrixD chi2 = diffT * invmx * diff;
//   return chi2[0][0];
// }

// double TemplateFitter::CalculateChi2(double par[4], bool plotBestfit, const char* outName) {
//   //std::cout << outName << std::endl;
//   TMatrixD prediction = BuildPredictionVector(par, plotBestfit, outName);
//   return CalcChi2Core(prediction);
// }

// double TemplateFitter::getChi2(const double* par) {
//   TMatrixD prediction = BuildPredictionVector(par);
//   double chi2_val = CalcChi2Core(prediction);

//   double penalty = 0.0;

//   // // Penalty on extreme dm2
//   // double epsilon = 1E-8;
//   double threshold = 0.05;
//   // if (par[3] < threshold) penalty = (threshold / (par[3] + epsilon)) * (threshold / (par[3] + epsilon));
//   if (par[3] < threshold) {
//     double log_dm2 = std::log10(par[3] + 1e-8);
//     double log_center = std::log10(threshold);
//     double sigma = 0.3;
//     penalty = std::pow((log_dm2 - log_center) / sigma, 2);
//   }
//   else if (par[3] > 500.) penalty = (par[3] - 500.) * (par[3] - 500.);

//   // Penalty on unphysical parameter regions
//   double lambda = 1E6;
//   double sumU = par[0] + par[1] + par[2];
//   if (sumU > 1.0) {
//     double excess = sumU - 1.0;
//     penalty += lambda * excess * excess;  // strong quadratic penalty
//   }

//   // std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << par[3] << "\t" << chi2_val << "\t" << chi2_val + penalty << "\n";

//   return chi2_val + penalty;
// }

// bool TemplateFitter::doFit(int nFunctionCalls, int nIterations, double Tolerance, const char* varName[4], double seed[4], double stepSize[4], std::vector<double>& parBf, std::vector<double>& parError, double& chi2, bool fixedUt42) {
//   ROOT::Math::Minimizer* fitter = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
//   fitter->SetMaxFunctionCalls(nFunctionCalls);
//   fitter->SetMaxIterations(nIterations);
//   fitter->SetTolerance(Tolerance);

//   for (int i = 0; i < 4; ++i) {
//     fitter->SetVariable(i, varName[i], seed[i], stepSize[i]);
//     fitter->SetVariableLimits(i, par_min[i], par_max[i]);
//   }
//   if (fixedUt42) fitter->FixVariable(2);

//   ROOT::Math::Functor functor(this, &TemplateFitter::getChi2, 4);
//   fitter->SetFunction(functor);
//   fitter->Minimize();

//   const double* bestfit = fitter->X();
//   const double* bestfitError = fitter->Errors();

//   parBf.resize(4);
//   parError.resize(4);

//   for (int i = 0; i < 4; ++i) {
//     parBf[i] = bestfit[i];
//     parError[i] = bestfitError[i];
//   }

//   chi2 = fitter->MinValue();
//   return true;
// }