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
#include <algorithm> 
#include "TMatrixDSym.h"
#include "TVectorD.h"
#include <iostream>
#include <iomanip>
using namespace std;

// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}

// Initialize TemplateFitter constructor
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



void TemplateFitter::plotOscillationComponents(TCanvas* c, const double par[4],
                                               TH1D* CC_tp_mm, TH1D* CC_tp_me, TH1D* CC_tp_ee, TH1D* CC_tp_em,
                                               TH1D* nue_tp_mm, TH1D* nue_tp_me, TH1D* nue_tp_ee, TH1D* nue_tp_em,
                                               bool plotBestfit) const
{
  TH1D* tgt = new TH1D("tgt", "", nbins, 0, nbins);
  TH1D* FCtgt = new TH1D("FCtgt", "", nbins, 0, nbins);
  tgt->SetDirectory(nullptr);
  FCtgt->SetDirectory(nullptr);
  for (int i = 0; i < nbins; ++i) {
      tgt->SetBinContent(i + 1, this->target(i, 0));
      FCtgt->SetBinContent(i + 1, this->FCtarget(i, 0));
  }


  TH1D* h_mm = new TH1D("h_mm", "", nbins, 0, nbins);
  TH1D* h_me = new TH1D("h_me", "", nbins, 0, nbins);
  TH1D* h_ee = new TH1D("h_ee", "", nbins, 0, nbins);
  TH1D* h_em = new TH1D("h_em", "", nbins, 0, nbins);

  h_mm->SetDirectory(nullptr);
  h_me->SetDirectory(nullptr);
  h_ee->SetDirectory(nullptr);
  h_em->SetDirectory(nullptr);

  // Fill histograms
  for (int i = 0; i < nbins; ++i) {
    if (i < nbinsCC) {
      h_mm->SetBinContent(i + 1, CC_tp_mm->GetBinContent(i + 1));
      h_em->SetBinContent(i + 1, CC_tp_em->GetBinContent(i + 1));
    }
    else if (i >= nbinsCC && i < 2 * nbinsCC) {
      int idx = i - nbinsCC;
      h_ee->SetBinContent(i + 1, CC_tp_ee->GetBinContent(idx + 1));
      h_me->SetBinContent(i + 1, CC_tp_me->GetBinContent(idx + 1));
    }
    else {
      int idx = i - 2 * nbinsCC;
      h_ee->SetBinContent(i + 1, nue_tp_ee->GetBinContent(idx + 1));
      h_em->SetBinContent(i + 1, nue_tp_em->GetBinContent(idx + 1));
      h_mm->SetBinContent(i + 1, nue_tp_mm->GetBinContent(idx + 1));
      h_me->SetBinContent(i + 1, nue_tp_me->GetBinContent(idx + 1));
    }
  }

  // Assign colors: same origin → same color
  auto color_mu = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 1);
  auto color_e  = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 2);

  h_mm->SetFillColor(color_mu);  // mu-origin survival
  h_me->SetFillColor(color_mu);  // mu-origin oscillated
  h_ee->SetFillColor(color_e);   // e-origin survival
  h_em->SetFillColor(color_e);   // e-origin oscillated

  // Create the stack
  THStack* h_stack = new THStack("h_stack", 
      Form("Oscillation Hypothesis: U_{e4}^{2} = %.2f, U_{#mu4}^{2} = %.2f, U_{#tau4}^{2} = %.2f, #Delta m^{2} = %.1f", 
      par[0], par[1], par[2], par[3]));

  // Add to stack in correct order: survival first, oscillated on top
  h_stack->Add(h_mm);  // mu→mu survival
  h_stack->Add(h_em);  // e→mu oscillated
  h_stack->Add(h_ee);  // e→e survival
  h_stack->Add(h_me);  // mu→e oscillated

  h_stack->SetMinimum(80.);
  h_stack->SetMaximum(h_stack->GetMaximum()*8.);
  
  double x1 = h_mm->GetXaxis()->GetBinUpEdge(nbinsCC);
  double x2 = h_mm->GetXaxis()->GetBinUpEdge(nbinsCC+nbinsCC);

  // Draw
  c->cd();
  c->Clear();
  c->SetLogy();
  h_stack->Draw("hist");
  c->Update();
  TLatex latex;
  latex.SetNDC(false);
  //latex.SetTextColor(kRed);
  latex.SetTextSize(0.05); 	// Set text size
  latex.SetTextAlign(22);	// Center alignment
  double yPosition = pow(10, (gPad->GetUymin() + gPad->GetUymax())/2.); 	// Align to the middle of the pad
  latex.DrawLatex(nbinsCC/2, yPosition, "#nu_{#mu}-CC"); 	// Position for CCm
  latex.DrawLatex(nbinsCC+nbinsCC/2, yPosition, "#nu_{e}-CC"); 	// Position for CCe
  latex.DrawLatex(nbinsCC+nbinsCC + nbins_nue/2, yPosition, "#nu+e");	// Position for nue
  h_stack->GetXaxis()->SetTitle("bin number");
  h_stack->GetYaxis()->SetTitle("entries (/yr.POT)");
  TLine *line1 = new TLine(x1, 0, x1, pow(10, gPad->GetUymax()));
  TLine *line2 = new TLine(x2, 0, x2, pow(10, gPad->GetUymax()));
  line1->SetLineColor(kBlack);
  line2->SetLineColor(kBlack);
  line1->SetLineWidth(2.0);
  line2->SetLineWidth(2.0);
  line1->Draw("same");
  line2->Draw("same");
  //TLegend *lg = new TLegend(0.65,0.75,0.9,0.9);
  TLegend *leg = MakeLegend(0.6, 0.65, 0.75, 0.8);
  leg->AddEntry(h_mm,"#nu_{#mu} originated");
  leg->AddEntry(h_ee,"#nu_{e} originated");
  
  if(plotBestfit){
    tgt->SetLineColor(dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 3));
    FCtgt->SetLineColor(dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, 4));

    tgt->SetLineWidth(2.);
    tgt->SetLineStyle(9);

    FCtgt->SetLineWidth(2.);
    FCtgt->SetLineStyle(9);

    // tgt->Draw("HIST SAME");
    // FCtgt->Draw("HIST SAME");

    // leg->AddEntry(tgt, "original target");
    // leg->AddEntry(FCtgt, "FC thrown target");
  }

  leg->Draw();
  dunestyle::CenterTitles(h_stack->GetHistogram());
  dunestyle::Simulation();
  
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


// Drawing matrix
void DrawMatrix(const TMatrixD& M,
                       const char* title = "Matrix;Column j;Row i",
                       const char* out   = nullptr,   // e.g. "matrix.png" (nullptr = don't save)
                       bool autoLogZ     = true)
{
  const int nRows = M.GetNrows();
  const int nCols = M.GetNcols();

  static int uid = 0;
  const TString cname = Form("c_matrix_%d", uid);
  const TString hname = Form("h_matrix_%d", uid);
  ++uid;

  gStyle->SetOptStat(0);

  auto* c = new TCanvas(cname, cname, 950, 800);
  c->SetRightMargin(0.15);

  auto* h = new TH2D(hname, title, nCols, -0.5, nCols - 0.5,
                                   nRows, -0.5, nRows - 0.5);

  double minVal = 0.0, maxVal = 0.0;
  bool first = true;
  for (int i = 0; i < nRows; ++i) {
    for (int j = 0; j < nCols; ++j) {
      const double v = M(i, j);
      h->SetBinContent(j + 1, i + 1, v);
      if (first) { minVal = maxVal = v; first = false; }
      else { if (v < minVal) minVal = v; if (v > maxVal) maxVal = v; }
    }
  }

  h->SetContour(255);
  if (autoLogZ && minVal > 0.0) c->SetLogz();

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
  TH1D* baseCC = (TH1D*)CC[0].mCC->Clone("baseCC"); baseCC->Reset();
  TH1D* CC_tp_m  = (TH1D*)baseCC->Clone("CC_tp_m");
  TH1D* CC_tp_e  = (TH1D*)baseCC->Clone("CC_tp_e");
  TH1D* CC_tp_me = (TH1D*)baseCC->Clone("CC_tp_me");
  TH1D* CC_tp_ee = (TH1D*)baseCC->Clone("CC_tp_ee");
  TH1D* CC_tp_em = (TH1D*)baseCC->Clone("CC_tp_em");
  TH1D* CC_tp_mm = (TH1D*)baseCC->Clone("CC_tp_mm");

  TH1D* baseNue = (TH1D*)nue[0].m_nue->Clone("baseNue"); baseNue->Reset();
  TH1D* nue_tp       = (TH1D*)baseNue->Clone("nue_tp");
  TH1D* nue_tp_me    = (TH1D*)baseNue->Clone("nue_tp_me");
  TH1D* nue_tp_em    = (TH1D*)baseNue->Clone("nue_tp_em");
  TH1D* nue_tp_ee    = (TH1D*)baseNue->Clone("nue_tp_ee");
  TH1D* nue_tp_mm    = (TH1D*)baseNue->Clone("nue_tp_mm");
  TH1D* nue_tp_os    = (TH1D*)baseNue->Clone("nue_tp_os");
  TH1D* nue_tp_unos  = (TH1D*)baseNue->Clone("nue_tp_unos");

  // 2. Loop over energy bins and compute probabilities
  for (int i = 10; i < nbins_Ev; ++i) {
    int iL = (i < 240) ? i / 10 : (i < 400 ? 24 + (i - 240) / 40 : 28);
    double ft_m[7], ft_e[7];
    for (int k = 0; k < 7; ++k) {
      ft_m[k] = LEfit.LEfit_m[iL][k];
      ft_e[k] = LEfit.LEfit_e[iL][k];
    }

    double me = 0, mt = 0, em = 0, ee = 0, mm = 0;
    for (int j = 0; j <= 1000; ++j) {
      double e = binEdges[i] + j * (binEdges[i+1] - binEdges[i]) / 1000.0;
      me += getAvgPme(e, par[0], par[1], par[2], par[3], ft_m);
      mt += getAvgPmt(e, par[0], par[1], par[2], par[3], ft_m);
      em += getAvgPme(e, par[0], par[1], par[2], par[3], ft_e);
      ee += getAvgPee(e, par[0], par[1], par[2], par[3], ft_e);
      mm += getAvgPmm(e, par[0], par[1], par[2], par[3], ft_m);
    }

    double Pme = me / 1001.0;
    double Pem = em / 1001.0;
    double Pee = ee / 1001.0;
    double Pmm = mm / 1001.0;
    double Pmt = mt / 1001.0;

    CC_tp_me->Add(CC[i].mCC_w, Pme);
    CC_tp_ee->Add(CC[i].eCC, Pee);
    CC_tp_em->Add(CC[i].eCC, Pem);
    CC_tp_mm->Add(CC[i].mCC, Pmm);

    nue_tp_me->Add(nue[i].m_nue_w, Pme);
    nue_tp_mm->Add(nue[i].m_nue, Pmm + Pmt);
    nue_tp_em->Add(nue[i].e_nue_w, Pem);
    nue_tp_ee->Add(nue[i].e_nue, Pee);
  }

  CC_tp_e->Add(CC_tp_me); CC_tp_e->Add(CC_tp_ee);
  CC_tp_m->Add(CC_tp_em); CC_tp_m->Add(CC_tp_mm);
  nue_tp_os->Add(nue_tp_me); nue_tp_os->Add(nue_tp_em);
  nue_tp_unos->Add(nue_tp_ee); nue_tp_unos->Add(nue_tp_mm);
  nue_tp->Add(nue_tp_os); nue_tp->Add(nue_tp_unos);

  // 3. Fill TMatrixD
  TMatrixD prediction(nbins, 1);
  for (int bx = 0; bx < nbins; ++bx) {
    if (bx < nbinsCC)
      prediction[bx][0] = CC_tp_m->GetBinContent(bx + 1);
    else if (bx < 2 * nbinsCC)
      prediction[bx][0] = CC_tp_e->GetBinContent(bx - nbinsCC + 1);
    else
      prediction[bx][0] = nue_tp->GetBinContent(bx - 2 * nbinsCC + 1);
  }

  if ( plotBestfit ) {
    //std::cout << plotName << std::endl;
    TCanvas* c = new TCanvas("c", "", 900, 600);
    plotOscillationComponents( c, par, CC_tp_mm, CC_tp_me, CC_tp_ee, CC_tp_em, nue_tp_mm, nue_tp_me, nue_tp_ee, nue_tp_em, plotBestfit );
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

  //std::cout << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << par[3] << "\n";

  return prediction;
}



double TemplateFitter::getTarget(double par[4], bool originalTgt) {
  // 1) Build targets (full size)
  target   = BuildPredictionVector(par);    // N×1
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
}



double TemplateFitter::CalcChi2Core(const TMatrixD& prediction) const {
  const TMatrixD pred = dropBins_.empty()
      ? prediction
      : ReduceMatrix(prediction, dropBins_);

  TMatrixD diff = pred - myTarget;
  TMatrixD diffT(TMatrixD::kTransposed, diff);
  TMatrixD chi2 = diffT * invmx * diff;
  return chi2(0,0);
}

double TemplateFitter::CalculateChi2(double par[4], bool plotBestfit, const char* outName) {
  //std::cout << outName << std::endl;
  TMatrixD prediction = BuildPredictionVector(par, plotBestfit, outName);
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

bool TemplateFitter::doFit(int nFunctionCalls, int nIterations, double Tolerance, const char* varName[4], double seed[4], double stepSize[4], double (&parBf)[4], double (&parError)[4], double& chi2, bool fixedUt42) {
  ROOT::Math::Minimizer* fitter = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
  fitter->SetMaxFunctionCalls(nFunctionCalls);
  fitter->SetMaxIterations(nIterations);
  fitter->SetTolerance(Tolerance);

  for (int i = 0; i < 4; ++i) {
    fitter->SetVariable(i, varName[i], seed[i], stepSize[i]);
    fitter->SetVariableLimits(i, par_min[i], par_max[i]);
  }
  if (fixedUt42) fitter->FixVariable(2);

  ROOT::Math::Functor functor(this, &TemplateFitter::getChi2, 4);
  fitter->SetFunction(functor);
  fitter->Minimize();

  const double* bestfit = fitter->X();
  const double* bestfitError = fitter->Errors();
  for (int i = 0; i < 4; ++i) {
    parBf[i] = bestfit[i];
    parError[i] = bestfitError[i];
  }

  chi2 = fitter->MinValue();
  return true;
}