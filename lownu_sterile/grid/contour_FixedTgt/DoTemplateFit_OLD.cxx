#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TStyle.h"
#include <TRandom.h>
#include <list>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <limits>
#include <mutex>
#include <omp.h>
#include "ROOT/RConfig.hxx"
#include "ROOT/TThreadExecutor.hxx"
//#include "ROOT/EnableThreadSafety.hxx"

static const char data_path[]  = "/exp/dune/app/users/qvuong/data/lownu";

// Getting L vs E distribution and its fitting parameters for all 29 energy bands
void GetLvsE( double fitPara_m[29][7], double fitPara_e[29][7], TH1D *LEdep_m[29], TH1D *LEdep_e[29] ){
  TFile *ftP_m = new TFile("fitPara_m.root","READ");
  TFile *ftP_e = new TFile("fitPara_e.root","READ");

  TTree *tree_m = (TTree*)ftP_m->Get("pardir");
  TTree *tree_e = (TTree*)ftP_e->Get("pardir");

  double para0_m,para1_m,para2_m,para3_m,para4_m,para5_m,norm_m;
  double para0_e,para1_e,para2_e,para3_e,para4_e,para5_e,norm_e;

  tree_m->SetBranchAddress("para0",&para0_m);
  tree_m->SetBranchAddress("para1",&para1_m);
  tree_m->SetBranchAddress("para2",&para2_m);
  tree_m->SetBranchAddress("para3",&para3_m);
  tree_m->SetBranchAddress("para4",&para4_m);
  tree_m->SetBranchAddress("para5",&para5_m);
  tree_m->SetBranchAddress("norm",&norm_m);

  tree_e->SetBranchAddress("para0",&para0_e);
  tree_e->SetBranchAddress("para1",&para1_e);
  tree_e->SetBranchAddress("para2",&para2_e);
  tree_e->SetBranchAddress("para3",&para3_e);
  tree_e->SetBranchAddress("para4",&para4_e);
  tree_e->SetBranchAddress("para5",&para5_e);
  tree_e->SetBranchAddress("norm",&norm_e);

  for(int i=0;i<tree_m->GetEntries();i++){
    tree_m->GetEntry(i);
    fitPara_m[i][0] = para0_m;
    fitPara_m[i][1] = para1_m;
    fitPara_m[i][2] = para2_m;
    fitPara_m[i][3] = para3_m;
    fitPara_m[i][4] = para4_m;
    fitPara_m[i][5] = para5_m;
    fitPara_m[i][6] = norm_m;
  }
  ftP_m->Close();

  for(int i=0;i<tree_e->GetEntries();i++){
    tree_e->GetEntry(i);
    fitPara_e[i][0] = para0_e;
    fitPara_e[i][1] = para1_e;
    fitPara_e[i][2] = para2_e;
    fitPara_e[i][3] = para3_e;
    fitPara_e[i][4] = para4_e;
    fitPara_e[i][5] = para5_e;
    fitPara_e[i][6] = norm_e;
  }
  ftP_e->Close();


  TFile *f = new TFile("LE_1112_2.root","READ");
  TH2D *LvsE_e = (TH2D*)f->Get("h_e");
  TH2D *LvsE_m = (TH2D*)f->Get("h_m");
  for(int i = 0; i < 29; i++) {
    LEdep_e[i] = (TH1D*)LvsE_e->ProjectionY(Form("LE_e_bin%d",i+1),i+1,i+1);
    LEdep_e[i]->Scale(1./LEdep_e[i]->Integral("width"));

    LEdep_m[i] = (TH1D*)LvsE_m->ProjectionY(Form("LE_m_bin%d",i+1),i+1,i+1);
    LEdep_m[i]->Scale(1./LEdep_m[i]->Integral("width"));
  }
}

std::vector<double> GetLogBinCenters(double min, double max, int nbins) {
  std::vector<double> centers;
  if (nbins < 1 || min <= 0 || max <= 0 || min >= max) return centers;

  double log_min = std::log10(min);
  double log_max = std::log10(max);

  std::vector<double> edges(nbins + 1);
  for (int i = 0; i <= nbins; ++i) {
    edges[i] = log_min + i * (log_max - log_min) / nbins;
  }

  for (int i = 0; i < nbins; ++i) {
    double center_log = (edges[i] + edges[i + 1]) / 2.0;
    centers.push_back(std::pow(10, center_log));
  }

  return centers;
}

void Scan2D(const std::string& mode, const double* tgtPar, TemplateFitter& tf, int nPoints, int job_id, int total_jobs, bool print, double tgtchi2) {
  std::ofstream out(Form("chi2_scan_%s.txt", mode.c_str()));

  auto log_scan = [&](double min, double max) {
    return GetLogBinCenters(min, max, nPoints);
  };

  auto ue42_vals = log_scan(1e-4, 0.7);
  auto um42_vals = log_scan(1e-4, 0.7);
  auto ut42_vals = log_scan(1e-4, 0.7);
  auto dm2_vals  = log_scan(0.1, 100.0);

  long long idx = 0;

  if (mode == "ue42_um42") {
    for (int i = 0; i < nPoints; ++i) {
      double ue42 = ue42_vals[i];
      for (int j = 0; j < nPoints; ++j) {
        double um42 = um42_vals[j];
        if (idx++ % total_jobs != job_id) continue;
        double testPar[4] = {ue42, um42, tgtPar[2], tgtPar[3]};
        double chi2 = tf.CalculateChi2(testPar);
        out << ue42 << "\t" << um42 << "\t" << testPar[2] << "\t" << testPar[3] << "\t" << chi2 << "\t" << tgtchi2 << "\n";
      }
    }
  } else if (mode == "ue42_dm2") {
    for (int i = 0; i < nPoints; ++i) {
      double ue42 = ue42_vals[i];
      for (int j = 0; j < nPoints; ++j) {
        double dm2 = dm2_vals[j];
        if (idx++ % total_jobs != job_id) continue;
        double testPar[4] = {ue42, tgtPar[1], tgtPar[2], dm2};
        double chi2 = tf.CalculateChi2(testPar);
        out << ue42 << "\t" << testPar[1] << "\t" << testPar[2] << "\t" << dm2 << "\t" << chi2 << "\t" << tgtchi2 << "\n";
      }
    }
  } else if (mode == "um42_dm2") {
    for (int i = 0; i < nPoints; ++i) {
      double um42 = um42_vals[i];
      for (int j = 0; j < nPoints; ++j) {
        double dm2 = dm2_vals[j];
        if (idx++ % total_jobs != job_id) continue;
        double testPar[4] = {tgtPar[0], um42, tgtPar[2], dm2};
        double chi2 = tf.CalculateChi2(testPar);
        out << testPar[0] << "\t" << um42 << "\t" << testPar[2] << "\t" << dm2 << "\t" << chi2 << "\t" << tgtchi2 << "\n";
      }
    }
  } else if (mode == "ut42_dm2") {
    for (int i = 0; i < nPoints; ++i) {
      double ut42 = ut42_vals[i];
      for (int j = 0; j < nPoints; ++j) {
        double dm2 = dm2_vals[j];
        if (idx++ % total_jobs != job_id) continue;
        double testPar[4] = {tgtPar[0], tgtPar[1], ut42, dm2};
        double chi2 = tf.CalculateChi2(testPar);
        out << testPar[0] << "\t" << testPar[1] << "\t" << ut42 << "\t" << dm2 << "\t" << chi2 << "\t" << tgtchi2 << "\n";
      }
    }
  }

  out.close();
}


int main(int argc, char const *argv[]) {
  bool print = false;
  bool originalTgt = true;
  bool plotBestfit = false;
  double tgtPar[4] = {0.0, 0.0, 0.0, 0.0};
  double testPar[4] = {0.0, 0.0, 0.0, 0.0};
  int nPoints = 100;
  bool hasTestPar = false;
  int job_id = 0;
  int total_jobs = 1;
  std::string scanMode = "";

  int i = 1;
  while (i < argc) {
    std::string arg = argv[i];

    if (arg == "--tgtPar" && i + 4 < argc) {
      tgtPar[0] = std::atof(argv[i + 1]);
      tgtPar[1] = std::atof(argv[i + 2]);
      tgtPar[2] = std::atof(argv[i + 3]);
      tgtPar[3] = std::atof(argv[i + 4]);
      i += 5;
    } else if (arg == "--testPar" && i + 4 < argc) {
      testPar[0] = std::atof(argv[i + 1]);
      testPar[1] = std::atof(argv[i + 2]);
      testPar[2] = std::atof(argv[i + 3]);
      testPar[3] = std::atof(argv[i + 4]);
      hasTestPar = true;
      i += 5;
    } else if (arg == "--nPoints" && i + 1 < argc) {
      nPoints = std::atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "--print") {
      print = true;
      ++i;
    } else if (arg == "--plotBestfit") {
      plotBestfit = true;
      ++i;
    } else if (arg == "--job" && i + 1 < argc) {
      job_id = std::atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "--njobs" && i + 1 < argc) {
      total_jobs = std::atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "--scanMode" && i + 1 < argc) {
      scanMode = argv[i + 1];
      i += 2;
    } else {
      std::cerr << "Unknown argument: " << arg << std::endl;
      ++i;
    }
  }


  // Readin input files
  TFile *CC_f  = new TFile(Form("CC_output_54bins.root"),"READ");
  TFile *nue_f = new TFile("nue_output_0430.root","READ");

  TH1D *hnom = new TH1D("hnom","",nbins,0,nbins);
  TH1D *nue_nom = (TH1D*)nue_f->Get("hElep");
  TH1D *CCm_nom = (TH1D*)CC_f->Get(Form("wgt_MaCCQE_mElep0_o"));
  TH1D *CCe_nom = (TH1D*)CC_f->Get(Form("wgt_MaCCQE_eElep0_o"));
  for(int i=0; i<nbins; i++){
    if (i < nbinsCC)        hnom->SetBinContent(i+1, CCm_nom->GetBinContent(i+1));
    else if (i < 2*nbinsCC) hnom->SetBinContent(i+1, CCe_nom->GetBinContent(i-nbinsCC+1));
    else                    hnom->SetBinContent(i+1, nue_nom->GetBinContent(i-2*nbinsCC+1));
  }


  TH2D* CC_hm    = (TH2D*)CC_f->Get(Form("mElepRecoVsEv4"));
  TH2D* CC_hm_nc = (TH2D*)CC_f->Get(Form("nc_mElepRecoVsEv4"));
  TH2D* CC_he    = (TH2D*)CC_f->Get(Form("eElepRecoVsEv4"));
  TH2D* nue_hm   = (TH2D*)nue_f->Get(Form("mElepRecoVsEv"));
  TH2D* nue_hm_w = (TH2D*)nue_f->Get(Form("mElepRecoVsEv_w"));
  TH2D* nue_he   = (TH2D*)nue_f->Get(Form("eElepRecoVsEv"));
  TH2D* nue_he_w = (TH2D*)nue_f->Get(Form("eElepRecoVsEv_w"));

  TH1D * CC_templates_m[nbins_Ev];
  TH1D * CC_templates_m_nc[nbins_Ev];
  TH1D * CC_templates_e[nbins_Ev];
  TH1D * nue_templates_m[nbins_Ev];
  TH1D * nue_templates_m_w[nbins_Ev];
  TH1D * nue_templates_e[nbins_Ev];
  TH1D * nue_templates_e_w[nbins_Ev];

  for(int i=0; i<nbins_Ev; i++) {
    CC_templates_m[i]    = (TH1D*)CC_hm->ProjectionY(Form("CC_m_bin%d",i+1),i+1,i+1);
    CC_templates_m_nc[i] = (TH1D*)CC_hm_nc->ProjectionY(Form("CC_nc_m_bin%d",i+1),i+1,i+1);
    CC_templates_e[i]    = (TH1D*)CC_he->ProjectionY(Form("CC_e_bin%d",i+1),i+1,i+1);
    nue_templates_m[i]    = (TH1D*)nue_hm->ProjectionY(Form("nue_m_bin%d",i+1),i+1,i+1);
    nue_templates_m_w[i]  = (TH1D*)nue_hm_w->ProjectionY(Form("nue_w_m_bin%d",i+1),i+1,i+1);
    nue_templates_e[i]    = (TH1D*)nue_he->ProjectionY(Form("nue_e_bin%d",i+1),i+1,i+1);
    nue_templates_e_w[i]  = (TH1D*)nue_he_w->ProjectionY(Form("nue_w_e_bin%d",i+1),i+1,i+1);
  }

  TemplateGroup_CC ccTemplates[nbins_Ev];
  TemplateGroup_nue nueTemplates[nbins_Ev];
  for (int i = 0; i < nbins_Ev; ++i) {
    ccTemplates[i] = {CC_templates_m[i], CC_templates_m_nc[i], CC_templates_e[i]};
    nueTemplates[i] = {nue_templates_m[i], nue_templates_m_w[i], nue_templates_e[i], nue_templates_e_w[i]};
  }

  LE_FitParameters LEPars;
  double LEpar_m[29][7], LEpar_e[29][7];
  TH1D *LEdep_m[29], *LEdep_e[29];
  GetLvsE(LEpar_m, LEpar_e, LEdep_m, LEdep_e);

  LE_FitParameters leParams;
  for (int i = 0; i < 29; ++i) {
    leParams.LEdep_m[i] = LEdep_m[i];
    leParams.LEdep_e[i] = LEdep_e[i];
    for (int j = 0; j < 7; ++j) {
      leParams.LEfit_m[i][j] = LEpar_m[i][j];
      leParams.LEfit_e[i][j] = LEpar_e[i][j];
    }
  }

  double energy_bins[nbins_Ev+1];
  for( int b = 0; b <= nbins_Ev; ++b ) {
    energy_bins[b] = CC_he->GetXaxis()->GetBinLowEdge(b+1);
  }


  TFile *f_fl = new TFile(Form("flux_covmtr_54bins.root"),"READ");
  TH2D *fl_cov = (TH2D*)f_fl->Get("hfrcv4");
  TFile *f_sig = new TFile(Form("xS_unc_54bins.root"), "READ");
  TH2D *sig_cov = (TH2D*)f_sig->Get("frhcv");
  TFile *f_det = new TFile(Form("det_covmtr_54bins.root"),"READ");
  TH2D *det_cov = (TH2D*)f_det->Get("hfrcov4");

  CovarianceMatrix inputCov;

  for (int i = 0; i < nbins; ++i) {
    // inputCov.FCWeights[i] = (*scales)(u,i);
    inputCov.FCWeights[i] = 0.;
    // if(print) std::cout << i << "\t" << (*scales)(u,i) << "\n";
    // std::cout << fl_cov->GetBinContent(i+1, i+1) << "\t" << det_cov->GetBinContent(i+1, i+1) << "\t" << sig_cov->GetBinContent(i+1, i+1) << "\n";
    for (int j = 0; j < nbins; ++j) {
      inputCov.FluxMx[i][j]         = fl_cov->GetBinContent(i+1, j+1);
      inputCov.CrossSectionMx[i][j] = sig_cov->GetBinContent(i+1, j+1);
      inputCov.DetectorMx[i][j]     = det_cov->GetBinContent(i+1, j+1);
    }
  }

  std::vector<int> dropBins;
  // 41..53 (inclusive)
  for (int i = 41; i <= 53; ++i) dropBins.push_back(i);
  // last 8 of N=116 -> 108..115
  // for (int i = 116 - 8; i < 116; ++i) dropBins.push_back(i);
  // Sort & dedupe, in case ranges overlap
  std::sort(dropBins.begin(), dropBins.end());
  dropBins.erase(std::unique(dropBins.begin(), dropBins.end()), dropBins.end());
  if(print) {  
    std::cout << "dropBins (" << dropBins.size() << "): ";
    for (size_t i = 0; i < dropBins.size(); ++i) {
        if (i) std::cout << ", ";
        std::cout << dropBins[i];
    }
    std::cout << "\n\n";
  }


  TemplateFitter tf(ccTemplates, nueTemplates, dropBins);
  tf.setEnergyBins( energy_bins );
  tf.setCovmtr( inputCov );
  tf.setLEParameters( leParams );
  double ratio = tf.getTarget( tgtPar, originalTgt );
  double tgtchi2 = tf.CalculateChi2( tgtPar ); 
  if(print) std::cout << "\n thrown/nominal ratio = " << ratio << "\t target chi2 = " << tgtchi2 << std::endl;


  if(hasTestPar){
    plotBestfit = true;
    double par_NO[4] = {0., 0., 0., 0.};
    double testchi2 = tf.CalculateChi2( testPar, plotBestfit, Form("testtarget_%d") );
    double NOchi2 = tf.CalculateChi2( par_NO, plotBestfit, Form("NOtarget_%d") ); 
    std::cout << testchi2 << "\t" << NOchi2 << "\n";
    //tf.doFit(100, 10, 0.1, varNames, par_bf.data(), step3.data(), tempPar, tempErr, chi2_penalized, fixedUt42);
  }

  else{
    if (!scanMode.empty()) {
      Scan2D(scanMode, tgtPar, tf, nPoints, job_id, total_jobs, print, tgtchi2);
    }
  }

  return(0);
}



