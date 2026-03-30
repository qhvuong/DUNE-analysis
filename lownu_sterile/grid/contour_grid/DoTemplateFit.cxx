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



// Getting L vs E distribution and its fitting parameters for all nbins_Ev_LE energy bands
void GetLvsE( double fitPara_m[nbins_Ev_LE][7], double fitPara_e[nbins_Ev_LE][7], TH1D *LEdep_m[nbins_Ev_LE], TH1D *LEdep_e[nbins_Ev_LE] ){
  TFile *ftP_m = new TFile(Form("fitPara_m_%dbins.root", nbins_Ev_LE),"READ");
  TFile *ftP_e = new TFile(Form("fitPara_e_%dbins.root", nbins_Ev_LE),"READ");

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


  // TFile *f = new TFile("LE_1112_2.root","READ");
  TFile *f = new TFile(Form("LE_%dbins.root", nbins_Ev_LE),"READ");
  TH2D *LvsE_e = (TH2D*)f->Get("h_e");
  TH2D *LvsE_m = (TH2D*)f->Get("h_m");
  for(int i = 0; i < nbins_Ev_LE; i++) {
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
  auto dm2_vals  = log_scan(0.01, 100.0);

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


void Scan4DvsNull(TemplateFitter& tf, int nPoints, int job_id, int total_jobs, bool print) {
  std::ofstream out(Form("chi2_scan4D_vs_null_noUt42_%d.txt", job_id));

  auto log_scan = [&](double min, double max) {
    return GetLogBinCenters(min, max, nPoints);
  };

  auto ue42_vals = log_scan(1e-4, 0.7);
  auto um42_vals = log_scan(1e-4, 0.7);
  auto dm2_vals  = log_scan(0.01, 100.0);
  std::vector<double> ut42_vals = {0.0};  // only 2 values

  long long idx = 0;
  double par_NO[4] = {0., 0., 0., 0.};

  for (double ut42 : ut42_vals) {
    for (double ue42 : ue42_vals) {
      for (double um42 : um42_vals) {
        for (double dm2 : dm2_vals) {
          if (idx++ % total_jobs != job_id) continue;

          double tgtPar[4] = {ue42, um42, ut42, dm2};
          TemplateFitter tf_copy = tf;
          tf_copy.getTarget(tgtPar, /*originalTgt=*/false); 
          double thisChi2 = tf_copy.CalculateChi2(tgtPar);
          double nullChi2 = tf_copy.CalculateChi2(par_NO);
          // nullChi2 = (nullChi2 < 1e-10) ? 0.0 : nullChi2;
          // double dChi2 = thisChi2 - nullChi2;

          out << ue42 << "\t" << um42 << "\t" << ut42 << "\t" << dm2
              << "\t" << thisChi2 << "\t" << nullChi2 << "\n";
        }
      }
    }
  }

  out.close();
}



int main(int argc, char const *argv[]) {
  bool print = false;
  int nPoints = 100;
  int point_idx = -1;
  double ut42 = 0.0;  // fixed
  bool noNuE = false;
  int nucut = 3;
  int nucutNuE = 5;
  double testPar[4] = {0.0, 0.0, 0.0, 0.0};
  bool hasTestPar = false;
  bool fixedL = false;

  int i = 1;
  while (i < argc) {
    std::string arg = argv[i];

    if (arg == "--point" && i + 1 < argc) {
      point_idx = std::atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "--testPar" && i + 4 < argc) {
      testPar[0] = std::atof(argv[i + 1]);
      testPar[1] = std::atof(argv[i + 2]);
      testPar[2] = std::atof(argv[i + 3]);
      testPar[3] = std::atof(argv[i + 4]);
      hasTestPar = true;
      i += 5;
    } else if (arg == "--ut42" && i + 1 < argc) {
      ut42 = std::atof(argv[i + 1]);
      i += 2;
    } else if (arg == "--nPoints" && i + 1 < argc) {
      nPoints = std::atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "--nucut" && i + 1 < argc) {
      nucut = std::atoi(argv[i + 1]);
      i += 2;
    } else if (arg == "--print") {
      print = true;
      ++i;
    } else if (arg == "--noNuE") {
      noNuE = true;
      i += 1;
    } else if (arg == "--fixedL") {              
      fixedL = true;
      i += 1;
    } else {
      std::cerr << "Unknown argument: " << arg << std::endl;
      ++i;
    }
  }

  int total_points = nPoints * nPoints * nPoints;
  if (point_idx >= total_points) {
    std::cerr << "[ERROR] point_idx can't be >= " << total_points << std::endl;
    return 1;
  }




  // Readin input files
  TFile *CC_f = new TFile("CCTemplates_NoBKG_2MeV_FINAL.root","READ");
  TFile *nue_f = new TFile("nueTemplates_NoBKG_2MeV_FINAL.root","READ");

  // TH1D *hnom = new TH1D("hnom","",nbins,0,nbins);
  // TH1D *nue_nom = (TH1D*)nue_f->Get("hElep");
  // TH1D *CCm_nom = (TH1D*)CC_f->Get(Form("wgt_MaCCQE_mElep0_o"));
  // TH1D *CCe_nom = (TH1D*)CC_f->Get(Form("wgt_MaCCQE_eElep0_o"));
  // for(int i=0; i<nbins; i++){
  //   if (i < nbinsCC)        hnom->SetBinContent(i+1, CCm_nom->GetBinContent(i+1));
  //   else if (i < 2*nbinsCC) hnom->SetBinContent(i+1, CCe_nom->GetBinContent(i-nbinsCC+1));
  //   else                    hnom->SetBinContent(i+1, nue_nom->GetBinContent(i-2*nbinsCC+1));
  // }


  // CCvu templates
  TH2D* CC_hm_cc = (TH2D*)CC_f->Get(Form("mElepRecoVsEv_CC_nucut%d",nucut));
  TH2D* CC_hm_nc = (TH2D*)CC_f->Get(Form("mElepRecoVsEv_NC_nucut%d",nucut));
  TH2D* CC_hm_w  = (TH2D*)CC_f->Get(Form("mElepRecoVsEv_w_nucut%d",nucut));

  // CCve templates
  TH2D* CC_he_cc                = (TH2D*)CC_f->Get(Form("eElepRecoVsEv_CC_nucut%d",nucut));
  TH2D* CC_he_nc                = (TH2D*)CC_f->Get(Form("eElepRecoVsEv_NC_nucut%d",nucut));
  TH2D* CC_he_inverseEtheta2_cc = (TH2D*)CC_f->Get(Form("eElepRecoVsEv_inverseEtheta2_CC_nucut%d",nucutNuE));
  TH2D* CC_he_inverseEtheta2_nc = (TH2D*)CC_f->Get(Form("eElepRecoVsEv_inverseEtheta2_NC_nucut%d",nucutNuE));

  // v-e templates
  TH2D* nue_hm                  = (TH2D*)nue_f->Get(Form("mElepRecoVsEv"));
  TH2D* nue_hm_w                = (TH2D*)nue_f->Get(Form("mElepRecoVsEv_w"));
  TH2D* nue_he                  = (TH2D*)nue_f->Get(Form("eElepRecoVsEv"));
  TH2D* nue_he_w                = (TH2D*)nue_f->Get(Form("eElepRecoVsEv_w"));
  TH2D* nue_hm_inverseEtheta2   = (TH2D*)nue_f->Get(Form("mElepRecoVsEv_inverseEtheta2"));
  TH2D* nue_hm_inverseEtheta2_w = (TH2D*)nue_f->Get(Form("mElepRecoVsEv_inverseEtheta2_w"));
  TH2D* nue_he_inverseEtheta2   = (TH2D*)nue_f->Get(Form("eElepRecoVsEv_inverseEtheta2"));
  TH2D* nue_he_inverseEtheta2_w = (TH2D*)nue_f->Get(Form("eElepRecoVsEv_inverseEtheta2_w"));

  // std::cout << nue_hm->GetNbinsX() << "\t" << nue_hm->GetNbinsY() << "\n";

  TH1D * CCm_tp_cc[nbins_Ev];
  TH1D * CCm_tp_nc[nbins_Ev];
  TH1D * CCm_tp_w[nbins_Ev];

  TH1D * CCe_tp_cc[nbins_Ev];
  TH1D * CCe_tp_nc[nbins_Ev];

  TH1D * CCe_tp_inverseEtheta2_cc[nbins_Ev];
  TH1D * CCe_tp_inverseEtheta2_nc[nbins_Ev];

  TH1D * nue_tp_m[nbins_Ev];
  TH1D * nue_tp_m_w[nbins_Ev];
  TH1D * nue_tp_e[nbins_Ev];
  TH1D * nue_tp_e_w[nbins_Ev];
  TH1D * nue_tp_m_inverseEtheta2[nbins_Ev];
  TH1D * nue_tp_m_inverseEtheta2_w[nbins_Ev];
  TH1D * nue_tp_e_inverseEtheta2[nbins_Ev];
  TH1D * nue_tp_e_inverseEtheta2_w[nbins_Ev];

  for(int i=0; i<nbins_Ev; i++) {
    CCm_tp_cc[i]                = (TH1D*)CC_hm_cc               ->ProjectionY(Form("CC_mbin_cc%d",i+1),i+1,i+1);
    CCm_tp_nc[i]                = (TH1D*)CC_hm_nc               ->ProjectionY(Form("CC_mbin_nc%d",i+1),i+1,i+1);
    CCm_tp_w[i]                 = (TH1D*)CC_hm_w                ->ProjectionY(Form("CC_mbin_w%d",i+1),i+1,i+1);
    CCe_tp_cc[i]                = (TH1D*)CC_he_cc               ->ProjectionY(Form("CC_ebin_cc%d",i+1),i+1,i+1);
    CCe_tp_nc[i]                = (TH1D*)CC_he_nc               ->ProjectionY(Form("CC_ebin_nc%d",i+1),i+1,i+1);

    CCe_tp_inverseEtheta2_cc[i] = (TH1D*)CC_he_inverseEtheta2_cc->ProjectionY(Form("CC_ebin_inv_cc%d",i+1),i+1,i+1);
    CCe_tp_inverseEtheta2_nc[i] = (TH1D*)CC_he_inverseEtheta2_nc->ProjectionY(Form("CC_ebin_inv_nc%d",i+1),i+1,i+1);

    nue_tp_m[i]   = (TH1D*)nue_hm  ->ProjectionY(Form("nue_mbin%d",i+1),i+1,i+1);
    nue_tp_m_w[i] = (TH1D*)nue_hm_w->ProjectionY(Form("nue_mbin_w%d",i+1),i+1,i+1);
    nue_tp_e[i]   = (TH1D*)nue_he  ->ProjectionY(Form("nue_ebin%d",i+1),i+1,i+1);
    nue_tp_e_w[i] = (TH1D*)nue_he_w->ProjectionY(Form("nue_ebin_w%d",i+1),i+1,i+1);

    nue_tp_m_inverseEtheta2[i]   = (TH1D*)nue_hm_inverseEtheta2  ->ProjectionY(Form("nue_mbin_inv%d", i+1), i+1, i+1);
    nue_tp_m_inverseEtheta2_w[i] = (TH1D*)nue_hm_inverseEtheta2_w->ProjectionY(Form("nue_mbin_inv_w%d", i+1), i+1, i+1);
    nue_tp_e_inverseEtheta2[i]   = (TH1D*)nue_he_inverseEtheta2  ->ProjectionY(Form("nue_ebin_inv%d", i+1), i+1, i+1);
    nue_tp_e_inverseEtheta2_w[i] = (TH1D*)nue_he_inverseEtheta2_w->ProjectionY(Form("nue_ebin_inv_w%d", i+1), i+1, i+1);
  }


  TemplateGroup_CC ccTemplates[nbins_Ev];
  TemplateGroup_nue nueTemplates[nbins_Ev];
  for (int i = 0; i < nbins_Ev; ++i) {
    ccTemplates[i] = {CCm_tp_cc[i], CCm_tp_nc[i], CCm_tp_w[i], CCe_tp_cc[i], CCe_tp_nc[i], nue_tp_m[i], nue_tp_m_w[i], nue_tp_e[i], nue_tp_e_w[i]};
    nueTemplates[i] = {nue_tp_m_inverseEtheta2[i], nue_tp_m_inverseEtheta2_w[i], nue_tp_e_inverseEtheta2[i], nue_tp_e_inverseEtheta2_w[i], CCe_tp_inverseEtheta2_cc[i], CCe_tp_inverseEtheta2_nc[i]};
  }

  LE_FitParameters LEPars;
  double LEpar_m[nbins_Ev_LE][7], LEpar_e[nbins_Ev_LE][7];
  TH1D *LEdep_m[nbins_Ev_LE], *LEdep_e[nbins_Ev_LE];
  GetLvsE(LEpar_m, LEpar_e, LEdep_m, LEdep_e);

  LE_FitParameters leParams;
  for (int i = 0; i < nbins_Ev_LE; ++i) {
    leParams.LEdep_m[i] = LEdep_m[i];
    leParams.LEdep_e[i] = LEdep_e[i];
    for (int j = 0; j < 7; ++j) {
      leParams.LEfit_m[i][j] = LEpar_m[i][j];
      leParams.LEfit_e[i][j] = LEpar_e[i][j];
    }
  }

  double energy_bins[nbins_Ev+1];
  for( int b = 0; b <= nbins_Ev; ++b ) {
    energy_bins[b] = CC_he_cc->GetXaxis()->GetBinLowEdge(b+1);
  }


  // Get systematic uncertainties
  TFile *f_xS  = new TFile(Form("xS_covmtr_wBKG_2MeV_FINAL.root"), "READ");
  TFile *f_fl  = new TFile(Form("flux_covmtr_wBKG_2MeV_FINAL.root"),"READ");
  TFile *f_det = new TFile(Form("det_covmtr_wBKG_2MeV_FINAL.root"),"READ");

  TH2D *hxS  = (TH2D*)f_xS->Get(Form("FracCov_total_nucut%d",nucut));
  TH2D *hfl  = (TH2D*)f_fl->Get(Form("hfrcv_nucut%d",nucut));
  TH2D *hdet = (TH2D*)f_det->Get(Form("FracCov_total_nucut%d",nucut));

  // TFile *scales_f = new TFile(Form("FC%s_wBKG_2MeV.root",FCsample.c_str()),"READ");
  // TMatrixD *scales = (TMatrixD*)scales_f->Get("hscales");

  CovarianceMatrix inputCov;

  for (int i = 0; i < nbins; ++i) {
    inputCov.FCWeights[i] = 0.;
    for (int j = 0; j < nbins; ++j) {
      inputCov.FluxMx[i][j]         = hfl->GetBinContent(i+1, j+1);
      inputCov.CrossSectionMx[i][j] = hxS->GetBinContent(i+1, j+1);
      inputCov.DetectorMx[i][j]     = hdet->GetBinContent(i+1, j+1);
    }
  }

  std::vector<int> dropBins;
  // 41..53 (inclusive)
  for (int i = 41; i <= 53; ++i) dropBins.push_back(i);
  // last 7 of N=115 -> 108..114
  if(noNuE) {
    for (int i = nbins - nbins_nue; i < nbins; ++i) dropBins.push_back(i);
  } 
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

  TemplateFitter tf(ccTemplates, nueTemplates, dropBins, fixedL);
  tf.setEnergyBins( energy_bins );
  tf.setCovmtr( inputCov );
  tf.setLEParameters( leParams );

  if(print) std::cout << nPoints << "\n";

  if(point_idx >= 0) {
    std::vector<double> ue42_vals = GetLogBinCenters(1e-4, 0.7, nPoints);
    std::vector<double> um42_vals = GetLogBinCenters(1e-4, 0.7, nPoints);
    std::vector<double> dm2_vals  = GetLogBinCenters(0.1, 100.0, nPoints);
    

    int i_ue = point_idx / (nPoints * nPoints);
    int j_um = (point_idx / nPoints) % nPoints;
    int k_dm = point_idx % nPoints;

    double ue42 = ue42_vals[i_ue];
    double um42 = um42_vals[j_um];
    double dm2  = dm2_vals[k_dm];


    double tgtPar[4] = {ue42, um42, ut42, dm2};
    if(hasTestPar){
      for(int i=0; i<4; i++) tgtPar[i] = testPar[i];
    }
    // double tgtPar[4] = {0.04, 0.01, 0.2, 6.0};
    TemplateFitter tf_copy = tf;  // isolate per-target
    tf_copy.getTarget(tgtPar, false);

    double nullPar[4] = {0., 0., 0., 0.};
    double nullChi2  = tf_copy.CalculateChi2(nullPar);
    double thisChi2  = tf_copy.CalculateChi2(tgtPar);

    std::cout << point_idx << "\t" << ue42 << "\t" << um42 << "\t" << ut42 << "\t" << dm2 << "\t" << thisChi2 << "\t" << nullChi2 << "\n";

    std::ofstream out(Form("output.txt"));
    out << ue42 << "\t" << um42 << "\t" << ut42 << "\t" << dm2 << "\t" << thisChi2 << "\t" << nullChi2 << "\n";
    out.close();
  }

  return(0);
}



