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

#include <string>
#include <cmath>
#include <algorithm>

static const char data_path[]  = "/exp/dune/app/users/qvuong/data/lownu";
static const char* kEtheta2_cut = "_2MeV";
// Low-ν (GeV) cuts
static const double kNuCuts[] = {50., 1., 0.5, 0.3, 0.1, 0.08};
// static constexpr int kNNuCut = sizeof(kNuCuts)/sizeof(kNuCuts[0]);

static inline bool close(double a, double b) {
  // abs + relative tolerance; tweak eps if needed
  const double eps_abs = 1e-6;
  const double eps_rel = 1e-12;
  return std::fabs(a - b) <= eps_abs + eps_rel * std::max(std::fabs(a), std::fabs(b));
}

static std::string dm2_tag(double dm2) {
  if (close(dm2, 0.1))  return "smallDm2";
  if (close(dm2, 6.0))  return "medDm2";
  if (close(dm2, 80.0)) return "largeDm2";
  // Fallback if a different value sneaks in:
  char buf[64];
  std::snprintf(buf, sizeof(buf), "dm2_%g", dm2);
  return std::string(buf);
}

// Getting L vs E distribution and its fitting parameters for all 29 energy bands
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

void Scan2D(const std::string& mode, const double* tgtPar, TemplateFitter& tf, int nPoints, int job_id, int total_jobs, bool print, double tgtchi2, const std::string& outname ) {

  // std::ofstream out(Form("%s_%s.txt", outname, mode.c_str()));

  const std::string fname = mode + "_" + outname + ".txt";
  std::ofstream out(fname);

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
  int nucut = 3;
  int nucutNuE = 5;
  bool noNuE = false;
  bool fixedL = false;

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

  std::cout << "noNuE = " << noNuE << "\n";

  // Readin input files
  TFile *CC_f  = new TFile(Form("%s/input_dfiles/CCTemplates_NoBKG%s_FINAL.root", data_path, kEtheta2_cut),"READ");
  TFile *nue_f = new TFile(Form("%s/input_dfiles/nueTemplates_NoBKG%s_FINAL.root", data_path, kEtheta2_cut),"READ");

  // TH1D *hnom = new TH1D("hnom","",nbins,0,nbins);

  // TH1D *CCm_sig   = (TH1D*)CC_f->Get(Form("mElep%d",inu));
  // TH1D *CCe_sig   = (TH1D*)CC_f->Get(Form("eElep%d",inu));
  // TH1D *CCe_bkg_m = (TH1D*)nue_f->Get("mElep");
  // TH1D *CCe_bkg_e = (TH1D*)nue_f->Get("eElep");

  // TH1D *nue_sig = (TH1D*)nue_f->Get("hElep_inverseEtheta2");
  // TH1D *nue_bkg = (TH1D*)CC_f->Get(Form("eElep%d_inverseEtheta2",inu));

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

  std::cout << nue_hm->GetNbinsX() << "\t" << nue_hm->GetNbinsY() << "\n";

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

  std::cout << CCm_tp_cc[0]->GetNbinsX() << "\t" << CCe_tp_inverseEtheta2_cc[0]->GetNbinsX() << "\n";
  std::cout << nue_tp_m[0]->GetNbinsX() << "\t" << nue_tp_m_inverseEtheta2[0]->GetNbinsX() << "\n";


  // TH2D* CCm        = (TH2D*)CC_f->Get(Form("mElepRecoVsEv%d",inu));
  // TH2D* CCm_w      = (TH2D*)CC_f->Get(Form("nc_mElepRecoVsEv%d",inu));
  // TH2D* CCe        = (TH2D*)CC_f->Get(Form("eElepRecoVsEv%d",inu));
  // TH2D* CCe_BKGm   = (TH2D*)nue_f->Get(Form("mElepRecoVsEv"));
  // TH2D* CCe_BKGm_w = (TH2D*)nue_f->Get(Form("mElepRecoVsEv_w"));
  // TH2D* CCe_BKGe   = (TH2D*)nue_f->Get(Form("eElepRecoVsEv"));
  // TH2D* CCe_BKGe_w = (TH2D*)nue_f->Get(Form("eElepRecoVsEv_w"));

  // TH2D* nue_m   = (TH2D*)nue_f->Get(Form("mElepRecoVsEv_inverseEtheta2"));
  // TH2D* nue_m_w = (TH2D*)nue_f->Get(Form("mElepRecoVsEv_inverseEtheta2_w"));
  // TH2D* nue_e   = (TH2D*)nue_f->Get(Form("eElepRecoVsEv_inverseEtheta2"));
  // TH2D* nue_e_w = (TH2D*)nue_f->Get(Form("eElepRecoVsEv_inverseEtheta2_w"));
  // TH2D* nue_BKGe = (TH2D*)CC_f->Get(Form("eElepRecoVsEv%d_inverseEtheta2",inu));

  // TH1D* tp_CCm[nbins_Ev];
  // TH1D* tp_CCm_w[nbins_Ev];
  // TH1D* tp_CCe[nbins_Ev];
  // TH1D* tp_CCe_BKGm[nbins_Ev];
  // TH1D* tp_CCe_BKGm_w[nbins_Ev];
  // TH1D* tp_CCe_BKGe[nbins_Ev];
  // TH1D* tp_CCe_BKGe_w[nbins_Ev];

  // TH1D* tp_nue_m[nbins_Ev];
  // TH1D* tp_nue_m_w[nbins_Ev];
  // TH1D* tp_nue_e[nbins_Ev];
  // TH1D* tp_nue_e_w[nbins_Ev];
  // TH1D* tp_nue_BKGe[nbins_Ev];

  // for(int i=0; i<nbins_Ev; i++) {
  //   const int bx = i + 1;

  //   // CC μ channel
  //   tp_CCm[i]         = (TH1D*)CCm        ->ProjectionY(Form("tp_CCm_inu%d_bin%d",         inu, bx), bx, bx);
  //   tp_CCm_w[i]       = (TH1D*)CCm_w      ->ProjectionY(Form("tp_CCm_w_inu%d_bin%d",       inu, bx), bx, bx);

  //   // CC e channel (signal + ν–e backgrounds)
  //   tp_CCe[i]         = (TH1D*)CCe        ->ProjectionY(Form("tp_CCe_inu%d_bin%d",         inu, bx), bx, bx);
  //   tp_CCe_BKGm[i]    = (TH1D*)CCe_BKGm   ->ProjectionY(Form("tp_CCe_BKGm_inu%d_bin%d",    inu, bx), bx, bx);
  //   tp_CCe_BKGm_w[i]  = (TH1D*)CCe_BKGm_w ->ProjectionY(Form("tp_CCe_BKGm_w_inu%d_bin%d",  inu, bx), bx, bx);
  //   tp_CCe_BKGe[i]    = (TH1D*)CCe_BKGe   ->ProjectionY(Form("tp_CCe_BKGe_inu%d_bin%d",    inu, bx), bx, bx);
  //   tp_CCe_BKGe_w[i]  = (TH1D*)CCe_BKGe_w ->ProjectionY(Form("tp_CCe_BKGe_w_inu%d_bin%d",  inu, bx), bx, bx);

  //   // ν–e sample (signal + CCe→νe background)
  //   tp_nue_m[i]       = (TH1D*)nue_m      ->ProjectionY(Form("tp_nue_m_inu%d_bin%d",       inu, bx), bx, bx);
  //   tp_nue_m_w[i]     = (TH1D*)nue_m_w    ->ProjectionY(Form("tp_nue_m_w_inu%d_bin%d",     inu, bx), bx, bx);
  //   tp_nue_e[i]       = (TH1D*)nue_e      ->ProjectionY(Form("tp_nue_e_inu%d_bin%d",       inu, bx), bx, bx);
  //   tp_nue_e_w[i]     = (TH1D*)nue_e_w    ->ProjectionY(Form("tp_nue_e_w_inu%d_bin%d",     inu, bx), bx, bx);
  //   tp_nue_BKGe[i]    = (TH1D*)nue_BKGe   ->ProjectionY(Form("tp_nue_BKGe_inu%d_bin%d",    inu, bx), bx, bx);

  //   // CC_templates_m[i]    = (TH1D*)CC_hm->ProjectionY(Form("CC_m_bin%d",i+1),i+1,i+1);
  //   // CC_templates_m_nc[i] = (TH1D*)CC_hm_nc->ProjectionY(Form("CC_nc_m_bin%d",i+1),i+1,i+1);
  //   // CC_templates_e[i]    = (TH1D*)CC_he->ProjectionY(Form("CC_e_bin%d",i+1),i+1,i+1);
  //   // nue_templates_m[i]    = (TH1D*)nue_hm->ProjectionY(Form("nue_m_bin%d",i+1),i+1,i+1);
  //   // nue_templates_m_w[i]  = (TH1D*)nue_hm_w->ProjectionY(Form("nue_w_m_bin%d",i+1),i+1,i+1);
  //   // nue_templates_e[i]    = (TH1D*)nue_he->ProjectionY(Form("nue_e_bin%d",i+1),i+1,i+1);
  //   // nue_templates_e_w[i]  = (TH1D*)nue_he_w->ProjectionY(Form("nue_w_e_bin%d",i+1),i+1,i+1);
  // }

  // TemplateGroup_CC ccTemplates[nbins_Ev];
  // TemplateGroup_nue nueTemplates[nbins_Ev];
  // for (int i = 0; i < nbins_Ev; ++i) {
  //   ccTemplates[i] = {tp_CCm[i], tp_CCm_w[i], tp_CCe[i], tp_CCe_BKGm[i], tp_CCe_BKGm_w[i], tp_CCe_BKGe[i], tp_CCe_BKGe_w[i]};
  //   nueTemplates[i] = {tp_nue_m[i], tp_nue_m_w[i], tp_nue_e[i], tp_nue_e_w[i], tp_nue_BKGe[i]};
  // }

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
    // std::cout << energy_bins[b] << "\n";
  }


  TFile *f_fl = new TFile(Form("%s/uncertainties/flux_covmtr/flux_covmtr_wBKG%s_FINAL.root",data_path,kEtheta2_cut),"READ");
  TH2D *fl_cov = (TH2D*)f_fl->Get(Form("hfrcv_nucut%d",nucut));
  TFile *f_xS = new TFile(Form("%s/uncertainties/xS_covmtr/xS_covmtr_wBKG%s_FINAL.root",data_path,kEtheta2_cut), "READ");
  TH2D *xS_cov = (TH2D*)f_xS->Get(Form("FracCov_total_nucut%d",nucut));
  TFile *f_det = new TFile(Form("%s/uncertainties/det_covmtr/det_covmtr_wBKG%s_FINAL.root",data_path,kEtheta2_cut),"READ");
  TH2D *det_cov = (TH2D*)f_det->Get(Form("FracCov_total_nucut%d",nucut));

  fl_cov->SetMaximum(0.02);
  xS_cov->SetMaximum(0.02);
  det_cov->SetMaximum(0.02);

  fl_cov->SetMinimum(0.);
  xS_cov->SetMinimum(0.);
  det_cov->SetMinimum(0.);

  fl_cov->SetTitle(Form("Fractional Flux Covariance Matrix (#nu < %.1f GeV)",kNuCuts[nucut]));
  xS_cov->SetTitle(Form("Fractional Cross-section Covariance Matrix (#nu < %.1f GeV)",kNuCuts[nucut]));
  det_cov->SetTitle(Form("Fractional Detector Covariance Matrix (#nu < %.1f GeV)",kNuCuts[nucut]));

  fl_cov->GetXaxis()->SetTitle("Bin number"); fl_cov->GetYaxis()->SetTitle("Bin number");
  xS_cov->GetXaxis()->SetTitle("Bin number"); xS_cov->GetYaxis()->SetTitle("Bin number");
  det_cov->GetXaxis()->SetTitle("Bin number"); det_cov->GetYaxis()->SetTitle("Bin number");

  gStyle->SetNumberContours(999); TColor::InvertPalette(); 
  // gStyle->SetPalette(kColorPrintableOnGrey); 
  TCanvas *c = new TCanvas("c","",2400,600);
  c->Divide(3,1);
  c->cd(1);
  fl_cov->Draw("colz");
  c->cd(2);
  xS_cov->Draw("colz");
  c->cd(3);
  det_cov->Draw("colz");
  c->SaveAs("matrices.png");

  CovarianceMatrix inputCov;

  for (int i = 0; i < nbins; ++i) {
    // inputCov.FCWeights[i] = (*scales)(u,i);
    inputCov.FCWeights[i] = 0.;
    // if(print) std::cout << i << "\t" << (*scales)(u,i) << "\n";
    // std::cout << fl_cov->GetBinContent(i+1, i+1) << "\t" << det_cov->GetBinContent(i+1, i+1) << "\t" << sig_cov->GetBinContent(i+1, i+1) << "\n";
    for (int j = 0; j < nbins; ++j) {
      inputCov.FluxMx[i][j]         = fl_cov->GetBinContent(i+1, j+1);
      inputCov.CrossSectionMx[i][j] = xS_cov->GetBinContent(i+1, j+1);
      inputCov.DetectorMx[i][j]     = det_cov->GetBinContent(i+1, j+1);
    }
  }

  // std::vector<int> dropBins;
  // // 41..53 (inclusive)
  // for (int i = 41; i <= 53; ++i) dropBins.push_back(i);
  // // last 8 of N=116 -> 108..115
  // // for (int i = 116 - 8; i < 116; ++i) dropBins.push_back(i);
  // // Sort & dedupe, in case ranges overlap
  // std::sort(dropBins.begin(), dropBins.end());
  // dropBins.erase(std::unique(dropBins.begin(), dropBins.end()), dropBins.end());
  // if(print) {  
  //   std::cout << "dropBins (" << dropBins.size() << "): ";
  //   for (size_t i = 0; i < dropBins.size(); ++i) {
  //       if (i) std::cout << ", ";
  //       std::cout << dropBins[i];
  //   }
  //   std::cout << "\n\n";
  // }

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
  double ratio = tf.getTarget( tgtPar, originalTgt );
  double tgtchi2 = tf.CalculateChi2( tgtPar ); 
  if(print) std::cout << "\n thrown/nominal ratio = " << ratio << "\t target chi2 = " << tgtchi2 << std::endl;


  if(hasTestPar){
    plotBestfit = true;
    double par_NO[4] = {0., 0., 0., 0.};
    double testchi2 = tf.CalculateChi2( testPar, plotBestfit, Form("testtarget_nucut%d",nucut) );
    double NOchi2 = tf.CalculateChi2( par_NO, plotBestfit, Form("NOtarget_nucut%d",nucut) ); 
    std::cout << testchi2 << "\t" << NOchi2 << "\n";
    //tf.doFit(100, 10, 0.1, varNames, par_bf.data(), step3.data(), tempPar, tempErr, chi2_penalized, fixedUt42);
  }

  else{
    if (!scanMode.empty()) {
      std::string outname = dm2_tag(tgtPar[3]);         // dm2 is index 3
      // outname += "_nucut" + std::to_string(nucut);
      if (noNuE) outname += "_noNuE";
      Scan2D(scanMode, tgtPar, tf, nPoints, job_id, total_jobs, print, tgtchi2, outname);
    }
  }

  return(0);
}



