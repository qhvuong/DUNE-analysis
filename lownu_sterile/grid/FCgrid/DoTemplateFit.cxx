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

using namespace FitUtils;

static const char data_path[] = "/exp/dune/app/users/qvuong/data/lownu";
const char* varNames[4] = {"ue4", "um4", "ut4", "dm2"};

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

double RunThreeStageFit(TemplateFitter& tf, std::vector<double>& par_bf, std::vector<double>& par_err,
                        int nPoints_dm2, bool fixedUt42, bool noLimit, bool print) {
  double bestChi2 = std::numeric_limits<double>::max();
  std::vector<double> tempPar(4), tempErr(4);
  double chi2_penalized, chi2;

  if (print) std::cout << "\n\n STAGE 1 FIT... \n";
  auto seedSets = GetStage1SeedSets(noLimit);
  auto step1 = GetStage1StepSizes();

  #pragma omp parallel
  {
    if (print) {
      #pragma omp single
      std::cout << "OpenMP using " << omp_get_num_threads() << " threads\n";
    }

    std::vector<double> tempPar_local(4), tempErr_local(4);
    double bestChi2_local = std::numeric_limits<double>::max();
    std::vector<double> bestPar_local(4), bestErr_local(4);

    #pragma omp for
    for (int i = 0; i < static_cast<int>(seedSets.size()); ++i) {
      std::vector<double> seed = seedSets[i];

      if (print) {
        #pragma omp critical
        {
          std::cout << "[Stage 1: Seed " << i << "] Seeding: ";
          for (double s : seed) std::cout << s << " ";
          std::cout << std::endl;
        }
      }

      TemplateFitter tf_local = tf;
      double chi2_penalized_local, chi2_local;

      tf_local.doFit(Stage1MaxCalls, Stage1MaxIters, Stage1Tolerance,
                     varNames, seed.data(), step1.data(),
                     tempPar_local, tempErr_local,
                     chi2_penalized_local, fixedUt42);

      chi2_local = tf_local.CalculateChi2(tempPar_local.data());

      if (print) {
        #pragma omp critical
        {
          std::cout << "  Result: ";
          for (double val : tempPar_local) std::cout << val << " ";
          std::cout << " | chi2 = " << chi2_local
                    << " (penalized chi2 = " << chi2_penalized_local << ")\n";
        }
      }

      if (chi2_local < bestChi2_local) {
        bestChi2_local = chi2_local;
        bestPar_local = tempPar_local;
        bestErr_local = tempErr_local;
      }
    }

    #pragma omp critical
    {
      if (bestChi2_local < bestChi2) {
        bestChi2 = bestChi2_local;
        par_bf = bestPar_local;
        par_err = bestErr_local;
      }
    }
  }

  if (print) {
    std::cout << "\n>>> After Stage 1 Best Fit:\n";
    for (int i = 0; i < 4; ++i)
      std::cout << "  par[" << i << "] = " << par_bf[i] << " ± " << par_err[i] << "\n";
    std::cout << "  Best chi2 = " << bestChi2 << "\n\n\n";
  }

  if (print) std::cout << "\n\n STAGE 2 FIT... \n";

  auto dm2Scan = GetDM2SeedsAndSteps(nPoints_dm2);
  std::vector<double> bestPar_local(4), bestErr_local(4);
  double bestChi2_local = std::numeric_limits<double>::max();

  #pragma omp parallel
  {
    TemplateFitter tf_local = tf;
    std::vector<double> thread_bestPar(4), thread_bestErr(4);
    double thread_bestChi2 = std::numeric_limits<double>::max();

    #pragma omp for
    for (int idx = 0; idx < static_cast<int>(dm2Scan.size()); ++idx) {
      auto [dm2, step] = dm2Scan[idx];

      std::vector<double> seed = {par_bf[0], par_bf[1], par_bf[2], dm2};
      double stepSize[4] = {0.01, 0.01, 0.01, step};
      std::vector<double> localPar(4), localErr(4);
      double chi2_penalized_local;

      if (print) {
        #pragma omp critical
        {
          std::cout << "[Stage 2: dm2 = " << dm2 << "] Seeding: ";
          for (double s : seed) std::cout << s << " ";
          std::cout << std::endl;
        }
      }

      tf_local.doFit(Stage2MaxCalls, Stage2MaxIters, Stage2Tolerance,
                     varNames, seed.data(), stepSize,
                     localPar, localErr, chi2_penalized_local, fixedUt42);

      double chi2_local = tf_local.CalculateChi2(localPar.data());

      if (print) {
        #pragma omp critical
        {
          std::cout << "  Result: ";
          for (double val : localPar) std::cout << val << " ";
          std::cout << " | chi2 = " << chi2_local
                    << " (penalized chi2 = " << chi2_penalized_local << ")\n";
        }
      }

      if (chi2_local < thread_bestChi2) {
        thread_bestChi2 = chi2_local;
        thread_bestPar = localPar;
        thread_bestErr = localErr;
      }
    }

    #pragma omp critical
    {
      if (thread_bestChi2 < bestChi2_local) {
        bestChi2_local = thread_bestChi2;
        bestPar_local = thread_bestPar;
        bestErr_local = thread_bestErr;
      }
    }
  }

  bestChi2 = bestChi2_local;
  par_bf = bestPar_local;
  par_err = bestErr_local;

  if (print) {
    std::cout << "\n>>> After Stage 2 Best Fit:\n";
    for (int i = 0; i < 4; ++i)
      std::cout << "  par[" << i << "] = " << par_bf[i] << " ± " << par_err[i] << "\n";
    std::cout << "  Best chi2 = " << bestChi2 << "\n\n\n";
  }

  if (print) std::cout << "\n\n STAGE 3 FIT... \n";

  auto step3 = GetStage3StepSizes();
  tf.doFit(Stage3MaxCalls, Stage3MaxIters, Stage3Tolerance,
           varNames, par_bf.data(), step3.data(),
           tempPar, tempErr, chi2_penalized, fixedUt42);

  chi2 = tf.CalculateChi2(tempPar.data());

  if (chi2 < bestChi2) {
    bestChi2 = chi2;
    par_bf = tempPar;
    par_err = tempErr;
  }

  if (print) {
    std::cout << "\n>>> After Stage 3 Best Fit:\n";
    for (int i = 0; i < 4; ++i)
      std::cout << "  par[" << i << "] = " << par_bf[i] << " ± " << par_err[i] << "\n";
    std::cout << "  Best chi2 = " << bestChi2 << "\n\n\n";
  }

  return bestChi2;
}



int main(int argc, char const *argv[])
{
  ROOT::EnableThreadSafety();


  int u = -1;                               // Default value
  bool fixedUt42 = false;                   // Default: Ut42 is free
  bool print = false;                       // Default is false, not printing anything. If true, print out all fitting steps for debugging
  bool originalTgt = false;                 // If called, perform the fit with the original target, else do it with the FC thrown one
  bool plotBestfit = false;
  //bool plotTarget = false;
  double tgtPar[4] = {0.0, 0.0, 0.0, 0.0};  // Default array
  double testPar[4] = {0.0, 0.0, 0.0, 0.0}; // Default array
  int nPoints_dm2 = 8;                      // Default value is 8, possible values: 4, 8, 16
  bool hasTestPar = false;
  std::string FCsample = "tot";             // Default is throwing all systematics: stat+flux+det+sig
                                            // Other allowed names: "stat, flux, det, xS"
  // std::string testName = "AsimovS";         // Default is testing at Asimov sensitivity
  //                                           // Other allowed names: "AsimovM, AsimovL, phyS, phyM, phyL"
  bool grid = false;
  bool trueSeeding = false;
  bool noNuE = false;
  int nucut = 3;
  int nucutNuE = 5;
  bool noUpperLimits = false;

  int i = 1;
  while (i < argc) {
    std::string arg = argv[i];

    if (arg == "--u" && i + 1 < argc) {
      u = std::atoi(argv[i + 1]);
      i += 2;
    }
    else if (arg == "--tgtPar" && i + 4 < argc) {
      tgtPar[0] = std::atof(argv[i + 1]);
      tgtPar[1] = std::atof(argv[i + 2]);
      tgtPar[2] = std::atof(argv[i + 3]);
      tgtPar[3] = std::atof(argv[i + 4]);
      i += 5;
    }
    else if (arg == "--testPar" && i + 4 < argc) {
      testPar[0] = std::atof(argv[i + 1]);
      testPar[1] = std::atof(argv[i + 2]);
      testPar[2] = std::atof(argv[i + 3]);
      testPar[3] = std::atof(argv[i + 4]);
      hasTestPar = true;
      i += 5;
    }
    else if (arg == "--nPoints_dm2" && i + 1 < argc) {
      nPoints_dm2 = std::atoi(argv[i + 1]);
      i += 2;
    }
    else if (arg == "--nucut" && i + 1 < argc) {
      nucut = std::atoi(argv[i + 1]);
      i += 2;
    }
    else if (arg == "--nucutNuE" && i + 1 < argc) {
      nucutNuE = std::atoi(argv[i + 1]);
      i += 2;
    }
    else if (arg == "--fixedUt42") {
      fixedUt42 = true;
      i += 1;
    }
    else if (arg == "--trueSeeding") {
      trueSeeding = true;
      i += 1;
    }
    else if (arg == "--print") {
      print = true;
      i += 1;
    }
    else if (arg == "--grid") {
      grid = true;
      i += 1;
    }
    else if (arg == "--originalTgt") {              
      originalTgt = true;
      i += 1;
    }
    else if (arg == "--noUpperLimits") {              
      noUpperLimits = true;
      i += 1;
    }
    else if (arg == "--plotBestfit") {              
      plotBestfit = true;
      i += 1;
    }
    // else if (arg == "--plotTarget") {              
    //   plotTarget = true;
    //   i += 1;
    // }
    else if (arg == "--FCsample" && i + 1 < argc) {
      FCsample = argv[i + 1];
      i += 2;
    }
    else if (arg == "--noNuE") {
      noNuE = true;
      i += 1;
    }
    // else if (arg == "--testName" && i + 1 < argc) {
    //   testName = argv[i + 1];
    //   i += 2;
    // }
    else {
      std::cerr << "Unknown argument: " << arg << std::endl;
      ++i;
    }

  }

  if(print) std::cout << "\n nPoints_dm2 = " << nPoints_dm2 << "\n\n";

  // Readin input files
  // TFile *CC_f  = new TFile(Form("CC_output_54bins.root"),"READ");
  // TFile *nue_f = new TFile("nue_output_0430.root","READ");

  TFile *CC_f = new TFile("CCTemplates_NoBKG_2MeV_FINAL.root","READ");
  TFile *nue_f = new TFile("nueTemplates_NoBKG_2MeV_FINAL.root","READ");


  // TH1D *hnom = new TH1D("hnom","",nbins,0,nbins);
  // TH1D *nue_nom = (TH1D*)nue_f->Get("hElep");

  // TH1D *CCm_nom = (TH1D*)CC_f->Get(Form("mElep_nucut%d",nucut));
  // TH1D *CCe_nom = (TH1D*)CC_f->Get(Form("eElep_nucut%d",nucut));

  // for(int i=0; i<nbins; i++){
  //   if (i < nbinsCC)        hnom->SetBinContent(i+1, CCm_nom->GetBinContent(i+1));
  //   else if (i < 2*nbinsCC) hnom->SetBinContent(i+1, CCe_nom->GetBinContent(i-nbinsCC+1));
  //   else                    hnom->SetBinContent(i+1, nue_nom->GetBinContent(i-2*nbinsCC+1));
  // }


  // CCvu templates
  TH2D* CC_hm_cc = (TH2D*)CC_f->Get(Form("mElepRecoVsEv_CC_nucut%d",nucut));
  TH2D* CC_hm_nc = (TH2D*)CC_f->Get(Form("mElepRecoVsEv_NC_nucut%d",nucut));
  TH2D* CC_hm_w  = (TH2D*)CC_f->Get(Form("mElepRecoVsEv_w_nucut%d",nucut));

  // TH2D* histo = (TH2D*)CC_f->Get(Form("mElepRecoVsEv_NC_nucut4"));
  // TCanvas *c = new TCanvas("c","",800,600);
  // TH1D *hx = (TH1D*)histo->ProjectionY("hx",  /*firstY*/1, /*lastY*/histo->GetNbinsX());
  // hx->SetDirectory(nullptr); // detach from current file/dir
  // hx->Draw();
  // c->SaveAs("hx.png");

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

  TFile *scales_f = new TFile(Form("FC%s_wBKG_2MeV_FINAL.root",FCsample.c_str()),"READ");
  TMatrixD *scales = (TMatrixD*)scales_f->Get("hscales");

  CovarianceMatrix inputCov;

  for (int i = 0; i < nbins; ++i) {
    inputCov.FCWeights[i] = (*scales)(u,i);
    // if(print) std::cout << i << "\t" << (*scales)(u,i) << "\n";
    // std::cout << fl_cov->GetBinContent(i+1, i+1) << "\t" << det_cov->GetBinContent(i+1, i+1) << "\t" << sig_cov->GetBinContent(i+1, i+1) << "\n";
    for (int j = 0; j < nbins; ++j) {
      inputCov.FluxMx[i][j]         = hfl->GetBinContent(i+1, j+1);
      inputCov.CrossSectionMx[i][j] = hxS->GetBinContent(i+1, j+1);
      inputCov.DetectorMx[i][j]     = hdet->GetBinContent(i+1, j+1);
    }

    // std::cout << inputCov.FluxMx[i][i] << "\t" << inputCov.DetectorMx[i][i] << "\t" << inputCov.CrossSectionMx[i][i] << "\n";
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

  TemplateFitter tf(ccTemplates, nueTemplates, dropBins, noUpperLimits);
  tf.setEnergyBins( energy_bins );
  tf.setCovmtr( inputCov );
  tf.setLEParameters( leParams );
  double ratio = tf.getTarget( tgtPar, originalTgt );
  double tgtchi2 = tf.CalculateChi2( tgtPar ); 
  if(print) {
    std::cout << "\n thrown/nominal ratio = " << ratio << "\t target chi2 = " << tgtchi2 << std::endl;
    std::cout << "target Para: " << tgtPar[0] << "\t" << tgtPar[1] << "\t" << tgtPar[2] << "\t" << tgtPar[3] << "\n";
  }

  if(hasTestPar){
    plotBestfit = true;
    double par_NO[4] = {0., 0., 0., 0.};
    double testchi2 = tf.CalculateChi2( testPar, plotBestfit, Form("testtarget_%d",u) );
    double NOchi2 = tf.CalculateChi2( par_NO, plotBestfit, Form("NOtarget_%d",u) ); 
    std::cout << "testchi2 = " << testchi2 << "\t\t NOchi2 = " << NOchi2 << "\n";
  }

  else {
    std::vector<double> bestfit(4), bestfitError(4);
    double bfchi2, chi2_penalized;

    if (trueSeeding) {
      auto step3 = GetStage3StepSizes();
      tf.doFit(Stage3MaxCalls, Stage3MaxIters, Stage3Tolerance,
              varNames, tgtPar, step3.data(),
              bestfit, bestfitError, chi2_penalized, fixedUt42);
      bfchi2 = tf.CalculateChi2(bestfit.data());
    }

    else {  
      bfchi2 = RunThreeStageFit(tf, bestfit, bestfitError, nPoints_dm2, fixedUt42, noUpperLimits, print);

      tf.CalculateChi2(bestfit.data(), plotBestfit,
                      Form("bestfit_%d", u));
    }

    if (!grid) {
      std::ofstream myfile(Form("job_output_FAILED_%d.txt", u));
      myfile << u << "\t"
            << bestfit[0] << "\t" << bestfit[1] << "\t"
            << bestfit[2] << "\t" << bestfit[3] << "\t"
            << bfchi2 << "\t" << tgtchi2 << "\t" << ratio << std::endl;
      myfile.close(); 
    }
    else {
      std::cout << u << "\t"
                << bestfit[0] << "\t" << bestfit[1] << "\t"
                << bestfit[2] << "\t" << bestfit[3] << "\t"
                << bfchi2 << "\t" << tgtchi2 << "\t" << ratio << std::endl;
    }
  }

  CC_f->Close();  delete CC_f;
  nue_f->Close(); delete nue_f;
  f_xS->Close();  delete f_xS;
  f_fl->Close();  delete f_fl;
  f_det->Close(); delete f_det;

  return(0);
}



