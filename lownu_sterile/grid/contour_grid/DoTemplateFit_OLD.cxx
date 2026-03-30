#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TStyle.h"
#include <TRandom.h>
#include <list>
#include <TMath.h>
#include <iostream>
#include <fstream>
#include <DUNEStyle.h>

using namespace std;

int main(int argc, char const *argv[])
{

  float Ue42, Um42, Ut42, dm2;
  int point;
  int i = 0;
  while( i < argc ) {
    if( argv[i] == std::string("--ue42") ) {
      Ue42 = atof(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--um42") ) {
      Um42 = atof(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--ut42") ) {
      Ut42 = atof(argv[i+1]);
      i += 2;
    } else if( argv[i] == std::string("--dm2") ) {
      dm2 = atof(argv[i+1]);
      i += 2;
    }
      else i += 1;
  }

  TFile *ftP_m = new TFile(Form("fitPara_m.root"),"READ");
  TFile *ftP_e = new TFile(Form("fitPara_e.root"),"READ");

  TTree *tree_m = (TTree*)ftP_m->Get("pardir");
  TTree *tree_e = (TTree*)ftP_e->Get("pardir");

  double para0_m,para1_m,para2_m,para3_m,para4_m,para5_m,norm_m,fitPara_m[29][7];
  double para0_e,para1_e,para2_e,para3_e,para4_e,para5_e,norm_e,fitPara_e[29][7];

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


  TFile *f = new TFile(Form("LE_1112_2.root"),"READ");
  TH2D *LvsE_e = (TH2D*)f->Get("h_e");
  TH2D *LvsE_m = (TH2D*)f->Get("h_m");
  TH1D *LEdep_e[29], *LEdep_m[29];
  for(int i = 0; i < 29; i++) {
    LEdep_e[i] = (TH1D*)LvsE_e->ProjectionY(Form("LE_e_bin%d",i+1),i+1,i+1);
    LEdep_e[i]->Scale(1./LEdep_e[i]->Integral("width"));

    LEdep_m[i] = (TH1D*)LvsE_m->ProjectionY(Form("LE_m_bin%d",i+1),i+1,i+1);
    LEdep_m[i]->Scale(1./LEdep_m[i]->Integral("width"));
  }


  char var[20] = "ElepReco";
  int para, nuCut;
  para = 2;
  nuCut = 4;

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


  TFile *scales_f = new TFile(Form("FC%s.root",FCsample.c_str()),"READ");
  TMatrixD *scales = (TMatrixD*)scales_f->Get("hscales");

  CovarianceMatrix inputCov;

  for (int i = 0; i < nbins; ++i) {
    inputCov.FCWeights[i] = (*scales)(u,i);
    // if(print) std::cout << i << "\t" << (*scales)(u,i) << "\n";
    // std::cout << fl_cov->GetBinContent(i+1, i+1) << "\t" << det_cov->GetBinContent(i+1, i+1) << "\t" << sig_cov->GetBinContent(i+1, i+1) << "\n";
    for (int j = 0; j < nbins; ++j) {
      inputCov.FluxMx[i][j]         = fl_cov->GetBinContent(i+1, j+1);
      inputCov.CrossSectionMx[i][j] = sig_cov->GetBinContent(i+1, j+1);
      inputCov.DetectorMx[i][j]     = det_cov->GetBinContent(i+1, j+1);
    }

    // std::cout << inputCov.FluxMx[i][i] << "\t" << inputCov.DetectorMx[i][i] << "\t" << inputCov.CrossSectionMx[i][i] << "\n";
  }

  tf.setEnergyBins( energy_bins );
  tf.setCovmtr( fl_bins, sig_bins, det_bins );

  double par_tgt[4], par_no[4];
  par_tgt[0] = Ue42;
  par_tgt[1] = Um42;
  par_tgt[2] = Ut42;
  par_tgt[3] = dm2;

  for(int ii = 0; ii < 4; ii++) {
    par_no[ii] = 0.;
  }

  tf.setPara( var, nuCut, fitPara_m, fitPara_e );

  tf.getTarget( par_tgt );
  
  double nochi2 = tf.bfChi2(par_no[0], par_no[1], par_no[2], par_no[3]);

  //std::cout << nochi2 << "\n";


  ofstream myfile;
  myfile.open("output.txt");
  //myfile << Ue42 << "\t" << Um42 << "\t" << Ut42 << "\t" << dm2 << "\t" << bf_Ue42 << "\t" << bf_Um42 << "\t" << bf_Ut42 << "\t" << bf_dm2 << "\t" << bf_chi2 << "\t" << nochi2 << "\n";
  myfile << Ue42 << "\t" << Um42 << "\t" << Ut42 << "\t" << dm2 << "\t" << nochi2 << "\n";
  myfile.close();


  return(0);

}


