#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TStyle.h"
#include <TRandom.h>
#include <list>
#include <iostream>
#include <fstream>

int main()
{
  TFile *ftP_m = new TFile("/dune/app/users/qvuong/lownu/LEdep/fitPara_m.root","READ");
  TFile *ftP_e = new TFile("/dune/app/users/qvuong/lownu/LEdep/fitPara_e.root","READ");

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


  TFile *f = new TFile("/dune/app/users/qvuong/lownu/LEdep/LE_1112_2.root","READ");
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
  
  //for(int N = 300;N<901;N+=100){
  int N = 10000;
  for(int para=1; para<3; para++){
  for(int nuCut=0; nuCut<4; nuCut+=3){

  //int para, nuCut;
  //para = 2;
  //nuCut = 0;

  TFile *CC_f  = new TFile("/dune/app/users/qvuong/data/lownu/CC_output.root","READ");
  TFile *nue_f = new TFile("/dune/app/users/qvuong/data/lownu/nue_output.root","READ");

  TH2D* CC_hm    = (TH2D*)CC_f->Get(Form("m_h%sVsEv%d",var,nuCut));
  TH2D* CC_hm_nc = (TH2D*)CC_f->Get(Form("nc_m_h%sVsEv%d",var,nuCut));
  TH2D* CC_he    = (TH2D*)CC_f->Get(Form("e_h%sVsEv%d",var,nuCut));
  TH2D* nue_hm   = (TH2D*)nue_f->Get(Form("m_h%sVsEv0",var));
  TH2D* nue_hm_w = (TH2D*)nue_f->Get(Form("m_h%sVsEv0_w",var));
  TH2D* nue_he   = (TH2D*)nue_f->Get(Form("e_h%sVsEv0",var));
  TH2D* nue_he_w = (TH2D*)nue_f->Get(Form("e_h%sVsEv0_w",var));

  std::cout << CC_hm->GetNbinsX() << "\t" << CC_hm->GetNbinsY() << "\n";
  std::cout << CC_he->GetNbinsX() << "\t" << CC_he->GetNbinsY() << "\n";
  std::cout << nue_hm->GetNbinsX() << "\t" << nue_hm->GetNbinsY() << "\n";


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

  TemplateFitter tf( CC_templates_m, CC_templates_m_nc, CC_templates_e, nue_templates_m, nue_templates_m_w, nue_templates_e, nue_templates_e_w, LEdep_m, LEdep_e );


  std::cout << CC_templates_m[4]->GetNbinsX() << "\n";
  std::cout << CC_templates_e[4]->GetNbinsX() << "\n";
  std::cout << nue_templates_m[4]->GetNbinsX() << "\n";

  double energy_bins[nbins_Ev+1];
  for( int b = 0; b <= nbins_Ev; ++b ) {
    energy_bins[b] = CC_he->GetXaxis()->GetBinLowEdge(b+1);
  } 

  TFile *cov_f = new TFile(Form("../FC_stat_fl_3sig_%d_%d.root",nuCut,N),"READ");
  TH2D *cov = (TH2D*)cov_f->Get("cv");

  double cov_bins[nbins+1][nbins+1];
  for(int i=0; i<nbins; i++) {
    for(int j=0; j<nbins; j++) {
      cov_bins[i][j]  = cov->GetBinContent(i+1, j+1);
    }
  }

  tf.setEnergyBins( energy_bins );
  tf.setCovmtr( cov_bins );

  double oscpar[3], seed[3];
  if( para == 1 ) {
  oscpar[0] = 0.01;
  oscpar[1] = 0.0016;
  oscpar[2] = 1.3;
  }
  if( para == 2 ) {
  oscpar[0] = 0.04;
  oscpar[1] = 0.01;
  oscpar[2] = 6.0;
  }
  for(int ii = 0; ii < 3; ii++) {
    seed[ii] = oscpar[ii]; }

  tf.setPara( var, oscpar, para, nuCut, seed, fitPara_m, fitPara_e );

  tf.getTarget( oscpar );

  double bf_dm2, bf_Uee2, bf_Umm2;
  double par[3];
  bool isOK = tf.doFit( bf_Uee2, bf_Umm2 , bf_dm2);
  par[0] = bf_Uee2;
  par[1] = bf_Umm2;
  par[2] = bf_dm2;
  //double chi2 = tf.bfChi2 ( par );
  tf.Draw( par );
/*
  std::ofstream out;
  out.open(Form("FC_chi2_%d%d_%d.txt",para,nuCut,N));
  out << chi2 << "\n";
  out.close();
*/
  }
  }
  //}
}



