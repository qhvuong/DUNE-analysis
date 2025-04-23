#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TStyle.h"
#include <TRandom.h>
#include <list>
#include <iostream>
#include <fstream>

using namespace std;

int main(int argc, char const *argv[])
//int main()
{

  int u;
  int i = 0;
  while( i < argc ) {
    if( argv[i] == std::string("--u") ) {
      u = atof(argv[i+1]);
      i += 2;
    }

    else i += 1;
  }
  
  TFile *ftP_m = new TFile("fitPara_m.root","READ");
  TFile *ftP_e = new TFile("fitPara_e.root","READ");

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


  TFile *f = new TFile("LE_1112_2.root","READ");
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
  
  int para = 2;
  int nuCut = 3;

  TFile *CC_f  = new TFile("CC_output_56bins.root","READ");
  TFile *nue_f = new TFile("nue_output_8bins.root","READ");

  TH2D* CC_hm    = (TH2D*)CC_f->Get(Form("m_h%sVsEv%d",var,nuCut));
  TH2D* CC_hm_nc = (TH2D*)CC_f->Get(Form("nc_m_h%sVsEv%d",var,nuCut));
  TH2D* CC_he    = (TH2D*)CC_f->Get(Form("e_h%sVsEv%d",var,nuCut));
  TH2D* nue_hm   = (TH2D*)nue_f->Get(Form("m_h%sVsEv0",var));
  TH2D* nue_hm_w = (TH2D*)nue_f->Get(Form("m_h%sVsEv0_w",var));
  TH2D* nue_he   = (TH2D*)nue_f->Get(Form("e_h%sVsEv0",var));
  TH2D* nue_he_w = (TH2D*)nue_f->Get(Form("e_h%sVsEv0_w",var));


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

  double energy_bins[nbins_Ev+1];
  for( int b = 0; b <= nbins_Ev; ++b ) {
    energy_bins[b] = CC_he->GetXaxis()->GetBinLowEdge(b+1);
  } 


  TFile *f_fl = new TFile(Form("flux_covmtr%d_120.root",nuCut),"READ");
  TH2D *fl_cov = (TH2D*)f_fl->Get("hcvfr");
  TFile *f_sys = new TFile(Form("total_sigmtr%d_5sig_120.root",nuCut), "READ");
  TH2D *sys_cov = (TH2D*)f_sys->Get("frhcv");
  //TFile *scales_f = new TFile(Form("FCstat%d_120_10000.root",nuCut),"READ");
  TFile *scales_f = new TFile(Form("FCflux%d_120_10000.root",nuCut),"READ");
  TMatrixD *scales = (TMatrixD*)scales_f->Get("hscales");

  TCanvas *c = new TCanvas("c","",800,300);
  c->Divide(2,1);
  c->cd(1);
  fl_cov->Draw("colz");
  c->cd(2);
  sys_cov->Draw("colz");
  c->SaveAs("cov.png");
  

  TH1D *hist_scales = new TH1D("hist_scales","",100,0.0,2.0);
/*
  for(int u=0; u<10000; u++){
*/
  double wgt[nbins+1];
  double sys_bins[nbins+1][nbins+1], fl_bins[nbins+1][nbins+1];
  for(int i=0; i<nbins; i++) {
    wgt[i] = (*scales)(u,i);
    for(int j=0; j<nbins; j++) {
      fl_bins[i][j]    = fl_cov->GetBinContent(i+1, j+1);
      sys_bins[i][j]   = sys_cov->GetBinContent(i+1, j+1);
    }
  }
  

  tf.setEnergyBins( energy_bins );
  tf.setCovmtr( fl_bins, sys_bins, wgt );

  tf.setPara( var, nuCut, fitPara_m, fitPara_e, u );

  double par_tgt[3];
  par_tgt[0] = par_tgt[1] = par_tgt[2] = 0.;

  double sc = tf.getTarget( par_tgt );
  //std::cout << "ratio=\t" << sc << "\n";
/*
  hist_scales->Fill(sc);

  if(u%100==0) std::cout << u*100./10000 << "percent...\n";
  }

  TCanvas *c = new TCanvas("c","",800,600);
  c->SetGrid();
  hist_scales->Draw();
  hist_scales->SetTitle("stats only throw");
  hist_scales->GetXaxis()->SetTitle("throw/original area ratio");
  hist_scales->GetYaxis()->SetTitle("entries");
  c->SaveAs("hist_scales_stats.png");
*/

  double bf_dm2, bf_Uee2, bf_Umm2, nochi2, bfchi2;

  const char data_path[] = "/pnfs/dune/scratch/users/qvuong/output";
  const char out_path[] = "/exp/dune/app/users/qvuong/data/lownu/Feldman_Cousins";

  ifstream fstat_true(Form("%s/FCstat_new/seeding_true/output_%d.txt",data_path,u));
  ifstream fstat_alg(Form("%s/FCstat_new/seeding_algorithm/output_%d.txt",data_path,u));
  ifstream fflux_true(Form("%s/FCflux/seeding_true/output_%d.txt",data_path,u));
  ifstream fflux_alg(Form("%s/FCflux/seeding_algorithm/output_%d.txt",data_path,u));


  fstat_true >> bf_Uee2 >> bf_Umm2 >> bf_dm2 >> bfchi2 >> nochi2; 
  //fstat_alg >> bf_Uee2 >> bf_Umm2 >> bf_dm2 >> bfchi2 >> nochi2; 
  //fflux_alg >> bf_Uee2 >> bf_Umm2 >> bf_dm2 >> bfchi2 >> nochi2; 
  //fflux_true >> bf_Uee2 >> bf_Umm2 >> bf_dm2 >> bfchi2 >> nochi2; 
  
  double bf_chi2 = tf.bfChi2( bf_Uee2, bf_Umm2 , bf_dm2 );
  //double bf_chi2 = tf.bfChi2( 0.04, 0.01 , 6. );
  printf( "FINE nue Best-fit Uee2 = %f, Umm2 = %f, dm2 = %f, nochi2 = %f, bfchi2 = %f, bf_chi2 = %f\n", bf_Uee2, bf_Umm2, bf_dm2, nochi2, bfchi2, bf_chi2);



}



