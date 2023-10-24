#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TStyle.h"
#include <TRandom.h>
#include <list>

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
  int para, nuCut;
  para = 2;
  nuCut = 0;

  TFile *CC_f  = new TFile("/dune/app/users/qvuong/data/lownu/CC_output.root","READ");
  TFile *nue_f = new TFile("/dune/app/users/qvuong/data/lownu/nue_output.root","READ");

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

  TFile *f_fl = new TFile(Form("../../../flux_covmtr/flux_covmtr%d_10000.root",nuCut),"READ");
  TH2D *fl_cov = (TH2D*)f_fl->Get("hcv");
  TFile *f_sys = new TFile(Form("../../../xS_covmtr/total_sigmtr%d_5sig.root",nuCut), "READ");
  TH2D *sys_cov = (TH2D*)f_sys->Get("hcv");

  TCanvas *c = new TCanvas("c","",800,600);
  fl_cov->Draw("colz");
  c->SaveAs("fl.png");
  sys_cov->Draw("colz");
  c->SaveAs("sys.png");

  double cov_bins[nbins+1][nbins+1], sys_bins[nbins+1][nbins+1], fl_bins[nbins+1][nbins+1];
  for(int i=0; i<nbins; i++) {
    for(int j=0; j<nbins; j++) {
      fl_bins[i][j]  = fl_cov->GetBinContent(i+1, j+1);
      sys_bins[i][j] = sys_cov->GetBinContent(i+1, j+1);
      cov_bins[i][j] = fl_bins[i][j] + sys_bins[i][j];
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
  double bfc2 = tf.bfChi2 ( par );

  std::cout << "bfchi2: " << par[0] << "\t" << par[1] << "\t" << par[2] << "\t" << bfc2 << "\n";

/*
  double oscpar[3], oscpar_max[3], oscpar_min[3], stepsize[3], seed[3];
  int N = 20;

  oscpar_max[0] = 1E-4; oscpar_min[0] = 1E-4;	//Uee2 range
  oscpar_max[1] = 1E-2; oscpar_min[1] = 1E-4;	//Umm2 range
  oscpar_max[2] = 10.0; oscpar_min[2] = 1.0;	//dm2 range

  TH2D *h0 = new TH2D("h0","",N,oscpar_min[1],oscpar_max[1],N,oscpar_min[2],oscpar_max[2]);

  {
  for(int j = 1; j <= N; j++) {
    oscpar[1] = h0->GetXaxis()->GetBinCenter(j);
    for(int k = 1; k <= N; k++) {
      oscpar[2] = h0->GetYaxis()->GetBinCenter(k);

      oscpar[0] = oscpar_min[0];

      printf("%f percent...\n", ((j-1)*N+k)*100./(N*N) );
    
      for(int ii = 0; ii < 3; ii++) {
        seed[ii] = oscpar[ii]; 
      }

      tf.setPara( var, oscpar, nuCut, EvCut, seed, fitPara_m, fitPara_e );

      tf.getTarget( oscpar );

      double bf_dm2, bf_Uee2, bf_Umm2;
      double par[3];
      bool isOK = tf.doFit( bf_Uee2, bf_Umm2 , bf_dm2);
      par[0] = bf_Uee2;
      par[1] = bf_Umm2;
      par[2] = bf_dm2;
      double chi2 = tf.bfChi2( par );
      double nochi2 = tf.noChi2( par );
      double sen = sqrt(std::abs(nochi2 - chi2));
      printf( "nue Best-fit Uee2 = %.2f, Umm2 = %.2f, dm2 = %.2f, chi2 = %f, noChi2 = %f, sen = %f\n", bf_Uee2, bf_Umm2, bf_dm2, chi2, nochi2, sen );

      //double s2mue2 = (par[0]*par[0]) * (par[1]*par[1]);

      //h1->Fill(s2mue2,oscpar[2],chi2);
      h0->Fill(oscpar[1],oscpar[2],sen);
    }
    }


    gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
    gStyle->SetNumberContours(999);

    TCanvas *c0 = new TCanvas("c0","",600,600);
    c0->SetLogx();
    c0->SetLogz();
    h0->SetStats(0);
    h0->SetMaximum(8);
    h0->SetMinimum(0.1);
    h0->Draw("colz");
    h0->SetTitle(Form("Sensitivity Contour (Uee2=%f)",oscpar[0]));
    h0->GetXaxis()->SetTitle("Umm2");
    h0->GetYaxis()->SetTitle("dm2");
    c0->SaveAs(Form("contour_Logx.png"));

    TCanvas *c1 = new TCanvas("c1","",600,600);
    //c0->SetLogz();
    h0->SetStats(0);
    h0->SetMaximum(8);
    h0->SetMinimum(0.1);
    h0->Draw("colz");
    h0->SetTitle(Form("Sensitivity Contour (Uee2=%f)",oscpar[0]));
    h0->GetXaxis()->SetTitle("Umm2");
    h0->GetYaxis()->SetTitle("dm2");
    c1->SaveAs(Form("contour.png"));

    gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
    gStyle->SetNumberContours(999);

    TFile* out = new TFile("contour.root","RECREATE");
    h0->Write();
    out->Close();
  }
*/
}



