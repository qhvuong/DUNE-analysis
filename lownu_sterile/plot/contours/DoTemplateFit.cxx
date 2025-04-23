#include "TemplateFitter.cxx"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TStyle.h"
#include <TRandom.h>
#include <list>
#include <TMath.h>
#include "DUNEStyle.h"

int main()
{
  const char data_path[] = "/exp/dune/app/users/qvuong/data/lownu";

  TFile *ftP_m = new TFile(Form("%s/LEdep/fitPara_m.root", data_path),"READ");
  TFile *ftP_e = new TFile(Form("%s/LEdep/fitPara_e.root", data_path),"READ");

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

  TFile *f = new TFile(Form("%s/LEdep/LE_1112_2.root", data_path),"READ");
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
  nuCut = 3;

  TFile *CC_f  = new TFile(Form("%s/input_dfiles/CC_output_56bins.root",data_path),"READ");
  TFile *nue_f = new TFile(Form("%s/input_dfiles/nue_output_8bins.root",data_path),"READ");

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

  TFile *f_fl = new TFile(Form("%s/flux_covmtr/flux_covmtr%d_120.root",data_path,nuCut),"READ");
  TH2D *fl_cov = (TH2D*)f_fl->Get("hcv");
  TFile *f_sig = new TFile(Form("%s/xS_covmtr/total_sigmtr%d_5sig_120.root",data_path,nuCut), "READ");
  TH2D *sig_cov = (TH2D*)f_sig->Get("hcv");

  double sig_bins[nbins+1][nbins+1], fl_bins[nbins+1][nbins+1];
  for(int i=0; i<nbins; i++) {
    for(int j=0; j<nbins; j++) {
      fl_bins[i][j]    = fl_cov->GetBinContent(i+1, j+1);
      sig_bins[i][j]   = sig_cov->GetBinContent(i+1, j+1);
    }
  }

  tf.setEnergyBins( energy_bins );
  tf.setCovmtr( fl_bins, sig_bins );

  double seed_set[4][3], chi2_LowerLimit = 5E-3;
  seed_set[0][0] = 0.;
  seed_set[0][1] = 0.;
  seed_set[0][2] = 0.;
  seed_set[1][0] = 0.16;
  seed_set[1][1] = 0.24;
  seed_set[1][2] = 0.;
  seed_set[2][0] = 0.;
  seed_set[2][1] = 0.;
  seed_set[2][2] = 100.;
  seed_set[3][0] = 0.16;
  seed_set[3][1] = 0.24;
  seed_set[3][2] = 100.;

  double oscpar[3], oscpar_max[3], oscpar_min[3];
  int N = 2;

  oscpar_max[0] = 0.3; oscpar_min[0] = 1E-4;	//Uee2 range
  oscpar_max[1] = 0.3; oscpar_min[1] = 1E-4;	//Umm2 range
  oscpar_max[2] = 6.0; oscpar_min[2] = 6.0;	//dm2 range

  double log_min = TMath::Log10(oscpar_min[0]);
  double log_max = TMath::Log10(oscpar_max[0]);
  double binWidth = (log_max-log_min)/N;

  std::vector<double> binEdges(N + 1, 0);
  for (int i = 0; i <= N; ++i) {
    binEdges[i] = TMath::Power(10, log_min + i * binWidth);
  }

  for(int i=0; i<N; i++){
    double x = (binEdges[i] + binEdges[i+1])/2.;
    std::cout << x << "\n";}

/*
  TH2D *h = new TH2D("h","",N,&binEdges[0],N,&binEdges[0]);
      
  for(int j = 1; j <= N; j++) {
    oscpar[1] = h->GetXaxis()->GetBinCenter(j);
    for(int k = 1; k <= N; k++) {
      oscpar[0] = h->GetYaxis()->GetBinCenter(k);
      oscpar[2] = oscpar_min[2];
      std::cout << "\noscpar \t" << oscpar[0] << "\t" << oscpar[1] << "\t" << oscpar[2] << "\n";
      //if(j*k%50 == 0) printf("%f percent...\n", ((j-1)*N+k)*100./(N*N) );

      tf.setPara( var, nuCut, fitPara_m, fitPara_e );
      tf.getTarget( oscpar );

      double bf_dm2, bf_Uee2, bf_Umm2, bf_chi2=1E9, par_bf[3], seed[3];
  
      for(int run=0; run<4; run++){
	seed[0] = seed_set[run][0];
	seed[1] = seed_set[run][1];
	seed[2] = seed_set[run][2];
	std::cout << "\ncoarse\t" << run+1 << "\t" << seed[0] << "\t" << seed[1] << "\t" << seed[2] << "\n";
	bool isOK = tf.doFitCoarse( seed, bf_Uee2, bf_Umm2 , bf_dm2);
	double c2 = tf.bfChi2( bf_Uee2, bf_Umm2 , bf_dm2 );
	printf( "COARSE nue Best-fit Uee2 = %f, Umm2 = %f, dm2 = %f, chi2 = %f\n", bf_Uee2, bf_Umm2, bf_dm2, c2);
        if(bf_chi2>c2) {
	  bf_chi2 = c2;
	  par_bf[0] = bf_Uee2;
	  par_bf[1] = bf_Umm2;
	  par_bf[2] = bf_dm2;
        }
        if(bf_chi2<chi2_LowerLimit) break;
      }  
      printf( "FINAL COARSE nue Best-fit Uee2 = %f, Umm2 = %f, dm2 = %f, chi2 = %f\n", par_bf[0], par_bf[1], par_bf[2], bf_chi2);

      seed[0] = par_bf[0];
      seed[1] = par_bf[1];
      seed[2] = 0.;
     
      if(bf_chi2>chi2_LowerLimit){
        bool fitFine = tf.doFitFine1( par_bf, bf_Uee2, bf_Umm2 , bf_dm2);
        double c2 = tf.bfChi2( bf_Uee2, bf_Umm2 , bf_dm2 );
        if(bf_chi2 > c2) {
	  bf_chi2 = c2;
	  par_bf[0] = bf_Uee2;
	  par_bf[1] = bf_Umm2;
	  par_bf[2] = bf_dm2;
	}

      do{  
        if(bf_chi2<chi2_LowerLimit) break;
	std::cout << "\n FINE \t" << seed[0] << "\t" << seed[1] << "\t" << seed[2] << "\n";
	bool isOK1 = tf.doFitFine1( seed, bf_Uee2, bf_Umm2 , bf_dm2);
	double fine_chi2 = tf.bfChi2( bf_Uee2, bf_Umm2 , bf_dm2 );
	printf( "FINE nue Best-fit Uee2 = %f, Umm2 = %f, dm2 = %f, chi2 = %f\n", bf_Uee2, bf_Umm2, bf_dm2, fine_chi2);
	if(bf_chi2 > fine_chi2) {
	  bf_chi2 = fine_chi2;
	  par_bf[0] = bf_Uee2;
	  par_bf[1] = bf_Umm2;
	  par_bf[2] = bf_dm2;
	}
        if(fine_chi2<chi2_LowerLimit) break;
	seed[2] = (1.0+seed[2])*2.; 
      } while(seed[2]<200.0);
      printf( "FINAL FINE1 nue Best-fit Uee2 = %f, Umm2 = %f, dm2 = %f, chi2 = %f\n", par_bf[0], par_bf[1], par_bf[2], bf_chi2);
      }     
 
      bool isOK2 = tf.doFitFine2( par_bf, bf_Uee2, bf_Umm2 , bf_dm2 );
      double fine_chi2 = tf.bfChi2( bf_Uee2, bf_Umm2 , bf_dm2 );
      printf( "FINE nue Best-fit Uee2 = %f, Umm2 = %f, dm2 = %f, chi2 = %f\n", bf_Uee2, bf_Umm2, bf_dm2, fine_chi2);
      if(bf_chi2 > fine_chi2) {
        bf_chi2 = fine_chi2;
        par_bf[0] = bf_Uee2;
        par_bf[1] = bf_Umm2;
        par_bf[2] = bf_dm2;
      }
      

      double nopar[3];
      nopar[0] = nopar[1] = nopar[2] = 0.;

      double nochi2 = tf.bfChi2( nopar[0], nopar[1], nopar[2] );
      double dchi2 = nochi2 - bf_chi2;
      printf( "FINAL FINE2 nue Best-fit Uee2 = %f, Umm2 = %f, dm2 = %f, chi2 = %f, nochi2 = %f\n", par_bf[0], par_bf[1], par_bf[2], bf_chi2, nochi2);

      h->Fill(oscpar[1],oscpar[0],dchi2);

    }
    }

    gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
    gStyle->SetNumberContours(999);

    h->SetStats(0);
    //h->SetMaximum(8);
    //h->SetMinimum(0.1);
    h->SetTitle(Form("Sensitivity Contour (#Deltam^{2}=%.2f)",oscpar[2]));
    h->GetYaxis()->SetTitle("U_{e4}^{2}");
    h->GetXaxis()->SetTitle("U_{#mu4}^{2}");

    TCanvas *c0 = new TCanvas("c0","",600,600);
    c0->SetLogx();
    c0->SetLogy();
    c0->SetLogz();
    c0->SetGrid();
    h->Draw("colz");
    c0->SaveAs(Form("contour%d.png",N));
*/

}



