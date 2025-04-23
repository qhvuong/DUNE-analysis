#include "DUNEStyle.h"


static const int N = 1; // number of universes
static const int nbins_CC = 56;

static const int nbins_nue = 8;
static const int nbins = 2*nbins_CC + nbins_nue;

static const char data_path[] = "/exp/dune/app/users/qvuong/data/lownu";


TMatrixD ECovars_m2_m2   ( nbins_CC, nbins_CC );
TMatrixD ECovars_m2_e2   ( nbins_CC, nbins_CC );
TMatrixD ECovars_m2_nue  ( nbins_CC, nbins_nue );

TMatrixD ECovars_e2_m2   ( nbins_CC, nbins_CC );
TMatrixD ECovars_e2_e2   ( nbins_CC, nbins_CC );
TMatrixD ECovars_e2_nue  ( nbins_CC, nbins_nue );

TMatrixD ECovars_nue_m2  ( nbins_nue, nbins_CC );
TMatrixD ECovars_nue_e2  ( nbins_nue, nbins_CC );
TMatrixD ECovars_nue_nue ( nbins_nue, nbins_nue );


TMatrixD frECovars_m2_m2   ( nbins_CC, nbins_CC );
TMatrixD frECovars_m2_e2   ( nbins_CC, nbins_CC );
TMatrixD frECovars_m2_nue  ( nbins_CC, nbins_nue );

TMatrixD frECovars_e2_m2   ( nbins_CC, nbins_CC );
TMatrixD frECovars_e2_e2   ( nbins_CC, nbins_CC );
TMatrixD frECovars_e2_nue  ( nbins_CC, nbins_nue );

TMatrixD frECovars_nue_m2  ( nbins_nue, nbins_CC );
TMatrixD frECovars_nue_e2  ( nbins_nue, nbins_CC );
TMatrixD frECovars_nue_nue ( nbins_nue, nbins_nue );



TMatrixD ECovars_m4_m4   ( nbins_CC, nbins_CC );
TMatrixD ECovars_m4_e4   ( nbins_CC, nbins_CC );
TMatrixD ECovars_m4_nue  ( nbins_CC, nbins_nue );

TMatrixD ECovars_e4_m4   ( nbins_CC, nbins_CC );
TMatrixD ECovars_e4_e4   ( nbins_CC, nbins_CC );
TMatrixD ECovars_e4_nue  ( nbins_CC, nbins_nue );

TMatrixD ECovars_nue_m4  ( nbins_nue, nbins_CC );
TMatrixD ECovars_nue_e4  ( nbins_nue, nbins_CC );


TMatrixD frECovars_m4_m4   ( nbins_CC, nbins_CC );
TMatrixD frECovars_m4_e4   ( nbins_CC, nbins_CC );
TMatrixD frECovars_m4_nue  ( nbins_CC, nbins_nue );

TMatrixD frECovars_e4_m4   ( nbins_CC, nbins_CC );
TMatrixD frECovars_e4_e4   ( nbins_CC, nbins_CC );
TMatrixD frECovars_e4_nue  ( nbins_CC, nbins_nue );

TMatrixD frECovars_nue_m4  ( nbins_nue, nbins_CC );
TMatrixD frECovars_nue_e4  ( nbins_nue, nbins_CC );



TMatrixD ECovars_mm   ( nbins_CC, nbins_CC );
TMatrixD ECovars_me   ( nbins_CC, nbins_CC );
TMatrixD ECovars_mnue ( nbins_CC, nbins_nue );

TMatrixD ECovars_em   ( nbins_CC, nbins_CC );
TMatrixD ECovars_ee   ( nbins_CC, nbins_CC );
TMatrixD ECovars_enue ( nbins_CC, nbins_nue );

TMatrixD ECovars_nuem   ( nbins_nue, nbins_CC );
TMatrixD ECovars_nuee   ( nbins_nue, nbins_CC );
TMatrixD ECovars_nuenue ( nbins_nue, nbins_nue );


TMatrixD frECovars_mm   ( nbins_CC, nbins_CC );
TMatrixD frECovars_me   ( nbins_CC, nbins_CC );
TMatrixD frECovars_mnue ( nbins_CC, nbins_nue );

TMatrixD frECovars_em   ( nbins_CC, nbins_CC );
TMatrixD frECovars_ee   ( nbins_CC, nbins_CC );
TMatrixD frECovars_enue ( nbins_CC, nbins_nue );

TMatrixD frECovars_nuem   ( nbins_nue, nbins_CC );
TMatrixD frECovars_nuee   ( nbins_nue, nbins_CC );
TMatrixD frECovars_nuenue ( nbins_nue, nbins_nue );


TMatrixD ECovars2 ( nbins, nbins );
TMatrixD ECovars  ( nbins, nbins );
TMatrixD ECovars4 ( nbins, nbins );
TMatrixD frECovars2 ( nbins, nbins );
TMatrixD frECovars  ( nbins, nbins );
TMatrixD frECovars4 ( nbins, nbins );




void covariance_pn(TH1D* hnom, TH1D* hp, TH1D *hn, TH2D* cov, TH2D* frCov, TCanvas* c, const char *name){
  int nbins = hnom->GetNbinsX();
  //std::cout << "nbins = " << nbins << "\n";

  for(int ii=0; ii<nbins; ii++){
    double n_i = hnom->GetBinContent(ii+1);
    double ps_i = hp->GetBinContent(ii+1);  //positive shift histogram (hpmu or hpel)
    double ns_i = hn->GetBinContent(ii+1);  //negative shift histogram (hnmu or hnel)

    double pdelta_i = ps_i - n_i;
    double ndelta_i = ns_i - n_i;

    for(int jj=0; jj<nbins; jj++){
      double n_j = hnom->GetBinContent(jj+1);
      double ps_j = hp->GetBinContent(jj+1);
      double ns_j = hn->GetBinContent(jj+1);

      double pdelta_j = ps_j - n_j;
      double ndelta_j = ns_j - n_j;

      double pfcov_ij = abs(pdelta_i) * abs(pdelta_j);  //positive shift covariance
      double nfcov_ij = abs(ndelta_i) * abs(ndelta_j);  //negative shift covariance
      double afcov_ij = (pfcov_ij+nfcov_ij) / 2;   //average covariance which ends up in the histograms
      cov->SetBinContent(ii+1, jj+1, afcov_ij);
      double cov_ij = afcov_ij/(n_i*n_j);
      frCov->SetBinContent(ii+1, jj+1, cov_ij); 
    }
  }
  
 
  //makes lines to divide matrix to show the three different samples
    double x = cov->GetXaxis()->GetBinUpEdge(nbins_CC);
    double y1 = cov->GetYaxis()->GetBinLowEdge(1);
    double y2 = cov->GetYaxis()->GetBinUpEdge(nbins);

    double y = cov->GetYaxis()->GetBinUpEdge(nbins_CC);
    double x1 = cov->GetXaxis()->GetBinLowEdge(1);
    double x2 = cov->GetXaxis()->GetBinUpEdge(nbins);

    double w = cov->GetXaxis()->GetBinUpEdge(2*nbins_CC);
    double v1 = cov->GetYaxis()->GetBinLowEdge(1);
    double v2 = cov->GetYaxis()->GetBinUpEdge(nbins);

    double v = cov->GetYaxis()->GetBinUpEdge(2*nbins_CC);
    double w1 = cov->GetXaxis()->GetBinLowEdge(1);
    double w2 = cov->GetXaxis()->GetBinUpEdge(nbins);

    TLine *l1 = new TLine(x, y1, x, y2);  //splits muCC from the rest in x axis 
    TLine *l2 = new TLine(x1, y, x2, y);  //splits muCC from the rest in y axis
    TLine *l3 = new TLine(w, v1, w, v2);  //splits nu+e from the rest in x axis
    TLine *l4 = new TLine(w1, v, w2, v);  //splits nu+e from the rest in y axis

    cov->GetXaxis()->SetTitle("Bin number");
    cov->GetYaxis()->SetTitle("Bin number");
    frCov->GetXaxis()->SetTitle("Bin number");
    frCov->GetYaxis()->SetTitle("Bin number");

    c->Clear();
    frCov->SetMinimum(1E-5);
    frCov->SetMaximum(1.);
    frCov->SetStats(0);
    frCov->Draw("colz");
    l1->SetLineWidth(2);
    l1->SetLineColor(kBlack);
    l1->Draw("same");
    l2->SetLineWidth(2);
    l2->SetLineColor(kBlack);
    l2->Draw("same");
    l3->SetLineWidth(2);
    l3->SetLineColor(kBlack);
    l3->Draw("same");
    l4->SetLineWidth(2);
    l4->SetLineColor(kBlack);
    l4->Draw("same");
    c->SaveAs(Form("%s/uncertainties/xS_covmtr/xScov/fractional_%s.png",data_path,name));


}





void covariance(TH1D* hnom, TH1D* h, TH2D* cov, TH2D* frCov, TCanvas* c, const char *name){
  int nbins = hnom->GetNbinsX();
  //std::cout << "nbins = " << nbins << "\n";

  for(int ii=0; ii<nbins; ii++){
    double n_i = hnom->GetBinContent(ii+1);
    double s_i = h->GetBinContent(ii+1);  // shift histogram (hpmu or hpel)
    double delta_i = s_i - n_i;

    for(int jj=0; jj<nbins; jj++){
      double n_j = hnom->GetBinContent(jj+1);
      double s_j = h->GetBinContent(jj+1);
      double delta_j = s_j - n_j;
      double cov_ij = abs(delta_i) * abs(delta_j);  //positive shift covariance
      cov->SetBinContent(ii+1, jj+1, cov_ij);  //hmu_cov or hel_cov
      cov_ij = cov_ij/(n_i*n_j);
      frCov->SetBinContent(ii+1, jj+1, cov_ij); 
    }
  }
  

  //makes lines to divide matrix to show the three different samples
    double x = cov->GetXaxis()->GetBinUpEdge(nbins_CC);
    double y1 = cov->GetYaxis()->GetBinLowEdge(1);
    double y2 = cov->GetYaxis()->GetBinUpEdge(nbins);

    double y = cov->GetYaxis()->GetBinUpEdge(nbins_CC);
    double x1 = cov->GetXaxis()->GetBinLowEdge(1);
    double x2 = cov->GetXaxis()->GetBinUpEdge(nbins);

    double w = cov->GetXaxis()->GetBinUpEdge(2*nbins_CC);
    double v1 = cov->GetYaxis()->GetBinLowEdge(1);
    double v2 = cov->GetYaxis()->GetBinUpEdge(nbins);

    double v = cov->GetYaxis()->GetBinUpEdge(2*nbins_CC);
    double w1 = cov->GetXaxis()->GetBinLowEdge(1);
    double w2 = cov->GetXaxis()->GetBinUpEdge(nbins);

    TLine *l1 = new TLine(x, y1, x, y2);  //splits muCC from the rest in x axis 
    TLine *l2 = new TLine(x1, y, x2, y);  //splits muCC from the rest in y axis
    TLine *l3 = new TLine(w, v1, w, v2);  //splits nu+e from the rest in x axis
    TLine *l4 = new TLine(w1, v, w2, v);  //splits nu+e from the rest in y axis

    cov->GetXaxis()->SetTitle("Bin number");
    cov->GetYaxis()->SetTitle("Bin number");
    frCov->GetXaxis()->SetTitle("Bin number");
    frCov->GetYaxis()->SetTitle("Bin number");

    c->Clear();
    frCov->SetMinimum(1E-5);
    frCov->SetMaximum(1.);
    frCov->SetStats(0);
    frCov->Draw("colz");
    l1->SetLineWidth(2);
    l1->SetLineColor(kBlack);
    l1->Draw("same");
    l2->SetLineWidth(2);
    l2->SetLineColor(kBlack);
    l2->Draw("same");
    l3->SetLineWidth(2);
    l3->SetLineColor(kBlack);
    l3->Draw("same");
    l4->SetLineWidth(2);
    l4->SetLineColor(kBlack);
    l4->Draw("same");
    c->SaveAs(Form("%s/uncertainties/xS_covmtr/xScov/fractional_%s.png",data_path,name));
 
}






void ElepCov()
{

  TFile *f     = new TFile(Form("%s/input_dfiles/CCtest_output.root",data_path), "READ");
  TFile *f_nue = new TFile(Form("%s/input_dfiles/nue_output.root",data_path), "READ");

  const int N_nucut = 5;
  double nucut[N_nucut] = {50., 5., 1., 0.5, 0.3};

  TH1D *hn = new TH1D("hn","",nbins,0,nbins);
  TH1D *ho = new TH1D("ho","",nbins,0,nbins);
  TH1D *hp = new TH1D("hp","",nbins,0,nbins);
  
  TH1D *nue_n = (TH1D*)f_nue->Get("hSigma_n");
  TH1D *nue_o = (TH1D*)f_nue->Get("hElep");
  TH1D *nue_p = (TH1D*)f_nue->Get("hSigma_p");

  TCanvas *c = new TCanvas("c","",800,600);
  //gStyle->SetPalette(kColorPrintableOnGrey); 
  TColor::InvertPalette();
  c->SetLogz();

  std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};
  const char *name[] = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};

  //std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"};
  //const char *name[] = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"}; 

  //std::list <const char *> namelist = {"wgt_MaCCQE"};
  //const char *name[] = {"wgt_MaCCQE"}; 

  int N_wgt = namelist.size();

  std::cout << N_wgt << "\n";

  TH2D ***hcv_n = new TH2D**[N_nucut];
  TH2D ***hcv = new TH2D**[N_nucut];
  TH2D ***hcv_p = new TH2D**[N_nucut];

  TH2D ***frhcv_n = new TH2D**[N_nucut];
  TH2D ***frhcv = new TH2D**[N_nucut];
  TH2D ***frhcv_p = new TH2D**[N_nucut];

  for(int j=0; j<N_nucut; j++){
    hcv_n[j] = new TH2D*[N_wgt+1];
    hcv[j] = new TH2D*[N_wgt+1];
    hcv_p[j] = new TH2D*[N_wgt+1];

    frhcv_n[j] = new TH2D*[N_wgt+1];
    frhcv[j] = new TH2D*[N_wgt+1];
    frhcv_p[j] = new TH2D*[N_wgt+1];

    for(int iw=0; iw<N_wgt+1; iw++){
      if(iw<N_wgt){
        hcv_n[j][iw] = new TH2D(Form("%s_hcv%d_n",name[iw],j), Form("%s Covariance Matrix -1sigma (#nu < %.2f GeV)", name[iw], nucut[j]), nbins,0,nbins,nbins,0,nbins);
        hcv[j][iw] = new TH2D(Form("%s_hcv%d",name[iw],j), Form("%s Covariance Matrix (#nu < %.2f GeV)", name[iw], nucut[j]), nbins,0,nbins,nbins,0,nbins);
        hcv_p[j][iw] = new TH2D(Form("%s_hcv%d_p",name[iw],j), Form("%s Covariance Matrix +1sigma (#nu < %.2f GeV)", name[iw], nucut[j]), nbins,0,nbins,nbins,0,nbins);

        frhcv_n[j][iw] = new TH2D(Form("%s_frhcv%d_n",name[iw],j), Form("%s Fractional Covariance Matrix -1sigma (#nu < %.2f GeV)", name[iw], nucut[j]), nbins,0,nbins,nbins,0,nbins);
        frhcv[j][iw] = new TH2D(Form("%s_frhcv%d",name[iw],j), Form("%s Fractional Covariance Matrix (#nu < %.2f GeV)", name[iw], nucut[j]), nbins,0,nbins,nbins,0,nbins);
        frhcv_p[j][iw] = new TH2D(Form("%s_frhcv%d_p",name[iw],j), Form("%s Fractional Covariance Matrix +1sigma (#nu < %.2f GeV)", name[iw], nucut[j]), nbins,0,nbins,nbins,0,nbins);
      }

      else{
        hcv_n[j][iw] = new TH2D(Form("nue_hcv%d_n",j), Form("Nue-xS Covariance Matrix -1sigma (#nu < %.2f GeV)", nucut[j]), nbins,0,nbins,nbins,0,nbins);
        hcv[j][iw] = new TH2D(Form("nue_hcv%d",j), Form("Nue-xS Covariance Matrix (#nu < %.2f GeV)", nucut[j]), nbins,0,nbins,nbins,0,nbins);
        hcv_p[j][iw] = new TH2D(Form("nue_hcv%d_p",j), Form("Nue-xS Covariance Matrix +1sigma (#nu < %.2f GeV)", nucut[j]), nbins,0,nbins,nbins,0,nbins);

        frhcv_n[j][iw] = new TH2D(Form("nue_frhcv%d_n",j), Form("Nue-xS Fractional Covariance Matrix -1sigma (#nu < %.2f GeV)", nucut[j]), nbins,0,nbins,nbins,0,nbins);
        frhcv[j][iw] = new TH2D(Form("nue_frhcv%d",j), Form("Nue-xS Fractional Covariance Matrix (#nu < %.2f GeV)", nucut[j]), nbins,0,nbins,nbins,0,nbins);
        frhcv_p[j][iw] = new TH2D(Form("nue_frhcv%d_p",j), Form("Nue-xS Fractional Covariance Matrix +1sigma (#nu < %.2f GeV)", nucut[j]), nbins,0,nbins,nbins,0,nbins);
      }
    }
  }

  
  TH1D *mn = (TH1D*)f->Get(Form("wgt_MaCCQE_mElep0_n"));
  TH1D *mo = (TH1D*)f->Get(Form("wgt_MaCCQE_mElep0_o"));
  TH1D *mp = (TH1D*)f->Get(Form("wgt_MaCCQE_mElep0_p"));

  TH1D *en = (TH1D*)f->Get(Form("wgt_MaCCQE_eElep0_n"));
  TH1D *eo = (TH1D*)f->Get(Form("wgt_MaCCQE_eElep0_o"));
  TH1D *ep = (TH1D*)f->Get(Form("wgt_MaCCQE_eElep0_p"));


  for(int j=0; j<N_nucut; j++){
    if(j%2==0) continue;

    mo = (TH1D*)f->Get(Form("wgt_MaCCQE_mElep%d_o",j));
    eo = (TH1D*)f->Get(Form("wgt_MaCCQE_eElep%d_o",j));

    for(int iw=0; iw<N_wgt+1; iw++){
      if(iw%5==0) std::cout << "Nu cut " << j << ":\t Dials " << iw << "\n";

      if(iw<N_wgt){
    
        mn = (TH1D*)f->Get(Form("%s_mElep%d_n",name[iw],j));
        mp = (TH1D*)f->Get(Form("%s_mElep%d_p",name[iw],j));

        en = (TH1D*)f->Get(Form("%s_eElep%d_n",name[iw],j));
        ep = (TH1D*)f->Get(Form("%s_eElep%d_p",name[iw],j));
    
        for(int i=0; i<nbins; i++){
          if(i<nbins_CC){
            hn->SetBinContent(i+1, mn->GetBinContent(i+1));
            ho->SetBinContent(i+1, mo->GetBinContent(i+1));
            hp->SetBinContent(i+1, mp->GetBinContent(i+1));
          }
          else if(i>=nbins_CC && i<2*nbins_CC){
            hn->SetBinContent(i+1, en->GetBinContent(i-nbins_CC+1));
            ho->SetBinContent(i+1, eo->GetBinContent(i-nbins_CC+1));
            hp->SetBinContent(i+1, ep->GetBinContent(i-nbins_CC+1));
          }
          else{
            hn->SetBinContent(i+1, nue_o->GetBinContent(i-2*nbins_CC+1));
            ho->SetBinContent(i+1, nue_o->GetBinContent(i-2*nbins_CC+1));
            hp->SetBinContent(i+1, nue_o->GetBinContent(i-2*nbins_CC+1));
          }
        }

        covariance(ho, hn, hcv_n[j][iw], frhcv_n[j][iw], c, Form("%s_%d_n",name[iw],j));
        covariance(ho, hp, hcv_p[j][iw], frhcv_p[j][iw], c, Form("%s_%d_p",name[iw],j));
        covariance_pn(ho, hp, hn, hcv[j][iw], frhcv[j][iw], c, Form("%s_%d",name[iw],j));
      }

      else{
        for(int i=0; i<nbins; i++){
          if(i<nbins_CC){
            hn->SetBinContent(i+1, mo->GetBinContent(i+1));
            ho->SetBinContent(i+1, mo->GetBinContent(i+1));
            hp->SetBinContent(i+1, mo->GetBinContent(i+1));
          }
          else if(i>=nbins_CC && i<2*nbins_CC){
            hn->SetBinContent(i+1, eo->GetBinContent(i-nbins_CC+1));
            ho->SetBinContent(i+1, eo->GetBinContent(i-nbins_CC+1));
            hp->SetBinContent(i+1, eo->GetBinContent(i-nbins_CC+1));
          }
          else{
            hn->SetBinContent(i+1, nue_n->GetBinContent(i-2*nbins_CC+1));
            ho->SetBinContent(i+1, nue_o->GetBinContent(i-2*nbins_CC+1));
            hp->SetBinContent(i+1, nue_p->GetBinContent(i-2*nbins_CC+1));
          }
        }

        covariance(ho, hn, hcv_n[j][iw], frhcv_n[j][iw], c, Form("nue_%d_n",j));
        covariance(ho, hp, hcv_p[j][iw], frhcv_p[j][iw], c, Form("nue_%d_p",j));
        covariance_pn(ho, hp, hn, hcv[j][iw], frhcv[j][iw], c, Form("nue_%d",j));
      }
    }
  }

  c->Close();
/*
  TFile *out = new TFile(Form("%s/uncertainties/xS_covmtr/xS_covmtr_TEST.root",data_path), "RECREATE");
  for(int j=0; j<N_nucut; j++){
  for(int iw=0; iw<N_wgt+1; iw++){
    hcv_n[j][iw]->Write();
    hcv[j][iw]->Write();
    hcv_p[j][iw]->Write();

    frhcv_n[j][iw]->Write();
    frhcv[j][iw]->Write();
    frhcv_p[j][iw]->Write();
  }}

  out->Close();
*/
}
