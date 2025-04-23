#include "DUNEStyle.h"

static const int nbins_CC = 56;
static const int nbins_nue = 8;
static const int nbins = 2*nbins_CC + nbins_nue;


void plotCovfr(TH2D* cov, const char *name, TCanvas *c)
{ 
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

    c->Clear();
    cov->SetMinimum(1E-5);
    cov->SetMaximum(1.);
    cov->SetStats(0);
    cov->Draw("colz");
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
    c->RedrawAxis();
    c->SaveAs(Form("%s.png",name));
}


void plotCov(TH2D* cov, const char *name, TCanvas *c){
  //int nbins = hnom->GetNbinsX();
  //std::cout << "nbins = " << nbins << "\n";

  c->Clear();
  
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

    //TCanvas *c = new TCanvas("c", "", 800, 600);
    cov->SetMinimum(1E-3);
    cov->SetMaximum(1E11);
    cov->SetStats(0);
    //cov->Draw("colz");
    //c->SaveAs(Form("%s.png",name));
    
    cov->Draw("colz");
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
    c->RedrawAxis();
    c->SaveAs(Form("%s.png",name));
  
}





void total()
{
  TH2D *hcv  = new TH2D("hcv","Total Absolute Covariance Matrix (#nu < 0.3 GeV)",nbins,0,nbins,nbins,0,nbins);
  TH2D *frhcv  = new TH2D("frhcv","Total Fractional Covariance Matrix (#nu < 0.3 GeV)",nbins,0,nbins,nbins,0,nbins);

  //std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};
  //const char *name[] = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};

  std::list <const char *> namelist_avg = {"wgt_MaCCQE", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio"};
  std::list <const char *> namelist_pos = {"wgt_VecFFCCQEshape", "wgt_Theta_Delta2Npi", "wgt_CCQEPauliSupViaKF", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};
  std::list <const char *> namelist_neg = {"wgt_MFP_N", "wgt_Mnv2p2hGaussEnhancement"};

  const char *name_avg[] = {"wgt_MaCCQE", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio"};
  const char *name_pos[] = {"wgt_VecFFCCQEshape", "wgt_Theta_Delta2Npi", "wgt_CCQEPauliSupViaKF", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};
  const char *name_neg[] = {"wgt_MFP_N", "wgt_Mnv2p2hGaussEnhancement"};

  int N_wgt_avg = namelist_avg.size();
  int N_wgt_pos = namelist_pos.size();
  int N_wgt_neg = namelist_neg.size();
  std::cout << N_wgt_avg << "\t" << N_wgt_pos << "\t" << N_wgt_neg << "\n";
  int N_wgt = 68;

  TFile *f = new TFile("/exp/dune/app/users/qvuong/data/lownu/uncertainties/xS_covmtr/xS_covmtr_TEST.root","READ");

  for(int iw=0; iw<N_wgt_avg; iw++){
    std::cout << iw << "\n";
  
    TH2D *cv  = (TH2D*)f->Get(Form("%s_hcv4",name_avg[iw]));     
    TH2D *frcv  = (TH2D*)f->Get(Form("%s_frhcv4",name_avg[iw]));  
    hcv->Add(cv);  
    frhcv->Add(frcv);
  }

  for(int iw=0; iw<N_wgt_pos; iw++){
    std::cout << iw << "\n";

    TH2D *cv = (TH2D*)f->Get(Form("%s_hcv4_p",name_pos[iw]));   
    TH2D *frcv = (TH2D*)f->Get(Form("%s_frhcv4_p",name_pos[iw]));    
    hcv->Add(cv); 
    frhcv->Add(frcv);
  }

  for(int iw=0; iw<N_wgt_neg; iw++){
    std::cout << iw << "\n";

    TH2D *cv  = (TH2D*)f->Get(Form("%s_hcv4_n",name_neg[iw]));    
    TH2D *frcv = (TH2D*)f->Get(Form("%s_frhcv4_n",name_neg[iw]));    
    hcv->Add(cv);    
    frhcv->Add(frcv);
  }



  TColor::InvertPalette();

  TCanvas *c = new TCanvas("c","",800,600);
  c->SetLogz();
  plotCov(hcv, "xS", c);
  plotCovfr(frhcv, "frxS", c);

  TH2D *nue = (TH2D*)f->Get("nue_hcv4");
  TH2D *frnue = (TH2D*)f->Get("nue_frhcv4");
  hcv->Add(nue);
  frhcv->Add(frnue);

  plotCov(hcv, "xS_nue", c);
  plotCovfr(frhcv, "frxS_nue", c);


/*
  TFile *out = new TFile("xS_unc.root", "RECREATE");
  hcv->Write();
  frhcv->Write();
  out->Close();
*/


}
