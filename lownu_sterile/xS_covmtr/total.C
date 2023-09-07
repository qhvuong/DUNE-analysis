void total()
{
  TH2D *cv2 = new TH2D("cv2","",300,0,300,300,0,300);
  TH2D *cv  = new TH2D("cv","",300,0,300,300,0,300);
  TH2D *cv4 = new TH2D("cv4","",300,0,300,300,0,300);

  //std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};

  //const char *name[] = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};

  //std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};

  //const char *name[] = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};

  std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"};
  const char *name[] = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"};

  int N_wgt = namelist.size();
  std::cout << N_wgt << "\n";

  char name1[100] = "wgt_Mnv2p2hGaussEnhancement";
  char name2[100] = "wgt_MFP_N";


  for(int iw=0; iw<N_wgt; iw++){
    if(iw%10==0) std::cout << iw*100.0/N_wgt << "\n";

    TFile *f = new TFile(Form("new_mtr/%s_Covmtr32.root",name[iw]),"READ");

    TH2D *hcv2 = (TH2D*)f->Get("hcv2");    
    TH2D *hcv  = (TH2D*)f->Get("hcv");    
    //TH2D *hcv4 = (TH2D*)f->Get("hcv4");    

    //std::cout << hcv->GetNbinsX() << "\n";
    //std::cout << hcv->GetNbinsY() << "\n";

    //cv2->Add(hcv2);
    if(name[iw]==name1) cv->Add(hcv2);
    else cv->Add(hcv);
    //cv4->Add(hcv4);

/*
    for(int j=0; j<300; j++){
      for(int k=0; k<300; k++){
        cv2[j][k] = cv2[j][k] + hcv2[j][k];
        //cv[j+1][k+1]  = cv[j+1][k+1]  + hcv[j+1][k+1];
        //cv4[j+1][k+1] = cv4[j+1][k+1] + hcv4[j+1][k+1];
      }
    }
*/
  }

  //cv2->SetTitle("Cross Section Covariance Matrix (-1sigma)");
  cv->SetTitle("Cross Section Covariance Matrix");
  //cv4->SetTitle("Cross Section Covariance Matrix (+1sigma)");

  //cv2->SetStats(0);
  cv->SetStats(0);
  //cv4->SetStats(0);

  cv->SetMinimum(0);
  cv->SetMaximum(1);

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  TCanvas *c3 = new TCanvas("c3","",900,600);
  cv->Draw("colz");
  c3->SaveAs("sigma_5.png");

  TFile *out = new TFile("total_sigma_mtr_5.root","RECREATE");
  cv->Write();
  out->Close();

  std::cout << "test\n";
}
