static const int N = 1; // number of universes
static const int nbins = 250;
static const int nbins_CC = 100;


TMatrixD E_m(N, nbins_CC);
TMatrixD E_e(N, nbins_CC);
TMatrixD E_m_nue(N, nbins_CC);
TMatrixD E_e_nue(N, nbins_CC);
TMatrixD E_nue(N, nbins_CC);
  
TMatrixD ECovars_m2_m2   ( nbins_CC, nbins_CC );
TMatrixD ECovars_e2_e2   ( nbins_CC, nbins_CC );
TMatrixD ECovars_nue_nue ( nbins_CC, nbins_CC );

TMatrixD ECovars_m2_e2   ( nbins_CC, nbins_CC );
TMatrixD ECovars_m2_nue  ( nbins_CC, nbins_CC );
TMatrixD ECovars_e2_m2   ( nbins_CC, nbins_CC );
TMatrixD ECovars_e2_nue  ( nbins_CC, nbins_CC );
TMatrixD ECovars_nue_m2  ( nbins_CC, nbins_CC );
TMatrixD ECovars_nue_e2  ( nbins_CC, nbins_CC );


TMatrixD frECovars_m2_m2   ( nbins_CC, nbins_CC );
TMatrixD frECovars_e2_e2   ( nbins_CC, nbins_CC );
TMatrixD frECovars_nue_nue ( nbins_CC, nbins_CC );

TMatrixD frECovars_m2_e2   ( nbins_CC, nbins_CC );
TMatrixD frECovars_m2_nue  ( nbins_CC, nbins_CC );
TMatrixD frECovars_e2_m2   ( nbins_CC, nbins_CC );
TMatrixD frECovars_e2_nue  ( nbins_CC, nbins_CC );
TMatrixD frECovars_nue_m2  ( nbins_CC, nbins_CC );
TMatrixD frECovars_nue_e2  ( nbins_CC, nbins_CC );

TMatrixD ECovars2  ( nbins, nbins );
TMatrixD frECovars2  ( nbins, nbins );



TMatrixD ECovars_m4_m4   ( nbins_CC, nbins_CC );
TMatrixD ECovars_e4_e4   ( nbins_CC, nbins_CC );

TMatrixD ECovars_m4_e4   ( nbins_CC, nbins_CC );
TMatrixD ECovars_m4_nue  ( nbins_CC, nbins_CC );
TMatrixD ECovars_e4_m4   ( nbins_CC, nbins_CC );
TMatrixD ECovars_e4_nue  ( nbins_CC, nbins_CC );
TMatrixD ECovars_nue_m4  ( nbins_CC, nbins_CC );
TMatrixD ECovars_nue_e4  ( nbins_CC, nbins_CC );

TMatrixD frECovars_m4_m4   ( nbins_CC, nbins_CC );
TMatrixD frECovars_e4_e4   ( nbins_CC, nbins_CC );

TMatrixD frECovars_m4_e4   ( nbins_CC, nbins_CC );
TMatrixD frECovars_m4_nue  ( nbins_CC, nbins_CC );
TMatrixD frECovars_e4_m4   ( nbins_CC, nbins_CC );
TMatrixD frECovars_e4_nue  ( nbins_CC, nbins_CC );
TMatrixD frECovars_nue_m4  ( nbins_CC, nbins_CC );
TMatrixD frECovars_nue_e4  ( nbins_CC, nbins_CC );

TMatrixD ECovars4  ( nbins, nbins );
TMatrixD frECovars4  ( nbins, nbins );



TMatrixD ECovars_mm   ( nbins_CC, nbins_CC );
TMatrixD ECovars_ee   ( nbins_CC, nbins_CC );

TMatrixD ECovars_me   ( nbins_CC, nbins_CC );
TMatrixD ECovars_mnue  ( nbins_CC, nbins_CC );
TMatrixD ECovars_em   ( nbins_CC, nbins_CC );
TMatrixD ECovars_enue  ( nbins_CC, nbins_CC );
TMatrixD ECovars_nuem  ( nbins_CC, nbins_CC );
TMatrixD ECovars_nuee  ( nbins_CC, nbins_CC );

TMatrixD frECovars_mm   ( nbins_CC, nbins_CC );
TMatrixD frECovars_ee   ( nbins_CC, nbins_CC );

TMatrixD frECovars_me   ( nbins_CC, nbins_CC );
TMatrixD frECovars_mnue  ( nbins_CC, nbins_CC );
TMatrixD frECovars_em   ( nbins_CC, nbins_CC );
TMatrixD frECovars_enue  ( nbins_CC, nbins_CC );
TMatrixD frECovars_nuem  ( nbins_CC, nbins_CC );
TMatrixD frECovars_nuee  ( nbins_CC, nbins_CC );

TMatrixD ECovars  ( nbins, nbins );
TMatrixD frECovars  ( nbins, nbins );


void test()
{

  int cutNu=3;

  TFile *f     = new TFile("/dune/app/users/qvuong/data/lownu/CC_output.root");
  TFile *f_nue = new TFile("/dune/app/users/qvuong/data/lownu/nue_output.root");

  //std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};

  //const char *name[] = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"}; 

  std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"};
  const char *name[] = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"}; 
  //std::list <const char *> namelist = {"wgt_BeRPA_B", "wgt_BeRPA_D"};
  //const char *name[] = {"wgt_BeRPA_B", "wgt_BeRPA_D"}; 

  int N_wgt = namelist.size();

  std::cout << N_wgt << "\n";

  for(int iw=0; iw<N_wgt; iw++){

  TH1D *m1       = (TH1D*)f->Get(Form("%s_m_hElep%d_sigma1",name[iw],cutNu));
  TH1D *e1   	 = (TH1D*)f->Get(Form("%s_e_hElep%d_sigma1",name[iw],cutNu));
  TH1D *m2       = (TH1D*)f->Get(Form("%s_m_hElep%d_sigma2",name[iw],cutNu));
  TH1D *e2   	 = (TH1D*)f->Get(Form("%s_e_hElep%d_sigma2",name[iw],cutNu));
  TH1D *CC_m_nom = (TH1D*)f->Get(Form("%s_m_hElep%d_sigma3",name[iw],cutNu));
  TH1D *CC_e_nom = (TH1D*)f->Get(Form("%s_e_hElep%d_sigma3",name[iw],cutNu));
  TH1D *m4       = (TH1D*)f->Get(Form("%s_m_hElep%d_sigma4",name[iw],cutNu));
  TH1D *e4   	 = (TH1D*)f->Get(Form("%s_e_hElep%d_sigma4",name[iw],cutNu));
  TH1D *m5       = (TH1D*)f->Get(Form("%s_m_hElep%d_sigma5",name[iw],cutNu));
  TH1D *e5   	 = (TH1D*)f->Get(Form("%s_e_hElep%d_sigma5",name[iw],cutNu));
  TH1D *nue 	 = (TH1D*)f_nue->Get("hElep2");
  TH1D *nue_nom  = (TH1D*)f_nue->Get("hElep2");


  for( int i = 0; i < nbins_CC; ++i ) { // columns
    for( int j = 0; j < nbins_CC; ++j ) { // columns
      // compute column covariance
      double covar_m2_m2 = 0.;
      double covar_e2_e2 = 0.;
      double covar_nue_nue = 0.;
      double covar_m2_e2 = 0.;
      double covar_m2_nue = 0.;
      double covar_e2_m2 = 0.;
      double covar_e2_nue = 0.;
      double covar_nue_m2 = 0.;
      double covar_nue_e2 = 0.;

      double var_m2_i = 0.;
      double var_m2_j = 0.;
      double var_e2_i = 0.;
      double var_e2_j = 0.;
      double var_nue_i = 0.;
      double var_nue_j = 0.;


      double covar_m4_m4 = 0.;
      double covar_e4_e4 = 0.;
      double covar_m4_e4 = 0.;
      double covar_m4_nue = 0.;
      double covar_e4_m4 = 0.;
      double covar_e4_nue = 0.;
      double covar_nue_m4 = 0.;
      double covar_nue_e4 = 0.;

      double var_m4_i = 0.;
      double var_m4_j = 0.;
      double var_e4_i = 0.;
      double var_e4_j = 0.;


      double covar_mm = 0.;
      double covar_ee = 0.;
      double covar_me = 0.;
      double covar_mnue = 0.;
      double covar_em = 0.;
      double covar_enue = 0.;
      double covar_nuem = 0.;
      double covar_nuee = 0.;

      //for( int k = 0; k < N; ++k ) { // rows
      int k = 0; {
        covar_m2_m2   += (m2->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (m2->GetBinContent(j+1)  - CC_m_nom->GetBinContent(j+1));
        covar_e2_e2   += (e2->GetBinContent(i+1)  - CC_e_nom->GetBinContent(i+1)) * (e2->GetBinContent(j+1)  - CC_e_nom->GetBinContent(j+1)); 
        //covar_nue_nue += (nue->GetBinContent(i+1) - nue_nom->GetBinContent(i+1))  * (nue->GetBinContent(j+1) - nue_nom->GetBinContent(j+1)); 

        covar_m2_e2   += (m2->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (e2->GetBinContent(j+1)  - CC_e_nom->GetBinContent(j+1));
        //covar_m2_nue  += (m2->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (nue->GetBinContent(j+1) - nue_nom->GetBinContent(j+1));
        covar_e2_m2   += (e2->GetBinContent(i+1)  - CC_e_nom->GetBinContent(i+1)) * (m2->GetBinContent(j+1)  - CC_m_nom->GetBinContent(j+1)); 
        //covar_e2_nue  += (e2->GetBinContent(i+1)  - CC_e_nom->GetBinContent(i+1)) * (nue->GetBinContent(j+1) - nue_nom->GetBinContent(j+1)); 
        //covar_nue_m2  += (nue->GetBinContent(i+1) - nue_nom->GetBinContent(i+1))  * (m2->GetBinContent(j+1)  - CC_m_nom->GetBinContent(j+1)); 
        //covar_nue_e2  += (nue->GetBinContent(i+1) - nue_nom->GetBinContent(i+1))  * (e2->GetBinContent(j+1)  - CC_e_nom->GetBinContent(j+1)); 


        covar_m4_m4   += (m4->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (m4->GetBinContent(j+1)  - CC_m_nom->GetBinContent(j+1));
        covar_e4_e4   += (e4->GetBinContent(i+1)  - CC_e_nom->GetBinContent(i+1)) * (e4->GetBinContent(j+1)  - CC_e_nom->GetBinContent(j+1)); 

        covar_m4_e4   += (m4->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (e4->GetBinContent(j+1)  - CC_e_nom->GetBinContent(j+1));
        //covar_m4_nue  += (m4->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (nue->GetBinContent(j+1) - nue_nom->GetBinContent(j+1));
        covar_e4_m4   += (e4->GetBinContent(i+1)  - CC_e_nom->GetBinContent(i+1)) * (m4->GetBinContent(j+1)  - CC_m_nom->GetBinContent(j+1)); 
        //covar_e4_nue  += (e4->GetBinContent(i+1)  - CC_e_nom->GetBinContent(i+1)) * (nue->GetBinContent(j+1) - nue_nom->GetBinContent(j+1)); 
        //covar_nue_m4  += (nue->GetBinContent(i+1) - nue_nom->GetBinContent(i+1))  * (m4->GetBinContent(j+1)  - CC_m_nom->GetBinContent(j+1)); 
        //covar_nue_e4  += (nue->GetBinContent(i+1) - nue_nom->GetBinContent(i+1))  * (e4->GetBinContent(j+1)  - CC_e_nom->GetBinContent(j+1)); 


        covar_mm = (abs(covar_m2_m2) + abs(covar_m4_m4))/2;
        covar_ee = (abs(covar_e2_e2) + abs(covar_e4_e4))/2;

        covar_me = (abs(covar_m2_e2) + abs(covar_m4_e4))/2;
        //covar_mnue = (abs(covar_m2_nue) + abs(covar_m4_nue))/2;
        covar_em = (abs(covar_e2_m2) + abs(covar_e4_m4))/2;
        //covar_enue = (abs(covar_e2_nue) + abs(covar_e4_nue))/2;
        //covar_nuem = (abs(covar_nue_m2) + abs(covar_nue_m4))/2;
        //covar_nuee = (abs(covar_nue_e2) + abs(covar_nue_e4))/2;
      }
      ECovars_m2_m2[i][j] = covar_m2_m2;
      ECovars_e2_e2[i][j] = covar_e2_e2;
      ECovars_m2_e2[i][j] = covar_m2_e2;
      ECovars_e2_m2[i][j] = covar_e2_m2;
      
      ECovars_m4_m4[i][j] = covar_m4_m4;
      ECovars_e4_e4[i][j] = covar_e4_e4;
      ECovars_m4_e4[i][j] = covar_m4_e4;
      ECovars_e4_m4[i][j] = covar_e4_m4;


      ECovars_mm[i][j] = covar_mm;
      ECovars_ee[i][j] = covar_ee;
      ECovars_me[i][j] = covar_me;
      //ECovars_mnue[i][j]  = covar_mnue;
      ECovars_em[i][j] = covar_em;
      //ECovars_enue[i][j]  = covar_enue;
      //ECovars_nuem[i][j]  = covar_nuem;
      //ECovars_nuee[i][j]  = covar_nuee;

 
      frECovars_mm[i][j] = covar_mm / (N * CC_m_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      frECovars_ee[i][j] = covar_ee / (N * CC_e_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      frECovars_me[i][j]   = covar_me   /(N * CC_m_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      //frECovars_mnue[i][j]  = covar_mnue  /(N * CC_m_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      frECovars_em[i][j]   = covar_em   /(N * CC_e_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      //frECovars_enue[i][j]  = covar_enue  /(N * CC_e_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      //frECovars_nuem[i][j]  = covar_nuem  /(N * nue_nom->GetBinContent(i+1)  * CC_m_nom->GetBinContent(j+1));
      //frECovars_nuee[i][j]  = covar_nuee  /(N * nue_nom->GetBinContent(i+1)  * CC_e_nom->GetBinContent(j+1));
      
      frECovars_m2_m2[i][j] = covar_m2_m2 / (N * CC_m_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      frECovars_e2_e2[i][j] = covar_e2_e2 / (N * CC_e_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      frECovars_m2_e2[i][j] = covar_m2_e2 /(N * CC_m_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      frECovars_e2_m2[i][j] = covar_e2_m2 /(N * CC_e_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      
      frECovars_m4_m4[i][j] = covar_m4_m4 / (N * CC_m_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      frECovars_e4_e4[i][j] = covar_e4_e4 / (N * CC_e_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      frECovars_m4_e4[i][j] = covar_m4_e4 /(N * CC_m_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      frECovars_e4_m4[i][j] = covar_e4_m4 /(N * CC_e_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
    }
  }

  for(int i = 0; i < nbins; i++) {
    for(int j = 0; j < nbins; j++) {

      if(i<nbins_CC){
        if(j<nbins_CC)                  ECovars[i][j] = ECovars_mm[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars[i][j] = ECovars_me[i][j-nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  ECovars[i][j] = ECovars_em[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars[i][j] = ECovars_ee[i-nbins_CC][j-nbins_CC];}

      if(i>=2*nbins_CC || j>=2*nbins_CC) ECovars[i][j] = 0.;


      if(i<nbins_CC){
        if(j<nbins_CC)                  frECovars[i][j] = frECovars_mm[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars[i][j] = frECovars_me[i][j-nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  frECovars[i][j] = frECovars_em[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars[i][j] = frECovars_ee[i-nbins_CC][j-nbins_CC];}

      if(i>=2*nbins_CC || j>=2*nbins_CC) frECovars[i][j] = 0.;


      if(i<nbins_CC){
        if(j<nbins_CC)                  ECovars2[i][j] = ECovars_m2_m2[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars2[i][j] = ECovars_m2_e2[i][j-nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  ECovars2[i][j] = ECovars_e2_m2[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars2[i][j] = ECovars_e2_e2[i-nbins_CC][j-nbins_CC];}

      if(i>=2*nbins_CC || j>=2*nbins_CC) ECovars2[i][j] = 0.;


      if(i<nbins_CC){
        if(j<nbins_CC)                  frECovars2[i][j] = frECovars_m2_m2[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars2[i][j] = frECovars_m2_e2[i][j-nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  frECovars2[i][j] = frECovars_e2_m2[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars2[i][j] = frECovars_e2_e2[i-nbins_CC][j-nbins_CC];}

      if(i>=2*nbins_CC || j>=2*nbins_CC) frECovars2[i][j] = 0.;
      

      if(i<nbins_CC){
        if(j<nbins_CC)                  ECovars4[i][j] = ECovars_m4_m4[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars4[i][j] = ECovars_m4_e4[i][j-nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  ECovars4[i][j] = ECovars_e4_m4[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars4[i][j] = ECovars_e4_e4[i-nbins_CC][j-nbins_CC];}

      if(i>=2*nbins_CC || j>=2*nbins_CC) ECovars4[i][j] = 0.;


      if(i<nbins_CC){
        if(j<nbins_CC)                  frECovars4[i][j] = frECovars_m4_m4[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars4[i][j] = frECovars_m4_e4[i][j-nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  frECovars4[i][j] = frECovars_e4_m4[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars4[i][j] = frECovars_e4_e4[i-nbins_CC][j-nbins_CC];}

      if(i>=2*nbins_CC || j>=2*nbins_CC) frECovars4[i][j] = 0.;
    }
  }

  TH2D *hcv2 = new TH2D("hcv2","",nbins,0,nbins,nbins,0,nbins);
  TH2D *hcv  = new TH2D("hcv","",nbins,0,nbins,nbins,0,nbins);
  TH2D *hcv4 = new TH2D("hcv4","",nbins,0,nbins,nbins,0,nbins);
  TH2D *frhcv2 = new TH2D("frhcv2","",nbins,0,nbins,nbins,0,nbins);
  TH2D *frhcv  = new TH2D("frhcv","",nbins,0,nbins,nbins,0,nbins);
  TH2D *frhcv4 = new TH2D("frhcv4","",nbins,0,nbins,nbins,0,nbins);

  hcv2->SetTitle(Form("%s -1sigma",name[iw]));
  hcv->SetTitle(Form("%s 0sigma",name[iw]));
  hcv4->SetTitle(Form("%s +1sigma",name[iw]));
  frhcv2->SetTitle(Form("fractional %s -1sigma",name[iw]));
  frhcv->SetTitle(Form("fractional %s 0sigma",name[iw]));
  frhcv4->SetTitle(Form("fractional %s +1sigma",name[iw]));

  for(int i=0; i<nbins; i++) {
    for(int j=0; j<nbins; j++) {
      hcv2->SetBinContent(i+1, j+1, ECovars2[i][j]);
      hcv->SetBinContent(i+1, j+1, ECovars[i][j]);
      hcv4->SetBinContent(i+1, j+1, ECovars4[i][j]);
      frhcv2->SetBinContent(i+1, j+1, frECovars2[i][j]);
      frhcv->SetBinContent(i+1, j+1, frECovars[i][j]);
      frhcv4->SetBinContent(i+1, j+1, frECovars4[i][j]);
    }
  }


  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
  hcv2->SetStats(0);
  hcv->SetStats(0);
  hcv4->SetStats(0);
  frhcv2->SetStats(0);
  frhcv->SetStats(0);
  frhcv4->SetStats(0);
  TCanvas *c = new TCanvas("c","",1200,600);
  c->Divide(3,2);
  c->cd(1);
  hcv2->Draw("colz");
  c->cd(2);
  hcv->Draw("colz");
  c->cd(3);
  hcv4->Draw("colz");
  c->cd(4);
  frhcv2->Draw("colz");
  c->cd(5);
  frhcv->Draw("colz");
  c->cd(6);
  frhcv4->Draw("colz");
  c->SaveAs(Form("%s_sigmtr%d.png",name[iw],cutNu));
  
  TCanvas *c1 = new TCanvas("c1","",1200,600);
  c1->Divide(3,2);
  c1->cd(1);
  gPad->SetLogz();
  hcv2->Draw("colz");
  c1->cd(2);
  gPad->SetLogz();
  hcv->Draw("colz");
  c1->cd(3);
  gPad->SetLogz();
  hcv4->Draw("colz");
  c1->cd(4);
  frhcv2->Draw("colz");
  c1->cd(5);
  frhcv->Draw("colz");
  c1->cd(6);
  frhcv4->Draw("colz");
  c1->SaveAs(Form("%s_sigmtr%d_logz.png",name[iw],cutNu));


  TFile *out = new TFile(Form("%s_sigmtr%d.root",name[iw],cutNu), "RECREATE");
  hcv2->Write();
  hcv->Write();
  hcv4->Write();
  frhcv2->Write();
  frhcv->Write();
  frhcv4->Write();
  out->Close();

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  }
}

