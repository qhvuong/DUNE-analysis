static const int N = 1; // number of universes
static const int nbins_CC = 58;
static const int nbins_nue = 8;
static const int nbins = 2*nbins_CC + nbins_nue;


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



void ElepCov()
{

  TFile *f     = new TFile("/exp/dune/app/users/qvuong/data/lownu/CC_output_58.root");
  TFile *f_nue = new TFile("/exp/dune/app/users/qvuong/data/lownu/nue_output_test.root");
  int cutNu = 3;

  //std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};

  //const char *name[] = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"}; 

  std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"};
  const char *name[] = {"wgt_MaCCQE", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_Mnv2p2hGaussEnhancement"}; 

  //std::list <const char *> namelist = {"wgt_Mnv2p2hGaussEnhancement"};
  //const char *name[] = {"wgt_Mnv2p2hGaussEnhancement"}; 

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
  TH1D *nue 	 = (TH1D*)f_nue->Get("hElep0");
  TH1D *nue_nom  = (TH1D*)f_nue->Get("hElep0");

  TH1D *h1 = (TH1D*)e1->Clone();
  TH1D *h2 = (TH1D*)e2->Clone();
  TH1D *h4 = (TH1D*)e4->Clone();
  TH1D *h5 = (TH1D*)e5->Clone();
  

  for( int i = 0; i < nbins; ++i ) { // columns
    for( int j = 0; j < nbins; ++j ) { // columns
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
      

      if(i<nbins_CC) {
        if(j<nbins_CC)                  covar_m2_m2   += (m2->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (m2->GetBinContent(j+1)             - CC_m_nom->GetBinContent(j+1));
        if(j>=nbins_CC && j<2*nbins_CC) covar_m2_e2   += (m2->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (e2->GetBinContent(j-nbins_CC+1)    - CC_e_nom->GetBinContent(j-nbins_CC+1));
        if(j>=2*nbins_CC)               covar_m2_nue  += (m2->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (nue->GetBinContent(j-2*nbins_CC+1) - nue_nom->GetBinContent(j-2*nbins_CC+1));}

      if(i>=nbins_CC && i<2*nbins_CC) {
        if(j<nbins_CC)                  covar_e2_m2   += (e2->GetBinContent(i-nbins_CC+1)  - CC_e_nom->GetBinContent(i-nbins_CC+1)) * (m2->GetBinContent(j+1)             - CC_m_nom->GetBinContent(j+1)); 
        if(j>=nbins_CC && j<2*nbins_CC) covar_e2_e2   += (e2->GetBinContent(i-nbins_CC+1)  - CC_e_nom->GetBinContent(i-nbins_CC+1)) * (e2->GetBinContent(j-nbins_CC+1)    - CC_e_nom->GetBinContent(j-nbins_CC+1)); 
        if(j>=2*nbins_CC)               covar_e2_nue  += (e2->GetBinContent(i-nbins_CC+1)  - CC_e_nom->GetBinContent(i-nbins_CC+1)) * (nue->GetBinContent(j-2*nbins_CC+1) - nue_nom->GetBinContent(j-2*nbins_CC+1));}
    
      if(i>=2*nbins_CC) {
        if(j<nbins_CC)                  covar_nue_m2  += (nue->GetBinContent(i-2*nbins_CC+1) - nue_nom->GetBinContent(i-2*nbins_CC+1))  * (m2->GetBinContent(j+1)             - CC_m_nom->GetBinContent(j+1)); 
        if(j>=nbins_CC && j<2*nbins_CC) covar_nue_e2  += (nue->GetBinContent(i-2*nbins_CC+1) - nue_nom->GetBinContent(i-2*nbins_CC+1))  * (e2->GetBinContent(j-nbins_CC+1)    - CC_e_nom->GetBinContent(j-nbins_CC+1)); 
        if(j>=2*nbins_CC)               covar_nue_nue += (nue->GetBinContent(i-2*nbins_CC+1) - nue_nom->GetBinContent(i-2*nbins_CC+1))  * (nue->GetBinContent(j-2*nbins_CC+1) - nue_nom->GetBinContent(j-2*nbins_CC+1));} 

      

      if(i<nbins_CC)                  var_m2_i  += (m2->GetBinContent(i+1)             - CC_m_nom->GetBinContent(i+1))           * (m2->GetBinContent(i+1)             - CC_m_nom->GetBinContent(i+1));
      if(i>=nbins_CC && i<2*nbins_CC) var_e2_i  += (e2->GetBinContent(i-nbins_CC+1)    - CC_e_nom->GetBinContent(i-nbins_CC+1))  * (e2->GetBinContent(i-nbins_CC+1)    - CC_e_nom->GetBinContent(i-nbins_CC+1));
      if(i>=2*nbins_CC)               var_nue_i += (nue->GetBinContent(i-2*nbins_CC+1) - nue_nom->GetBinContent(i-2*nbins_CC+1)) * (nue->GetBinContent(i-2*nbins_CC+1) - nue_nom->GetBinContent(i-2*nbins_CC+1));
        
      if(j<nbins_CC)                  var_m2_j  += (m2->GetBinContent(j+1)             - CC_m_nom->GetBinContent(j+1))           * (m2->GetBinContent(j+1)             - CC_m_nom->GetBinContent(j+1));
      if(j>=nbins_CC && j<2*nbins_CC) var_e2_j  += (e2->GetBinContent(j-nbins_CC+1)    - CC_e_nom->GetBinContent(j-nbins_CC+1))  * (e2->GetBinContent(j-nbins_CC+1)    - CC_e_nom->GetBinContent(j-nbins_CC+1));
      if(j>=2*nbins_CC)               var_nue_j += (nue->GetBinContent(j-2*nbins_CC+1) - nue_nom->GetBinContent(j-2*nbins_CC+1)) * (nue->GetBinContent(j-2*nbins_CC+1) - nue_nom->GetBinContent(j-2*nbins_CC+1));



      if(i<nbins_CC) {
        if(j<nbins_CC)                  covar_m4_m4   += (m4->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (m4->GetBinContent(j+1)             - CC_m_nom->GetBinContent(j+1));
        if(j>=nbins_CC && j<2*nbins_CC) covar_m4_e4   += (m4->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (e4->GetBinContent(j-nbins_CC+1)    - CC_e_nom->GetBinContent(j-nbins_CC+1));
        if(j>=2*nbins_CC)               covar_m4_nue  += (m4->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (nue->GetBinContent(j-2*nbins_CC+1) - nue_nom->GetBinContent(j-2*nbins_CC+1));}

      if(i>=nbins_CC && i<2*nbins_CC) {
        if(j<nbins_CC)                  covar_e4_m4   += (e4->GetBinContent(i-nbins_CC+1)  - CC_e_nom->GetBinContent(i-nbins_CC+1)) * (m4->GetBinContent(j+1)             - CC_m_nom->GetBinContent(j+1)); 
        if(j>=nbins_CC && j<2*nbins_CC) covar_e4_e4   += (e4->GetBinContent(i-nbins_CC+1)  - CC_e_nom->GetBinContent(i-nbins_CC+1)) * (e4->GetBinContent(j-nbins_CC+1)    - CC_e_nom->GetBinContent(j-nbins_CC+1)); 
        if(j>=2*nbins_CC)               covar_e4_nue  += (e4->GetBinContent(i-nbins_CC+1)  - CC_e_nom->GetBinContent(i-nbins_CC+1)) * (nue->GetBinContent(j-2*nbins_CC+1) - nue_nom->GetBinContent(j-2*nbins_CC+1));} 

      if(i>=2*nbins_CC) {
        if(j<nbins_CC)                  covar_nue_m4  += (nue->GetBinContent(i-2*nbins_CC+1) - nue_nom->GetBinContent(i-2*nbins_CC+1))  * (m4->GetBinContent(j+1)           - CC_m_nom->GetBinContent(j+1)); 
        if(j>=nbins_CC && j<2*nbins_CC) covar_nue_e4  += (nue->GetBinContent(i-2*nbins_CC+1) - nue_nom->GetBinContent(i-2*nbins_CC+1))  * (e4->GetBinContent(j-nbins_CC+1)  - CC_e_nom->GetBinContent(j-nbins_CC+1));} 



      if(i<nbins_CC)                  var_m4_i += (m4->GetBinContent(i+1)          - CC_m_nom->GetBinContent(i+1))          * (m4->GetBinContent(i+1)          - CC_m_nom->GetBinContent(i+1));
      if(i>=nbins_CC && i<2*nbins_CC) var_e4_i += (e4->GetBinContent(i-nbins_CC+1) - CC_e_nom->GetBinContent(i-nbins_CC+1)) * (e4->GetBinContent(i-nbins_CC+1) - CC_e_nom->GetBinContent(i-nbins_CC+1));
        
      if(j<nbins_CC)                  var_m4_j += (m4->GetBinContent(j+1)          - CC_m_nom->GetBinContent(j+1))          * (m4->GetBinContent(j+1)          - CC_m_nom->GetBinContent(j+1));
      if(j>=nbins_CC && j<2*nbins_CC) var_e4_j += (e4->GetBinContent(j-nbins_CC+1) - CC_e_nom->GetBinContent(j-nbins_CC+1)) * (e4->GetBinContent(j-nbins_CC+1) - CC_e_nom->GetBinContent(j-nbins_CC+1));

      

      if(i<nbins_CC) {
        if(j<nbins_CC)                  covar_mm = (abs(covar_m2_m2) + abs(covar_m4_m4))/2;
        if(j>=nbins_CC && j<2*nbins_CC) covar_me = (abs(covar_m2_e2) + abs(covar_m4_e4))/2;
        if(j>=2*nbins_CC)               covar_mnue = (abs(covar_m2_nue) + abs(covar_m4_nue))/2;}

      if(i>=nbins_CC && i<2*nbins_CC) {
        if(j<nbins_CC)                  covar_em = (abs(covar_e2_m2) + abs(covar_e4_m4))/2;
        if(j>=nbins_CC && j<2*nbins_CC) covar_ee = (abs(covar_e2_e2) + abs(covar_e4_e4))/2;
        if(j>=2*nbins_CC)               covar_enue = (abs(covar_e2_nue) + abs(covar_e4_nue))/2;}

      if(i>=2*nbins_CC) {
        if(j<nbins_CC)                  covar_nuem = (abs(covar_nue_m2) + abs(covar_nue_m4))/2;
        if(j>=nbins_CC && j<2*nbins_CC) covar_nuee = (abs(covar_nue_e2) + abs(covar_nue_e4))/2;}



      // TOTAL COVARIANCE MATRICES (SIGMA -1)
      if(i<nbins_CC) {                  
        if(j<nbins_CC)                  ECovars_m2_m2[i][j]              = covar_m2_m2;
        if(j>=nbins_CC && j<2*nbins_CC) ECovars_m2_e2[i][j-nbins_CC]     = covar_m2_e2;
        if(j>=2*nbins_CC)               ECovars_m2_nue[i][j-2*nbins_CC]  = covar_m2_nue;}

      if(i>=nbins_CC && i<2*nbins_CC) {
        if(j<nbins_CC)                  ECovars_e2_m2[i-nbins_CC][j]              = covar_e2_m2;
	if(j>=nbins_CC && j<2*nbins_CC) ECovars_e2_e2[i-nbins_CC][j-nbins_CC]     = covar_e2_e2;
	if(j>=2*nbins_CC)               ECovars_e2_nue[i-nbins_CC][j-2*nbins_CC]  = covar_e2_nue;}

      if(i>=2*nbins_CC) {
        if(j<nbins_CC)                  ECovars_nue_m2[i-2*nbins_CC][j]             = covar_nue_m2;
	if(j>=nbins_CC && j<2*nbins_CC) ECovars_nue_e2[i-2*nbins_CC][j-nbins_CC]    = covar_nue_e2;
	if(j>=2*nbins_CC)               ECovars_nue_nue[i-2*nbins_CC][j-2*nbins_CC] = covar_nue_nue;}
      


      if(i<nbins_CC) {
        if(j<nbins_CC)                  frECovars_m2_m2[i][j]             = covar_m2_m2  /(N * CC_m_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
	if(j>=nbins_CC && j<2*nbins_CC) frECovars_m2_e2[i][j-nbins_CC]    = covar_m2_e2  /(N * CC_m_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j-nbins_CC+1));
	if(j>=2*nbins_CC)               frECovars_m2_nue[i][j-2*nbins_CC] = covar_m2_nue /(N * CC_m_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}
	
      if(i>=nbins_CC && i<2*nbins_CC) {
        if(j<nbins_CC)                  frECovars_e2_m2[i-nbins_CC][j]             = covar_e2_m2  /(N * CC_e_nom->GetBinContent(i-nbins_CC+1) * CC_m_nom->GetBinContent(j+1));
	if(j>=nbins_CC && j<2*nbins_CC) frECovars_e2_e2[i-nbins_CC][j-nbins_CC]    = covar_e2_e2  /(N * CC_e_nom->GetBinContent(i-nbins_CC+1) * CC_e_nom->GetBinContent(j-nbins_CC+1));
	if(j>=2*nbins_CC)               frECovars_e2_nue[i-nbins_CC][j-2*nbins_CC] = covar_e2_nue /(N * CC_e_nom->GetBinContent(i-nbins_CC+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}

      if(i>=2*nbins_CC) {	
        if(j<nbins_CC)                  frECovars_nue_m2[i-2*nbins_CC][j]             = covar_nue_m2  /(N * nue_nom->GetBinContent(i-2*nbins_CC+1) * CC_m_nom->GetBinContent(j+1));
	if(j>=nbins_CC && j<2*nbins_CC) frECovars_nue_e2[i-2*nbins_CC][j-nbins_CC]    = covar_nue_e2  /(N * nue_nom->GetBinContent(i-2*nbins_CC+1) * CC_e_nom->GetBinContent(j-nbins_CC+1));
	if(j>=2*nbins_CC)               frECovars_nue_nue[i-2*nbins_CC][j-2*nbins_CC] = covar_nue_nue /(N * nue_nom->GetBinContent(i-2*nbins_CC+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}

      


      // TOTAL COVARIANCE MATRICES (SIGMA +1)
      if(i<nbins_CC) {                  
        if(j<nbins_CC)                  ECovars_m4_m4[i][j]              = covar_m4_m4;
        if(j>=nbins_CC && j<2*nbins_CC) ECovars_m4_e4[i][j-nbins_CC]     = covar_m4_e4;
        if(j>=2*nbins_CC)               ECovars_m4_nue[i][j-2*nbins_CC]  = covar_m4_nue;}

      if(i>=nbins_CC && i<2*nbins_CC) {
        if(j<nbins_CC)                  ECovars_e4_m4[i-nbins_CC][j]              = covar_e4_m4;
	if(j>=nbins_CC && j<2*nbins_CC) ECovars_e4_e4[i-nbins_CC][j-nbins_CC]     = covar_e4_e4;
	if(j>=2*nbins_CC)               ECovars_e4_nue[i-nbins_CC][j-2*nbins_CC]  = covar_e4_nue;}

      if(i>=2*nbins_CC) {
        if(j<nbins_CC)                  ECovars_nue_m4[i-2*nbins_CC][j]             = covar_nue_m4;
	if(j>=nbins_CC && j<2*nbins_CC) ECovars_nue_e4[i-2*nbins_CC][j-nbins_CC]    = covar_nue_e4;
	if(j>=2*nbins_CC)               ECovars_nue_nue[i-2*nbins_CC][j-2*nbins_CC] = covar_nue_nue;}
      


      if(i<nbins_CC) {
        if(j<nbins_CC)                  frECovars_m4_m4[i][j]             = covar_m4_m4  /(N * CC_m_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
	if(j>=nbins_CC && j<2*nbins_CC) frECovars_m4_e4[i][j-nbins_CC]    = covar_m4_e4  /(N * CC_m_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j-nbins_CC+1));
	if(j>=2*nbins_CC)               frECovars_m4_nue[i][j-2*nbins_CC] = covar_m4_nue /(N * CC_m_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}
	
      if(i>=nbins_CC && i<2*nbins_CC) {
        if(j<nbins_CC)                  frECovars_e4_m4[i-nbins_CC][j]             = covar_e4_m4  /(N * CC_e_nom->GetBinContent(i-nbins_CC+1) * CC_m_nom->GetBinContent(j+1));
	if(j>=nbins_CC && j<2*nbins_CC) frECovars_e4_e4[i-nbins_CC][j-nbins_CC]    = covar_e4_e4  /(N * CC_e_nom->GetBinContent(i-nbins_CC+1) * CC_e_nom->GetBinContent(j-nbins_CC+1));
	if(j>=2*nbins_CC)               frECovars_e4_nue[i-nbins_CC][j-2*nbins_CC] = covar_e4_nue /(N * CC_e_nom->GetBinContent(i-nbins_CC+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}

      if(i>=2*nbins_CC) {	
        if(j<nbins_CC)                  frECovars_nue_m4[i-2*nbins_CC][j]             = covar_nue_m4  /(N * nue_nom->GetBinContent(i-2*nbins_CC+1) * CC_m_nom->GetBinContent(j+1));
	if(j>=nbins_CC && j<2*nbins_CC) frECovars_nue_e4[i-2*nbins_CC][j-nbins_CC]    = covar_nue_e4  /(N * nue_nom->GetBinContent(i-2*nbins_CC+1) * CC_e_nom->GetBinContent(j-nbins_CC+1));
	if(j>=2*nbins_CC)               frECovars_nue_nue[i-2*nbins_CC][j-2*nbins_CC] = covar_nue_nue /(N * nue_nom->GetBinContent(i-2*nbins_CC+1) * nue_nom->GetBinContent(j-2*nbins_CC+1));}



      // TOTAL COVARIANCE MATRICES
      if(i<nbins_CC){
        if(j<nbins_CC)                  ECovars_mm[i][j]              = covar_mm  /N;
        if(j>=nbins_CC && j<2*nbins_CC) ECovars_me[i][j-nbins_CC]     = covar_me  /N;
        if(j>=2*nbins_CC)               ECovars_mnue[i][j-2*nbins_CC] = covar_mnue/N;}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  ECovars_em[i-nbins_CC][j]              = covar_em  /N;
        if(j>=nbins_CC && j<2*nbins_CC) ECovars_ee[i-nbins_CC][j-nbins_CC]     = covar_ee  /N;
        if(j>=2*nbins_CC)               ECovars_enue[i-nbins_CC][j-2*nbins_CC] = covar_enue/N;}

      if(i>=2*nbins_CC){
        if(j<nbins_CC)                  ECovars_nuem[i-2*nbins_CC][j]              = covar_nuem  /N;
        if(j>=nbins_CC && j<2*nbins_CC) ECovars_nuee[i-2*nbins_CC][j-nbins_CC]     = covar_nuee  /N;
        if(j>=2*nbins_CC)               ECovars_nuenue[i-2*nbins_CC][j-2*nbins_CC] = covar_nue_nue/N;}



      if(i<nbins_CC){
        if(j<nbins_CC)                  frECovars_mm[i][j]              = covar_mm  /(N * CC_m_nom->GetBinContent(i+1)  * CC_m_nom->GetBinContent(j+1));
        if(j>=nbins_CC && j<2*nbins_CC) frECovars_me[i][j-nbins_CC]     = covar_me  /(N * CC_m_nom->GetBinContent(i+1)  * CC_e_nom->GetBinContent(j-nbins_CC+1));
        if(j>=2*nbins_CC)               frECovars_mnue[i][j-2*nbins_CC] = covar_mnue/(N * CC_m_nom->GetBinContent(i+1)  * nue_nom->GetBinContent(j-2*nbins_CC+1));}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  frECovars_em[i-nbins_CC][j]              = covar_em  /(N * CC_e_nom->GetBinContent(i-nbins_CC+1)  * CC_m_nom->GetBinContent(j+1));
        if(j>=nbins_CC && j<2*nbins_CC) frECovars_ee[i-nbins_CC][j-nbins_CC]     = covar_ee  /(N * CC_e_nom->GetBinContent(i-nbins_CC+1)  * CC_e_nom->GetBinContent(j-nbins_CC+1));
        if(j>=2*nbins_CC)               frECovars_enue[i-nbins_CC][j-2*nbins_CC] = covar_enue/(N * CC_e_nom->GetBinContent(i-nbins_CC+1)  * nue_nom->GetBinContent(j-2*nbins_CC+1));}

      if(i>=2*nbins_CC){
        if(j<nbins_CC)                  frECovars_nuem[i-2*nbins_CC][j]            = covar_nuem  /(N * nue_nom->GetBinContent(i-2*nbins_CC+1)  * CC_m_nom->GetBinContent(j+1));
        if(j>=nbins_CC && j<2*nbins_CC) frECovars_nuee[i-2*nbins_CC][j-nbins_CC]   = covar_nuee  /(N * nue_nom->GetBinContent(i-2*nbins_CC+1)  * CC_e_nom->GetBinContent(j-nbins_CC+1));
        if(j>=2*nbins_CC)               ECovars_nuenue[i-2*nbins_CC][j-2*nbins_CC] = covar_nue_nue/(N * nue_nom->GetBinContent(i-2*nbins_CC+1)  * nue_nom->GetBinContent(j-2*nbins_CC+1));}

    }
  }

  for(int i = 0; i < nbins; i++) {
    for(int j = 0; j < nbins; j++) {
      if(i<nbins_CC){
        if(j<nbins_CC)                  ECovars2[i][j] = ECovars_m2_m2[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars2[i][j] = ECovars_m2_e2[i][j-nbins_CC];
        if(j>=2*nbins_CC)               ECovars2[i][j] = ECovars_m2_nue[i][j-2*nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  ECovars2[i][j] = ECovars_e2_m2[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars2[i][j] = ECovars_e2_e2[i-nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               ECovars2[i][j] = ECovars_e2_nue[i-nbins_CC][j-2*nbins_CC];}

      if(i>=2*nbins_CC){
        if(j<nbins_CC)                  ECovars2[i][j] = ECovars_nue_m2[i-2*nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars2[i][j] = ECovars_nue_e2[i-2*nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               ECovars2[i][j] = ECovars_nue_nue[i-2*nbins_CC][j-2*nbins_CC];}
      


      if(i<nbins_CC){
        if(j<nbins_CC)                  frECovars2[i][j] = frECovars_m2_m2[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars2[i][j] = frECovars_m2_e2[i][j-nbins_CC];
        if(j>=2*nbins_CC)               frECovars2[i][j] = frECovars_m2_nue[i][j-2*nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  frECovars2[i][j] = frECovars_e2_m2[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars2[i][j] = frECovars_e2_e2[i-nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               frECovars2[i][j] = frECovars_e2_nue[i-nbins_CC][j-2*nbins_CC];}

      if(i>=2*nbins_CC){
        if(j<nbins_CC)                  frECovars2[i][j] = frECovars_nue_m2[i-2*nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars2[i][j] = frECovars_nue_e2[i-2*nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               frECovars2[i][j] = frECovars_nue_nue[i-2*nbins_CC][j-2*nbins_CC];}
      




      if(i<nbins_CC){
        if(j<nbins_CC)                  ECovars[i][j] = ECovars_mm[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars[i][j] = ECovars_me[i][j-nbins_CC];
        if(j>=2*nbins_CC)               ECovars[i][j] = ECovars_mnue[i][j-2*nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  ECovars[i][j] = ECovars_em[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars[i][j] = ECovars_ee[i-nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               ECovars[i][j] = ECovars_enue[i-nbins_CC][j-2*nbins_CC];}

      if(i>=2*nbins_CC){
        if(j<nbins_CC)                  ECovars[i][j] = ECovars_nuem[i-2*nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars[i][j] = ECovars_nuee[i-2*nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               ECovars[i][j] = ECovars_nuenue[i-2*nbins_CC][j-2*nbins_CC];}
      


      if(i<nbins_CC){
        if(j<nbins_CC)                  frECovars[i][j] = frECovars_mm[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars[i][j] = frECovars_me[i][j-nbins_CC];
        if(j>=2*nbins_CC)               frECovars[i][j] = frECovars_mnue[i][j-2*nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  frECovars[i][j] = frECovars_em[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars[i][j] = frECovars_ee[i-nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               frECovars[i][j] = frECovars_enue[i-nbins_CC][j-2*nbins_CC];}

      if(i>=2*nbins_CC){
        if(j<nbins_CC)                  frECovars[i][j] = frECovars_nuem[i-2*nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars[i][j] = frECovars_nuee[i-2*nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               frECovars[i][j] = frECovars_nuenue[i-2*nbins_CC][j-2*nbins_CC];}






      if(i<nbins_CC){
        if(j<nbins_CC)                  ECovars4[i][j] = ECovars_m4_m4[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars4[i][j] = ECovars_m4_e4[i][j-nbins_CC];
        if(j>=2*nbins_CC)               ECovars4[i][j] = ECovars_m4_nue[i][j-2*nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  ECovars4[i][j] = ECovars_e4_m4[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars4[i][j] = ECovars_e4_e4[i-nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               ECovars4[i][j] = ECovars_e4_nue[i-nbins_CC][j-2*nbins_CC];}

      if(i>=2*nbins_CC){
        if(j<nbins_CC)                  ECovars4[i][j] = ECovars_nue_m4[i-2*nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars4[i][j] = ECovars_nue_e4[i-2*nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               ECovars4[i][j] = ECovars_nue_nue[i-2*nbins_CC][j-2*nbins_CC];}
      


      if(i<nbins_CC){
        if(j<nbins_CC)                  frECovars4[i][j] = frECovars_m4_m4[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars4[i][j] = frECovars_m4_e4[i][j-nbins_CC];
        if(j>=2*nbins_CC)               frECovars4[i][j] = frECovars_m4_nue[i][j-2*nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  frECovars4[i][j] = frECovars_e4_m4[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars4[i][j] = frECovars_e4_e4[i-nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               frECovars4[i][j] = frECovars_e4_nue[i-nbins_CC][j-2*nbins_CC];}

      if(i>=2*nbins_CC){
        if(j<nbins_CC)                  frECovars4[i][j] = frECovars_nue_m4[i-2*nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) frECovars4[i][j] = frECovars_nue_e4[i-2*nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               frECovars4[i][j] = frECovars_nue_nue[i-2*nbins_CC][j-2*nbins_CC];}
    }
  }




  TH2D *hcv2 = new TH2D("hcv2","",nbins,0,nbins,nbins,0,nbins);
  TH2D *hcv  = new TH2D("hcv","",nbins,0,nbins,nbins,0,nbins);
  TH2D *hcv4 = new TH2D("hcv4","",nbins,0,nbins,nbins,0,nbins);
  TH2D *hfrcv2 = new TH2D("hfrcv2","",nbins,0,nbins,nbins,0,nbins);
  TH2D *hfrcv  = new TH2D("hfrcv","",nbins,0,nbins,nbins,0,nbins);
  TH2D *hfrcv4 = new TH2D("hfrcv4","",nbins,0,nbins,nbins,0,nbins);

  for(int i=0; i<nbins; i++) {
    for(int j=0; j<nbins; j++) {
      hcv2->SetBinContent(i+1, j+1, ECovars2[i][j]);
      hcv->SetBinContent(i+1, j+1, ECovars[i][j]);
      hcv4->SetBinContent(i+1, j+1, ECovars4[i][j]);
      hfrcv2->SetBinContent(i+1, j+1, frECovars2[i][j]);
      hfrcv->SetBinContent(i+1, j+1, frECovars[i][j]);
      hfrcv4->SetBinContent(i+1, j+1, frECovars4[i][j]);
    }
  }


  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  double nu, Ev;
  if(cutNu == 0)      nu = 10.0;
  else if(cutNu == 3) nu = 0.3;

  hcv2->SetStats(0);
  hcv->SetStats(0);
  hcv4->SetStats(0);
  hfrcv2->SetStats(0);
  hfrcv->SetStats(0);
  hfrcv4->SetStats(0);

  hcv2->SetTitle(Form("%s Total Covariance (nu<%.1fGeV && Etheta2<3.0MeV) (-1sigma)",name[iw],nu));
  hcv->SetTitle(Form("%s Total Covariance (nu<%.1fGeV && Etheta2<3.0MeV)",name[iw],nu));
  hcv4->SetTitle(Form("%s Total Covariance (nu<%.1fGeV && Etheta2<3.0MeV) (+1sigma)",name[iw],nu));
  hfrcv2->SetTitle(Form("%s Fractional Covariance (nu<%.1fGeV && Etheta2<3.0MeV) (-1sigma)",name[iw],nu));
  hfrcv->SetTitle(Form("%s Fractional Covariance (nu<%.1fGeV && Etheta2<3.0MeV)",name[iw],nu));
  hfrcv4->SetTitle(Form("%s Fractional Covariance (nu<%.1fGeV && Etheta2<3.0MeV) (+1sigma)",name[iw],nu));

  hcv2->SetMinimum(0);
  hcv->SetMinimum(0);
  hcv4->SetMinimum(0);
  hfrcv2->SetMinimum(0);
  hfrcv->SetMinimum(0);
  hfrcv4->SetMinimum(0);

  TCanvas *ccv = new TCanvas("ccv","",1800,500);
  ccv->Divide(3,1);
  ccv->cd(1);
  gPad->SetRightMargin(0.15);
  hcv2->Draw("colz");
  ccv->cd(2);
  gPad->SetRightMargin(0.15);
  hcv->Draw("colz");
  ccv->cd(3);
  gPad->SetRightMargin(0.15);
  hcv4->Draw("colz");
  ccv->SaveAs(Form("%s_Covmx%d.png",name[iw],cutNu));
  
  TCanvas *ccvfr = new TCanvas("ccvfr","",1800,500);
  ccvfr->Divide(3,1);
  ccvfr->cd(1);
  gPad->SetRightMargin(0.15);
  hfrcv2->Draw("colz");
  ccvfr->cd(2);
  gPad->SetRightMargin(0.15);
  hfrcv->Draw("colz");
  ccvfr->cd(3);
  gPad->SetRightMargin(0.15);
  hfrcv4->Draw("colz");
  ccvfr->SaveAs(Form("%s_frCovmx%d.png",name[iw],cutNu));
  
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();


  TFile *out = new TFile(Form("%s_covmtr%d.root",name[iw],cutNu),"RECREATE");
  hcv2->Write();
  hcv->Write();
  hcv4->Write();
  hfrcv2->Write();
  hfrcv->Write();
  hfrcv4->Write();
  out->Close();

  }


}

/*

      //for( int k = 0; k < N; ++k ) { // rows





        var_m_i = (abs(var_m2_i) + abs(var_m4_i))/2;
        var_m_j = (abs(var_m2_j) + abs(var_m4_j))/2;
        var_e_i = (abs(var_e2_i) + abs(var_e4_i))/2;
        var_e_j = (abs(var_e2_j) + abs(var_e4_j))/2;





	

	
      if(i<nbins_CC) {
        if(j<nbins_CC)                  ECorrel_m2_m2[i][j]             = covar_m2_m2  /sqrt(var_m2_i * var_m2_j);
	if(j>=nbins_CC && j<2*nbins_CC) ECorrel_m2_e2[i][j-nbins_CC]    = covar_m2_e2  /sqrt(var_m2_i * var_e2_j);
	if(j>=2*nbins_CC)               ECorrel_m2_nue[i][j-2*nbins_CC] = covar_m2_nue /sqrt(var_m2_i * var_nue_j);}

      if(i>=nbins_CC && i<2*nbins_CC) {
	if(j<nbins_CC)                  ECorrel_e2_m2[i-nbins_CC][j]             = covar_e2_m2  /sqrt(var_e2_i * var_m2_j);
	if(j>=nbins_CC && j<2*nbins_CC) ECorrel_e2_e2[i-nbins_CC][j-nbins_CC]    = covar_e2_e2  /sqrt(var_e2_i * var_e2_j);
	if(j>=2*nbins_CC)               ECorrel_e2_nue[i-nbins_CC][j-2*nbins_CC] = covar_e2_nue /sqrt(var_e2_i * var_nue_j);}

      if(i>=2*nbins_CC) {
	if(j<nbins_CC)                  ECorrel_nue_m2[i-2*nbins_CC][j]             = covar_nue_m2  /sqrt(var_nue_i * var_m2_j);
	if(j>=nbins_CC && j<2*nbins_CC) ECorrel_nue_e2[i-2*nbins_CC][j-nbins_CC]    = covar_nue_e2  /sqrt(var_nue_i * var_e2_j);
	if(j>=2*nbins_CC)               ECorrel_nue_nue[i-2*nbins_CC][j-2*nbins_CC] = covar_nue_nue /sqrt(var_nue_i * var_nue_j);}
	


	ECovars_m4_m4[i][j]   = covar_m4_m4;
	ECovars_e4_e4[i][j]   = covar_e4_e4;
	
	ECovars_m4_e4[i][j]   = covar_m4_e4;
	ECovars_m4_nue[i][j]  = covar_m4_nue;
	ECovars_e4_m4[i][j]   = covar_e4_m4;
	ECovars_e4_nue[i][j]  = covar_e4_nue;
	ECovars_nue_m4[i][j]  = covar_nue_m4;
	ECovars_nue_e4[i][j]  = covar_nue_e4;
	
	frECovars_m4_m4[i][j]   = covar_m4_m4   /(N * CC_m_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
	frECovars_e4_e4[i][j]   = covar_e4_e4   /(N * CC_e_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
	
	frECovars_m4_e4[i][j]   = covar_m4_e4   /(N * CC_m_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
	frECovars_m4_nue[i][j]  = covar_m4_nue  /(N * CC_m_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
	frECovars_e4_m4[i][j]   = covar_e4_m4   /(N * CC_e_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
	frECovars_e4_nue[i][j]  = covar_e4_nue  /(N * CC_e_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
	frECovars_nue_m4[i][j]  = covar_nue_m4  /(N * nue_nom->GetBinContent(i+1)  * CC_m_nom->GetBinContent(j+1));
	frECovars_nue_e4[i][j]  = covar_nue_e4  /(N * nue_nom->GetBinContent(i+1)  * CC_e_nom->GetBinContent(j+1));
	
	
	
	ECorrel_m4_m4[i][j]   = covar_m4_m4   /sqrt(var_m4_i  * var_m4_j);
	ECorrel_e4_e4[i][j]   = covar_e4_e4   /sqrt(var_e4_i  * var_e4_j);
	
	ECorrel_m4_e4[i][j]   = covar_m4_e4   /sqrt(var_m4_i  * var_e4_j);
	ECorrel_m4_nue[i][j]  = covar_m4_nue  /sqrt(var_m4_i  * var_nue_j);
	ECorrel_e4_m4[i][j]   = covar_e4_m4   /sqrt(var_e4_i  * var_m4_j);
	ECorrel_e4_nue[i][j]  = covar_e4_nue  /sqrt(var_e4_i  * var_nue_j);
	ECorrel_nue_m4[i][j]  = covar_nue_m4  /sqrt(var_nue_i * var_m4_j);
	ECorrel_nue_e4[i][j]  = covar_nue_e4  /sqrt(var_nue_i * var_e4_j);


      ECovars_mm[i][j] = covar_mm;
      ECovars_ee[i][j] = covar_ee;
      ECovars_me[i][j]   = covar_me;
      ECovars_mnue[i][j]  = covar_mnue;
      ECovars_em[i][j]   = covar_em;
      ECovars_enue[i][j]  = covar_enue;
      ECovars_nuem[i][j]  = covar_nuem;
      ECovars_nuee[i][j]  = covar_nuee;
 
      frECovars_mm[i][j] = covar_mm / (N * CC_m_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      frECovars_ee[i][j] = covar_ee / (N * CC_e_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      frECovars_me[i][j]   = covar_me   /(N * CC_m_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      frECovars_mnue[i][j]  = covar_mnue  /(N * CC_m_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      frECovars_em[i][j]   = covar_em   /(N * CC_e_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      frECovars_enue[i][j]  = covar_enue  /(N * CC_e_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      frECovars_nuem[i][j]  = covar_nuem  /(N * nue_nom->GetBinContent(i+1)  * CC_m_nom->GetBinContent(j+1));
      frECovars_nuee[i][j]  = covar_nuee  /(N * nue_nom->GetBinContent(i+1)  * CC_e_nom->GetBinContent(j+1));

    }
  }

  for(int i = 0; i < nbins; i++) {
    for(int j = 0; j < nbins; j++) {
      if(i<nbins_CC){
        if(j<nbins_CC)                  ECovars2[i][j] = ECovars_m2_m2[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars2[i][j] = ECovars_m2_e2[i][j-nbins_CC];
        if(j>=2*nbins_CC)               ECovars2[i][j] = ECovars_m2_nue[i][j-2*nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  ECovars2[i][j] = ECovars_e2_m2[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars2[i][j] = ECovars_e2_e2[i-nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               ECovars2[i][j] = ECovars_e2_nue[i-nbins_CC][j-2*nbins_CC];}

      if(i>=2*nbins_CC){
        if(j<nbins_CC)                  ECovars2[i][j] = ECovars_nue_m2[i-2*nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECovars2[i][j] = ECovars_nue_e2[i-2*nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               ECovars2[i][j] = ECovars_nue_nue[i-2*nbins_CC][j-2*nbins_CC];}



      if(i<nbins_CC){
        if(j<nbins_CC)                  ECorrel2[i][j] = ECorrel_m2_m2[i][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECorrel2[i][j] = ECorrel_m2_e2[i][j-nbins_CC];
        if(j>=2*nbins_CC)               ECorrel2[i][j] = ECorrel_m2_nue[i][j-2*nbins_CC];}

      if(i>=nbins_CC && i<2*nbins_CC){
        if(j<nbins_CC)                  ECorrel2[i][j] = ECorrel_e2_m2[i-nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECorrel2[i][j] = ECorrel_e2_e2[i-nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               ECorrel2[i][j] = ECorrel_e2_nue[i-nbins_CC][j-2*nbins_CC];}

      if(i>=2*nbins_CC){
        if(j<nbins_CC)                  ECorrel2[i][j] = ECorrel_nue_m2[i-2*nbins_CC][j];
        if(j>=nbins_CC && j<2*nbins_CC) ECorrel2[i][j] = ECorrel_nue_e2[i-2*nbins_CC][j-nbins_CC];
        if(j>=2*nbins_CC)               ECorrel2[i][j] = ECorrel_nue_nue[i-2*nbins_CC][j-2*nbins_CC];}




      ECovars4[i][j]           = ECovars_m4_m4[i][j];
      ECovars4[i][j+nbins_E]   = ECovars_m4_e4[i][j];
      ECovars4[i][j+2*nbins_E] = ECovars_m4_nue[i][j];

      ECovars4[i+nbins_E][j]           = ECovars_e4_m4[i][j];
      ECovars4[i+nbins_E][j+nbins_E]   = ECovars_e4_e4[i][j];
      ECovars4[i+nbins_E][j+2*nbins_E] = ECovars_e4_nue[i][j];

      ECovars4[i+2*nbins_E][j]           = ECovars_nue_m4[i][j];
      ECovars4[i+2*nbins_E][j+nbins_E]   = ECovars_nue_e4[i][j];
      ECovars4[i+2*nbins_E][j+2*nbins_E] = ECovars_nue_nue[i][j];

      ECorrel4[i][j]           = ECorrel_m4_m4[i][j];
      ECorrel4[i][j+nbins_E]   = ECorrel_m4_e4[i][j];
      ECorrel4[i][j+2*nbins_E] = ECorrel_m4_nue[i][j];

      ECorrel4[i+nbins_E][j]           = ECorrel_e4_m4[i][j];
      ECorrel4[i+nbins_E][j+nbins_E]   = ECorrel_e4_e4[i][j];
      ECorrel4[i+nbins_E][j+2*nbins_E] = ECorrel_e4_nue[i][j];

      ECorrel4[i+2*nbins_E][j]           = ECorrel_nue_m4[i][j];
      ECorrel4[i+2*nbins_E][j+nbins_E]   = ECorrel_nue_e4[i][j];
      ECorrel4[i+2*nbins_E][j+2*nbins_E] = ECorrel_nue_nue[i][j];



      ECovars[i][j]           = ECovars_mm[i][j];
      ECovars[i][j+nbins_E]   = ECovars_me[i][j];
      ECovars[i][j+2*nbins_E] = ECovars_mnue[i][j];

      ECovars[i+nbins_E][j]           = ECovars_em[i][j];
      ECovars[i+nbins_E][j+nbins_E]   = ECovars_ee[i][j];
      ECovars[i+nbins_E][j+2*nbins_E] = ECovars_enue[i][j];

      ECovars[i+2*nbins_E][j]           = ECovars_nuem[i][j];
      ECovars[i+2*nbins_E][j+nbins_E]   = ECovars_nuee[i][j];
      ECovars[i+2*nbins_E][j+2*nbins_E] = ECovars_nue_nue[i][j];

    }
  }

 



  const Int_t Number = 3;
  Double_t Red[Number]    = { 0.00, 1.00, 1.00};
  Double_t Green[Number]  = { 0.00, 1.00, 0.00};
  Double_t Blue[Number]   = { 1.00, 1.00, 0.00};
  Double_t Length[Number] = { 0.00, 0.50, 1.00 };
  Int_t nb=999;
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);

  //hcr->GetZaxis()->SetRangeUser(-1., 1.);
  //hcr->SetContour(999);
  //TCanvas *ccr = new TCanvas("ccr","",800,600);
  //hcr2->Draw("colz");
  //ccr->SaveAs(Form("new_mtr/%s_Cor%d%d.png",name[iw],cutNu,cutEv));

  //gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();
*/
