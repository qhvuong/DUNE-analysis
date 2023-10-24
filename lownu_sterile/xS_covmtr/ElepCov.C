static const int N = 1; // number of universes
static const int nbins = 52;
static const int nbins_E = 100;

const int n_mu = 19; // number of muon bins in covariance
const int n_e = 7; // number of electron bins in covariance

// bin edges -- fix these to be whatever they actually are
const double mubins[20] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,8.,12.,16.,20.,40.,100.};
const double ebins[8] = {0.,2.,4.,6.,8.,10.,20.,100.};
  

TMatrixD covmx( nbins, nbins );
TMatrixD scaleCovars( nbins, nbins ); // pairwise covariance of columns of random numbers
TMatrixD scales( N, nbins );

TMatrixD E_m(N, nbins_E);
TMatrixD E_e(N, nbins_E);
TMatrixD E_m_nue(N, nbins_E);
TMatrixD E_e_nue(N, nbins_E);
TMatrixD E_nue(N, nbins_E);
  
TMatrixD ECovars_m2_m2   ( nbins_E, nbins_E );
TMatrixD ECovars_e2_e2   ( nbins_E, nbins_E );
//TMatrixD ECovars_m2_nue  ( nbins_E, nbins_E );
//TMatrixD ECovars_e2_nue  ( nbins_E, nbins_E );
TMatrixD ECovars_nue_nue ( nbins_E, nbins_E );

TMatrixD ECovars_m2_e2   ( nbins_E, nbins_E );
TMatrixD ECovars_m2_nue  ( nbins_E, nbins_E );
TMatrixD ECovars_e2_m2   ( nbins_E, nbins_E );
TMatrixD ECovars_e2_nue  ( nbins_E, nbins_E );
TMatrixD ECovars_nue_m2  ( nbins_E, nbins_E );
TMatrixD ECovars_nue_e2  ( nbins_E, nbins_E );


TMatrixD frECovars_m2_m2   ( nbins_E, nbins_E );
TMatrixD frECovars_e2_e2   ( nbins_E, nbins_E );
//TMatrixD ECovars_m2_nue  ( nbins_E, nbins_E );
//TMatrixD ECovars_e2_nue  ( nbins_E, nbins_E );
TMatrixD frECovars_nue_nue ( nbins_E, nbins_E );

TMatrixD frECovars_m2_e2   ( nbins_E, nbins_E );
TMatrixD frECovars_m2_nue  ( nbins_E, nbins_E );
TMatrixD frECovars_e2_m2   ( nbins_E, nbins_E );
TMatrixD frECovars_e2_nue  ( nbins_E, nbins_E );
TMatrixD frECovars_nue_m2  ( nbins_E, nbins_E );
TMatrixD frECovars_nue_e2  ( nbins_E, nbins_E );



TMatrixD ECorrel_m2_m2   ( nbins_E, nbins_E );
TMatrixD ECorrel_e2_e2   ( nbins_E, nbins_E );
TMatrixD ECorrel_nue_nue ( nbins_E, nbins_E );

TMatrixD ECorrel_m2_e2   ( nbins_E, nbins_E );
TMatrixD ECorrel_m2_nue  ( nbins_E, nbins_E );
TMatrixD ECorrel_e2_m2   ( nbins_E, nbins_E );
TMatrixD ECorrel_e2_nue  ( nbins_E, nbins_E );
TMatrixD ECorrel_nue_m2  ( nbins_E, nbins_E );
TMatrixD ECorrel_nue_e2  ( nbins_E, nbins_E );



TMatrixD ECovars_m4_m4   ( nbins_E, nbins_E );
TMatrixD ECovars_e4_e4   ( nbins_E, nbins_E );
//TMatrixD ECovars_m4_nue  ( nbins_E, nbins_E );
//TMatrixD ECovars_e4_nue  ( nbins_E, nbins_E );

TMatrixD ECovars_m4_e4   ( nbins_E, nbins_E );
TMatrixD ECovars_m4_nue  ( nbins_E, nbins_E );
TMatrixD ECovars_e4_m4   ( nbins_E, nbins_E );
TMatrixD ECovars_e4_nue  ( nbins_E, nbins_E );
TMatrixD ECovars_nue_m4  ( nbins_E, nbins_E );
TMatrixD ECovars_nue_e4  ( nbins_E, nbins_E );

TMatrixD frECovars_m4_m4   ( nbins_E, nbins_E );
TMatrixD frECovars_e4_e4   ( nbins_E, nbins_E );
//TMatrixD ECovars_m4_nue  ( nbins_E, nbins_E );
//TMatrixD ECovars_e4_nue  ( nbins_E, nbins_E );

TMatrixD frECovars_m4_e4   ( nbins_E, nbins_E );
TMatrixD frECovars_m4_nue  ( nbins_E, nbins_E );
TMatrixD frECovars_e4_m4   ( nbins_E, nbins_E );
TMatrixD frECovars_e4_nue  ( nbins_E, nbins_E );
TMatrixD frECovars_nue_m4  ( nbins_E, nbins_E );
TMatrixD frECovars_nue_e4  ( nbins_E, nbins_E );



TMatrixD ECorrel_m4_m4   ( nbins_E, nbins_E );
TMatrixD ECorrel_e4_e4   ( nbins_E, nbins_E );

TMatrixD ECorrel_m4_e4   ( nbins_E, nbins_E );
TMatrixD ECorrel_m4_nue  ( nbins_E, nbins_E );
TMatrixD ECorrel_e4_m4   ( nbins_E, nbins_E );
TMatrixD ECorrel_e4_nue  ( nbins_E, nbins_E );
TMatrixD ECorrel_nue_m4  ( nbins_E, nbins_E );
TMatrixD ECorrel_nue_e4  ( nbins_E, nbins_E );




TMatrixD ECovars_mm   ( nbins_E, nbins_E );
TMatrixD ECovars_ee   ( nbins_E, nbins_E );
//TMatrixD ECovars_nue_nue ( nbins_E, nbins_E );

TMatrixD ECovars_me   ( nbins_E, nbins_E );
TMatrixD ECovars_mnue  ( nbins_E, nbins_E );
TMatrixD ECovars_em   ( nbins_E, nbins_E );
TMatrixD ECovars_enue  ( nbins_E, nbins_E );
TMatrixD ECovars_nuem  ( nbins_E, nbins_E );
TMatrixD ECovars_nuee  ( nbins_E, nbins_E );

TMatrixD frECovars_mm   ( nbins_E, nbins_E );
TMatrixD frECovars_ee   ( nbins_E, nbins_E );
//TMatrixD ECovars_nue_nue ( nbins_E, nbins_E );

TMatrixD frECovars_me   ( nbins_E, nbins_E );
TMatrixD frECovars_mnue  ( nbins_E, nbins_E );
TMatrixD frECovars_em   ( nbins_E, nbins_E );
TMatrixD frECovars_enue  ( nbins_E, nbins_E );
TMatrixD frECovars_nuem  ( nbins_E, nbins_E );
TMatrixD frECovars_nuee  ( nbins_E, nbins_E );




TMatrixD ECovars  ( 3*nbins_E, 3*nbins_E );
TMatrixD frECovars  ( 3*nbins_E, 3*nbins_E );
TMatrixD ECorrel  ( 3*nbins_E, 3*nbins_E );

TMatrixD ECovars2 ( 3*nbins_E, 3*nbins_E );
TMatrixD frECovars2 ( 3*nbins_E, 3*nbins_E );
TMatrixD ECorrel2 ( 3*nbins_E, 3*nbins_E );

TMatrixD ECovars4 ( 3*nbins_E, 3*nbins_E );
TMatrixD frECovars4 ( 3*nbins_E, 3*nbins_E );
TMatrixD ECorrel4 ( 3*nbins_E, 3*nbins_E );

void ElepCov()
{

  int para, cutNu, cutEv;

  para=1;

  //for(para = 1; para <3; para++) {
  //TFile *f     = new TFile("/dune/app/users/qvuong/lownu/gen_data/new_tgt/total_output_032123.root");
  TFile *f     = new TFile("/dune/app/users/qvuong/data/lownu/CC_output.root");
  TFile *f_nue = new TFile("/dune/app/users/qvuong/data/lownu/nue_output.root");
  cutNu = 3;

  //std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};

  //const char *name[] = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"}; 

  std::list <const char *> namelist = {"wgt_MaCCQE"};
  const char *name[] = {"wgt_MaCCQE"}; 

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

  TH1D *Ev_m1       = (TH1D*)f->Get(Form("%s_m_hEv%d_sigma1",name[iw],cutNu));
  TH1D *Ev_e1       = (TH1D*)f->Get(Form("%s_e_hEv%d_sigma1",name[iw],cutNu));
  TH1D *Ev_m2       = (TH1D*)f->Get(Form("%s_m_hEv%d_sigma2",name[iw],cutNu));
  TH1D *Ev_e2       = (TH1D*)f->Get(Form("%s_e_hEv%d_sigma2",name[iw],cutNu));
  TH1D *Ev_CC_m_nom = (TH1D*)f->Get(Form("%s_m_hEv%d_sigma3",name[iw],cutNu));
  TH1D *Ev_CC_e_nom = (TH1D*)f->Get(Form("%s_e_hEv%d_sigma3",name[iw],cutNu));
  TH1D *Ev_m4       = (TH1D*)f->Get(Form("%s_m_hEv%d_sigma4",name[iw],cutNu));
  TH1D *Ev_e4       = (TH1D*)f->Get(Form("%s_e_hEv%d_sigma4",name[iw],cutNu));
  TH1D *Ev_m5       = (TH1D*)f->Get(Form("%s_m_hEv%d_sigma5",name[iw],cutNu));
  TH1D *Ev_e5       = (TH1D*)f->Get(Form("%s_e_hEv%d_sigma5",name[iw],cutNu));
  //TH1D *Ev_nue      = (TH1D*)f_nue->Get(Form("hEv%d",cutEv));
  //TH1D *Ev_nue_nom  = (TH1D*)f_nue->Get(Form("hEv%d",cutEv));

  TH1D *h1 = (TH1D*)e1->Clone();
  TH1D *h2 = (TH1D*)e2->Clone();
  TH1D *h4 = (TH1D*)e4->Clone();
  TH1D *h5 = (TH1D*)e5->Clone();

  TH1D *Ev_h1 = (TH1D*)Ev_e1->Clone();
  TH1D *Ev_h2 = (TH1D*)Ev_e2->Clone();
  TH1D *Ev_h4 = (TH1D*)Ev_e4->Clone();
  TH1D *Ev_h5 = (TH1D*)Ev_e5->Clone();

/*
  TCanvas *c = new TCanvas("c","",900,1100);
  TPad *pad1 = new TPad("pad1","pad1", 0, 0.4, 1, 1.0);
  TPad *pad2 = new TPad("pad2","pad2", 0, 0.05, 1, 0.4);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->SetGridy();         // Vertical grid
  pad1->Draw();
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad2->SetGridx(); // vertical grid
  pad2->SetGridy();         // Vertical grid
  pad2->Draw();

  c->cd();
  pad1->cd();
  CC_e_nom->SetTitle(Form("CC_e Elep in %s universe",name[iw]));
  CC_e_nom->SetMarkerColor(kBlack);
  CC_e_nom->SetMarkerStyle(8);
  CC_e_nom->SetMarkerSize(0.8);
  CC_e_nom->SetMaximum(9E3);
  CC_e_nom->Draw("same");
  e1->SetMarkerColor(kBlue);
  e1->SetMarkerStyle(21);
  e1->SetMarkerSize(0.5);
  e1->Draw("same");
  e2->SetMarkerColor(kBlue);
  e2->SetMarkerStyle(21);
  e2->SetMarkerSize(0.5);
  e2->Draw("same");
  e4->SetMarkerColor(kRed);
  e4->SetMarkerStyle(21);
  e4->SetMarkerSize(0.5);
  e4->Draw("same");
  e5->SetMarkerColor(kRed);
  e5->SetMarkerStyle(21);
  e5->SetMarkerSize(0.5);
  e5->Draw("same");

  c->cd();
  pad2->cd();
  TH1D *ratio1 = (TH1D*)h1->Clone("ratio1");
  ratio1->Divide(CC_e_nom);
  ratio1->SetLineColor(kBlue);
  ratio1->SetMaximum(1.4);
  ratio1->SetMinimum(0.7);
  ratio1->SetStats(0); 
  ratio1->GetYaxis()->SetTitle("ratio Elep/Elep_nom");
  ratio1->GetXaxis()->SetTitle("Elep (GeV)");
  ratio1->Draw("same");
  TH1D *ratio2 = (TH1D*)h2->Clone("ratio2");
  ratio2->Divide(CC_e_nom);
  ratio2->SetLineColor(kBlue);
  ratio2->SetStats(0); 
  ratio2->Draw("same");
  TH1D *ratio4 = (TH1D*)h4->Clone("ratio4");
  ratio4->Divide(CC_e_nom);
  ratio4->SetLineColor(kRed);
  ratio4->SetStats(0); 
  ratio4->Draw("same");
  TH1D *ratio5 = (TH1D*)h5->Clone("ratio5");
  ratio5->Divide(CC_e_nom);
  ratio5->SetLineColor(kRed);
  ratio5->SetStats(0); 
  ratio5->Draw("same");
  c->SaveAs(Form("nue_%s_ratio.png",name[iw]));

  //TFile *out = new TFile(Form("nue_%s_ratio.root"),"RECREATE");
*/  

/*
  TH1D *ratio = (TH1D*)h1->Clone("ratio");
  ratio->SetTitle(Form("CC_m Elep in %d universes (-2sigma)",N_wgt));
  ratio->Divide(CC_m_nom);
  ratio->SetLineColor(iw+1);
  ratio->SetMaximum(1.4);
  ratio->SetMinimum(0.7);
  ratio->SetStats(0); 
  ratio->GetYaxis()->SetTitle("ratio Elep/Elep_nom");
  ratio->GetXaxis()->SetTitle("Elep (GeV)");
  ratio->Draw("same");
*/

  for( int i = 0; i < nbins_E; ++i ) { // columns
    for( int j = 0; j < nbins_E; ++j ) { // columns
      // compute column covariance
      double covar_m2_m2 = 0.;
      double covar_e2_e2 = 0.;
      //double covar_m2_nue = 0.;
      //double covar_e2_nue = 0.;
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
      //double covar_m4_nue = 0.;
      //double covar_e4_nue = 0.;
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
        covar_nue_nue += (nue->GetBinContent(i+1) - nue_nom->GetBinContent(i+1))  * (nue->GetBinContent(j+1) - nue_nom->GetBinContent(j+1)); 

        covar_m2_e2   += (m2->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (e2->GetBinContent(j+1)  - CC_e_nom->GetBinContent(j+1));
        covar_m2_nue  += (m2->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (nue->GetBinContent(j+1) - nue_nom->GetBinContent(j+1));
        covar_e2_m2   += (e2->GetBinContent(i+1)  - CC_e_nom->GetBinContent(i+1)) * (m2->GetBinContent(j+1)  - CC_m_nom->GetBinContent(j+1)); 
        covar_e2_nue  += (e2->GetBinContent(i+1)  - CC_e_nom->GetBinContent(i+1)) * (nue->GetBinContent(j+1) - nue_nom->GetBinContent(j+1)); 
        covar_nue_m2  += (nue->GetBinContent(i+1) - nue_nom->GetBinContent(i+1))  * (m2->GetBinContent(j+1)  - CC_m_nom->GetBinContent(j+1)); 
        covar_nue_e2  += (nue->GetBinContent(i+1) - nue_nom->GetBinContent(i+1))  * (e2->GetBinContent(j+1)  - CC_e_nom->GetBinContent(j+1)); 

        var_m2_i += (m2->GetBinContent(i+1) - CC_m_nom->GetBinContent(i+1)) * (m2->GetBinContent(i+1) - CC_m_nom->GetBinContent(i+1));
        var_m2_j += (m2->GetBinContent(j+1) - CC_m_nom->GetBinContent(j+1)) * (m2->GetBinContent(j+1) - CC_m_nom->GetBinContent(j+1));

        var_e2_i += (e2->GetBinContent(i+1) - CC_e_nom->GetBinContent(i+1)) * (e2->GetBinContent(i+1) - CC_e_nom->GetBinContent(i+1));
        var_e2_j += (e2->GetBinContent(j+1) - CC_e_nom->GetBinContent(j+1)) * (e2->GetBinContent(j+1) - CC_e_nom->GetBinContent(j+1));

        var_nue_i += (nue->GetBinContent(i+1) - nue_nom->GetBinContent(i+1)) * (nue->GetBinContent(i+1) - nue_nom->GetBinContent(i+1));
        var_nue_j += (nue->GetBinContent(j+1) - nue_nom->GetBinContent(j+1)) * (nue->GetBinContent(j+1) - nue_nom->GetBinContent(j+1));


        covar_m4_m4   += (m4->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (m4->GetBinContent(j+1)  - CC_m_nom->GetBinContent(j+1));
        covar_e4_e4   += (e4->GetBinContent(i+1)  - CC_e_nom->GetBinContent(i+1)) * (e4->GetBinContent(j+1)  - CC_e_nom->GetBinContent(j+1)); 

        covar_m4_e4   += (m4->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (e4->GetBinContent(j+1)  - CC_e_nom->GetBinContent(j+1));
        covar_m4_nue  += (m4->GetBinContent(i+1)  - CC_m_nom->GetBinContent(i+1)) * (nue->GetBinContent(j+1) - nue_nom->GetBinContent(j+1));
        covar_e4_m4   += (e4->GetBinContent(i+1)  - CC_e_nom->GetBinContent(i+1)) * (m4->GetBinContent(j+1)  - CC_m_nom->GetBinContent(j+1)); 
        covar_e4_nue  += (e4->GetBinContent(i+1)  - CC_e_nom->GetBinContent(i+1)) * (nue->GetBinContent(j+1) - nue_nom->GetBinContent(j+1)); 
        covar_nue_m4  += (nue->GetBinContent(i+1) - nue_nom->GetBinContent(i+1))  * (m4->GetBinContent(j+1)  - CC_m_nom->GetBinContent(j+1)); 
        covar_nue_e4  += (nue->GetBinContent(i+1) - nue_nom->GetBinContent(i+1))  * (e4->GetBinContent(j+1)  - CC_e_nom->GetBinContent(j+1)); 

        var_m4_i += (m4->GetBinContent(i+1) - CC_m_nom->GetBinContent(i+1)) * (m4->GetBinContent(i+1) - CC_m_nom->GetBinContent(i+1));
        var_m4_j += (m4->GetBinContent(j+1) - CC_m_nom->GetBinContent(j+1)) * (m4->GetBinContent(j+1) - CC_m_nom->GetBinContent(j+1));
       
        var_e4_i += (e4->GetBinContent(i+1) - CC_e_nom->GetBinContent(i+1)) * (e4->GetBinContent(i+1) - CC_e_nom->GetBinContent(i+1));
        var_e4_j += (e4->GetBinContent(j+1) - CC_e_nom->GetBinContent(j+1)) * (e4->GetBinContent(j+1) - CC_e_nom->GetBinContent(j+1));


        covar_mm = (abs(covar_m2_m2) + abs(covar_m4_m4))/2;
        covar_ee = (abs(covar_e2_e2) + abs(covar_e4_e4))/2;

        covar_me = (abs(covar_m2_e2) + abs(covar_m4_e4))/2;
        covar_mnue = (abs(covar_m2_nue) + abs(covar_m4_nue))/2;
        covar_em = (abs(covar_e2_m2) + abs(covar_e4_m4))/2;
        covar_enue = (abs(covar_e2_nue) + abs(covar_e4_nue))/2;
        covar_nuem = (abs(covar_nue_m2) + abs(covar_nue_m4))/2;
        covar_nuee = (abs(covar_nue_e2) + abs(covar_nue_e4))/2;

/*
        var_m_i = (abs(var_m2_i) + abs(var_m4_i))/2;
        var_m_j = (abs(var_m2_j) + abs(var_m4_j))/2;
        var_e_i = (abs(var_e2_i) + abs(var_e4_i))/2;
        var_e_j = (abs(var_e2_j) + abs(var_e4_j))/2;
*/
      }
      ECovars_m2_m2[i][j]   = covar_m2_m2;
      ECovars_e2_e2[i][j]   = covar_e2_e2;
      ECovars_nue_nue[i][j] = covar_nue_nue;

      ECovars_m2_e2[i][j]   = covar_m2_e2;
      ECovars_m2_nue[i][j]  = covar_m2_nue;
      ECovars_e2_m2[i][j]   = covar_e2_m2;
      ECovars_e2_nue[i][j]  = covar_e2_nue;
      ECovars_nue_m2[i][j]  = covar_nue_m2;
      ECovars_nue_e2[i][j]  = covar_nue_e2;

      frECovars_m2_m2[i][j]   = covar_m2_m2   /(N * CC_m_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      frECovars_e2_e2[i][j]   = covar_e2_e2   /(N * CC_e_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      frECovars_nue_nue[i][j] = covar_nue_nue /(N * nue_nom->GetBinContent(i+1)  * nue_nom->GetBinContent(j+1));

      frECovars_m2_e2[i][j]   = covar_m2_e2   /(N * CC_m_nom->GetBinContent(i+1) * CC_e_nom->GetBinContent(j+1));
      frECovars_m2_nue[i][j]  = covar_m2_nue  /(N * CC_m_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      frECovars_e2_m2[i][j]   = covar_e2_m2   /(N * CC_e_nom->GetBinContent(i+1) * CC_m_nom->GetBinContent(j+1));
      frECovars_e2_nue[i][j]  = covar_e2_nue  /(N * CC_e_nom->GetBinContent(i+1) * nue_nom->GetBinContent(j+1));
      frECovars_nue_m2[i][j]  = covar_nue_m2  /(N * nue_nom->GetBinContent(i+1)  * CC_m_nom->GetBinContent(j+1));
      frECovars_nue_e2[i][j]  = covar_nue_e2  /(N * nue_nom->GetBinContent(i+1)  * CC_e_nom->GetBinContent(j+1));

      ECorrel_m2_m2[i][j]   = covar_m2_m2   /sqrt(var_m2_i  * var_m2_j);
      ECorrel_e2_e2[i][j]   = covar_e2_e2   /sqrt(var_e2_i  * var_e2_j);
      ECorrel_nue_nue[i][j] = covar_nue_nue /sqrt(var_nue_i * var_nue_j);

      ECorrel_m2_e2[i][j]   = covar_m2_e2   /sqrt(var_m2_i  * var_e2_j);
      ECorrel_m2_nue[i][j]  = covar_m2_nue  /sqrt(var_m2_i  * var_nue_j);
      ECorrel_e2_m2[i][j]   = covar_e2_m2   /sqrt(var_e2_i  * var_m2_j);
      ECorrel_e2_nue[i][j]  = covar_e2_nue  /sqrt(var_e2_i  * var_nue_j);
      ECorrel_nue_m2[i][j]  = covar_nue_m2  /sqrt(var_nue_i * var_m2_j);
      ECorrel_nue_e2[i][j]  = covar_nue_e2  /sqrt(var_nue_i * var_e2_j);


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

  for(int i = 0; i < nbins_E; i++) {
    for(int j = 0; j < nbins_E; j++) {
      ECovars2[i][j]           = ECovars_m2_m2[i][j];
      ECovars2[i][j+nbins_E]   = ECovars_m2_e2[i][j];
      ECovars2[i][j+2*nbins_E] = ECovars_m2_nue[i][j];

      ECovars2[i+nbins_E][j]           = ECovars_e2_m2[i][j];
      ECovars2[i+nbins_E][j+nbins_E]   = ECovars_e2_e2[i][j];
      ECovars2[i+nbins_E][j+2*nbins_E] = ECovars_e2_nue[i][j];

      ECovars2[i+2*nbins_E][j]           = ECovars_nue_m2[i][j];
      ECovars2[i+2*nbins_E][j+nbins_E]   = ECovars_nue_e2[i][j];
      ECovars2[i+2*nbins_E][j+2*nbins_E] = ECovars_nue_nue[i][j];

      ECorrel2[i][j]           = ECorrel_m2_m2[i][j];
      ECorrel2[i][j+nbins_E]   = ECorrel_m2_e2[i][j];
      ECorrel2[i][j+2*nbins_E] = ECorrel_m2_nue[i][j];

      ECorrel2[i+nbins_E][j]           = ECorrel_e2_m2[i][j];
      ECorrel2[i+nbins_E][j+nbins_E]   = ECorrel_e2_e2[i][j];
      ECorrel2[i+nbins_E][j+2*nbins_E] = ECorrel_e2_nue[i][j];

      ECorrel2[i+2*nbins_E][j]           = ECorrel_nue_m2[i][j];
      ECorrel2[i+2*nbins_E][j+nbins_E]   = ECorrel_nue_e2[i][j];
      ECorrel2[i+2*nbins_E][j+2*nbins_E] = ECorrel_nue_nue[i][j];



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

  double max=0., max2=0., max4=0.;
  double min=1., min2=1., min4=1.;

  for(int i = 0; i < 3*nbins_E; i++) {
    for(int j = 0; j < 3*nbins_E; j++) {
      //if(ECovars2[i][j] < 0.) std::cout << "ECovars2 = " << ECovars2[i][j] << "\n";
      //if(ECovars4[i][j] < 0.) std::cout << "ECovars4 = " << ECovars4[i][j] << "\n";
      if(max2 > ECovars2[i][j]) max2 = ECovars2[i][j];
      if(max4 > ECovars4[i][j]) max4 = ECovars4[i][j];
      if(min2 < ECovars2[i][j]) min2 = ECovars2[i][j];
      if(min4 < ECovars4[i][j]) min4 = ECovars4[i][j];
      //ECovars[i][j] = (ECovars2[i][j] + ECovars4[i][j])/2;
      ECorrel[i][j] = (ECorrel2[i][j] + ECorrel4[i][j])/2;
    }
  }

 
  //std::cout << min2 << "\t" << max2 << "\n" << min4 << "\t" << max4 << "\n";

  TH2D *hcv  = new TH2D("hcv","",300,0,300,300,0,300);
  TH2D *hcv2 = new TH2D("hcv2","",300,0,300,300,0,300);
  TH2D *hcv4 = new TH2D("hcv4","",300,0,300,300,0,300);
  TH2D *hfrcv  = new TH2D("hfrcv","",300,0,300,300,0,300);
  TH2D *hfrcv2 = new TH2D("hfrcv2","",300,0,300,300,0,300);
  TH2D *hfrcv4 = new TH2D("hfrcv4","",300,0,300,300,0,300);
  TH2D *hcr  = new TH2D("hcr","",300,0,300,300,0,300);
  TH2D *hcr2 = new TH2D("hcr2","",300,0,300,300,0,300);
  TH2D *hcr4 = new TH2D("hcr4","",300,0,300,300,0,300);

  for(int i=0; i<300; i++) {
    for(int j=0; j<300; j++) {
      hcv->SetBinContent(i+1, j+1, ECovars[i][j]);
      hcv2->SetBinContent(i+1, j+1, ECovars2[i][j]);
      hcv4->SetBinContent(i+1, j+1, ECovars4[i][j]);
      hfrcv->SetBinContent(i+1, j+1, frECovars[i][j]);
      hfrcv2->SetBinContent(i+1, j+1, frECovars2[i][j]);
      hfrcv4->SetBinContent(i+1, j+1, frECovars4[i][j]);
      hcr->SetBinContent(i+1, j+1, ECorrel[i][j]);
      hcr2->SetBinContent(i+1, j+1, ECorrel2[i][j]);
      hcr4->SetBinContent(i+1, j+1, ECorrel4[i][j]);
    }
  }


  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  double nu, Ev;
  if(cutNu == 0)      nu = 10.0;
  else if(cutNu == 3) nu = 0.3;
  if(cutEv == 0)      Ev = 3.0;
  else if(cutEv == 1) Ev = 0.8;
  else if(cutEv == 2) Ev = 0.5;

  hcv->SetStats(0);
  hcv->SetTitle(Form("%s Cross Section Covariance (nu<%.1fGeV && Etheta2<%.1fMeV)",name[iw],nu,Ev));
  hcr->SetStats(0);
  hcr->SetTitle(Form("%s Cross Section Correlation (nu<%.1fGeV && Etheta2<%.1fMeV)",name[iw],nu,Ev));

  hcv2->SetStats(0);
  hcv2->SetTitle(Form("%s Cross Section Covariance (nu<%.1fGeV && Etheta2<%.1fMeV) (-1sigma)",name[iw],nu,Ev));
  hcr2->SetStats(0);
  hcr2->SetTitle(Form("%s Cross Section Correlation (nu<%.1fGeV && Etheta2<%.1fMeV) (-1sigma)",name[iw],nu,Ev));

  hcv4->SetStats(0);
  hcv4->SetTitle(Form("%s Cross Section Covariance (nu<%.1fGeV && Etheta2<%.1fMeV) (+1sigma)",name[iw],nu,Ev));
  hcr4->SetStats(0);
  hcr4->SetTitle(Form("%s Cross Section Correlation (nu<%.1fGeV && Etheta2<%.1fMeV) (+1sigma)",name[iw],nu,Ev));

  hcv2->SetMinimum(0);
  hcv->SetMinimum(0);
  hcv4->SetMinimum(0);
  hfrcv2->SetMinimum(0);
  hfrcv->SetMinimum(0);
  hfrcv4->SetMinimum(0);

  TCanvas *c2 = new TCanvas("c2","",800,600);
  hcv2->Draw("colz");
  //c2->SaveAs(Form("new_mtr/%s_sigma2Cov%d%d.png",name[iw],cutNu,cutEv));
  TCanvas *c3 = new TCanvas("c3","",800,600);
  hcv->Draw("colz");
  //c3->SaveAs(Form("new_mtr/%s_sigma3Cov%d%d.png",name[iw],cutNu,cutEv));
  TCanvas *c4 = new TCanvas("c4","",800,600);
  hcv4->Draw("colz");
  //c4->SaveAs(Form("new_mtr/%s_sigma4Cov%d%d.png",name[iw],cutNu,cutEv));
/*
  TCanvas *c2fr = new TCanvas("c2fr","",800,600);
  hfrcv2->Draw("colz");
  //c2->SaveAs(Form("new_mtr/%s_sigma2Cov%d%d.png",name[iw],cutNu,cutEv));
  TCanvas *c3fr = new TCanvas("c3fr","",800,600);
  hfrcv->Draw("colz");
  //c3->SaveAs(Form("new_mtr/%s_sigma3Cov%d%d.png",name[iw],cutNu,cutEv));
  TCanvas *c4fr = new TCanvas("c4fr","",800,600);
  hfrcv4->Draw("colz");
  //c4->SaveAs(Form("new_mtr/%s_sigma4Cov%d%d.png",name[iw],cutNu,cutEv));
*/


/*

  const Int_t Number = 3;
  Double_t Red[Number]    = { 0.00, 1.00, 1.00};
  Double_t Green[Number]  = { 0.00, 1.00, 0.00};
  Double_t Blue[Number]   = { 1.00, 1.00, 0.00};
  Double_t Length[Number] = { 0.00, 0.50, 1.00 };
  Int_t nb=999;
  TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
*/
  //hcr->GetZaxis()->SetRangeUser(-1., 1.);
  //hcr->SetContour(999);
  //TCanvas *ccr = new TCanvas("ccr","",800,600);
  //hcr2->Draw("colz");
  //ccr->SaveAs(Form("new_mtr/%s_Cor%d%d.png",name[iw],cutNu,cutEv));

  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  //TFile *out = new TFile(Form("/dune/app/users/qvuong/lownu/sigma_mtr/new_mtr/%s_Covmtr%d%d.root",name[iw],cutNu,cutEv),"RECREATE");
  TFile *out = new TFile(Form("%s_Covmtr%d%d.root",name[iw],cutNu,cutEv),"RECREATE");
  hcv2->Write();
  //hcr2->Write();
  hcv->Write();
  //hcr->Write();
  hcv4->Write();
  //hcr4->Write();
  out->Close();

  //}
  //}

  }

  //c->SaveAs("numu1_ratio.png");

}

/*
  TDecompSVD svd0(ElepCovars0);
  //svd0.Decompose();
  TMatrixD inv0 = svd0.Invert();
  TDecompSVD svd3(ElepCovars3);
  //svd3.Decompose();
  TMatrixD inv3 = svd3.Invert();

  TH2D *test0 = new TH2D(inv0);
  TH2D *test3 = new TH2D(inv3);
  TCanvas *c2 = new TCanvas("c2","",1800,1350);
  c2->Divide(2,1);
  c2->cd(1);
  test0->Draw("colz");
  c2->cd(2);
  test3->Draw("colz");
  c2->SaveAs(Form("test1_total_Cov_%d.png",N));
*/
/*
    if(u<20) {

      c->cd();
      pad1->cd();
      TH1D *CC_m0 = (TH1D*)nue->Clone();
      CC_m0->SetTitle("Nue Elep in 20 universes (nu<10.0)");
      CC_m0->GetYaxis()->SetTitle("Elep histograms in 20 universes");
      //CC_m0->GetYaxis()->SetTitleSize(20);
      CC_m0->SetLineColor(u+1);
      CC_m0->SetStats(0);
      CC_m0->Draw("same");
      
      c->cd();
      pad2->cd();
      TH1D *ratio_m0 = (TH1D*)nue->Clone("ratio");
      ratio_m0->SetLineColor(u+1);
      ratio_m0->SetMaximum(1.7);
      ratio_m0->SetMinimum(0.3);
      ratio_m0->Sumw2(); 
      ratio_m0->SetStats(0); 
      ratio_m0->Divide(nue_nom);
      ratio_m0->SetTitle("");
      ratio_m0->GetYaxis()->SetTitle("ratio Elep/Elep_nom");
      ratio_m0->GetXaxis()->SetTitle("Elep (GeV)");
      ratio_m0->Draw("same");

      c_m3->cd();
      pad1->cd();
      TH1D *CC_m3 = (TH1D*)m3->Clone();
      CC_m3->SetTitle("CC_m Elep in 20 universes (nu<0.3)");
      CC_m3->GetYaxis()->SetTitle("Elep histograms in 20 universes");
      //CC_m0->GetYaxis()->SetTitleSize(20);
      CC_m3->SetLineColor(u+1);
      CC_m3->SetStats(0);
      CC_m3->Draw("same");
      
      c_m3->cd();
      pad2->cd();
      TH1D *ratio_m3 = (TH1D*)m3->Clone("ratio");
      ratio_m3->SetLineColor(u+1);
      ratio_m3->SetMaximum(1.7);
      ratio_m3->SetMinimum(0.3);
      ratio_m3->Sumw2(); 
      ratio_m3->SetStats(0); 
      ratio_m3->Divide(CC_m3_nom);
      ratio_m3->SetTitle("");
      ratio_m3->GetYaxis()->SetTitle("ratio Elep/Elep_nom");
      ratio_m3->GetXaxis()->SetTitle("Elep (GeV)");
      ratio_m3->Draw("same");

    }
    else continue;
*/
