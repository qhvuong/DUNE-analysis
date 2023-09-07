#include <iostream>
#include <string>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TTree.h>
#include <TRandom3.h>
#include <THStack.h>
#include <TChain.h>
#include <TLegend.h>
#include <list>

int main()
{

  // Load the CAF file = Common Analysis Format, a standard TTree that we use in DUNE
  TChain * tree = new TChain( "cafTree", "cafTree" );
  TChain * meta = new TChain( "meta", "meta" );
  
  for(int i = 0; i<400; i++){
  tree->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%02d.root",i) );
  meta->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%02d.root",i) ); // make certain this is the exact same file(s)
  std::cout << "File number:" << i << "\n";
  } 

  double total_pot = 0.;
  double pot;
  meta->SetBranchAddress( "pot", &pot );
  const int Nfiles = meta->GetEntries();
  for( int ii = 0; ii < Nfiles; ++ii ) {
    meta->GetEntry(ii);
    total_pot += pot;
  }
  double yrPOT = 1.1E21;
  double scalePOT = yrPOT/total_pot;

  const Int_t nbinsX = 430; const Int_t nbinsY = 100;
  Double_t xEdges[nbinsX+1], yEdges[nbinsY+1];
  xEdges[0]=yEdges[0]=0;
  for(int i=0; i<nbinsX+1; i++)
  {
    if(i<200)                xEdges[i+1] = xEdges[i] + 0.02;
    else if(i>=200 && i<240) xEdges[i+1] = xEdges[i] + 0.1;
    else if(i>=240 && i<400) xEdges[i+1] = xEdges[i] + 0.2;
    else if(i>=400 && i<420) xEdges[i+1] = xEdges[i] + 1.0;
    else                     xEdges[i+1] = xEdges[i] + 4.0;
  }
  for(int i=0; i<nbinsY+1; i++)
  {
    yEdges[i+1]=yEdges[i]+0.16;
  }

  const double mubins[20] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,8.,12.,16.,20.,40.,100.};
  const double ebins[8] = {0.,2.,4.,6.,8.,10.,20.,100.};

  TH2D *m_hElepRecoVsEv0    = new TH2D("m_hElepRecoVsEv0","",   nbinsX,xEdges,nbinsY,yEdges);
  TH2D *m_hElepRecoVsEv3    = new TH2D("m_hElepRecoVsEv3","",   nbinsX,xEdges,nbinsY,yEdges);
  TH2D *nc_m_hElepRecoVsEv0 = new TH2D("nc_m_hElepRecoVsEv0","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *nc_m_hElepRecoVsEv3 = new TH2D("nc_m_hElepRecoVsEv3","",nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepRecoVsEv0    = new TH2D("e_hElepRecoVsEv0","",   nbinsX,xEdges,nbinsY,yEdges);
  TH2D *e_hElepRecoVsEv3    = new TH2D("e_hElepRecoVsEv3","",   nbinsX,xEdges,nbinsY,yEdges);

  TH2D *m_hElepRecoVsEv0_cov = new TH2D("m_hElepRecoVsEv0_cov","",19,mubins,100,0,16);
  TH2D *m_hElepRecoVsEv3_cov = new TH2D("m_hElepRecoVsEv3_cov","",19,mubins,100,0,16);
  TH2D *e_hElepRecoVsEv0_cov = new TH2D("e_hElepRecoVsEv0_cov","",7, ebins,100,0,16);
  TH2D *e_hElepRecoVsEv3_cov = new TH2D("e_hElepRecoVsEv3_cov","",7, ebins,100,0,16);

  TH2D *m_hElepRecoVsEv0_cov_sigma1 = new TH2D("m_hElepRecoVsEv0_cov_sigma1","",19,mubins,100,0,16);
  TH2D *m_hElepRecoVsEv3_cov_sigma1 = new TH2D("m_hElepRecoVsEv3_cov_sigma1","",19,mubins,100,0,16);
  TH2D *e_hElepRecoVsEv0_cov_sigma1 = new TH2D("e_hElepRecoVsEv0_cov_sigma1","",19,mubins,100,0,16);
  TH2D *e_hElepRecoVsEv3_cov_sigma1 = new TH2D("e_hElepRecoVsEv3_cov_sigma1","",19,mubins,100,0,16);

  TH2D *m_hElepRecoVsEv0_cov_sigma2 = new TH2D("m_hElepRecoVsEv0_cov_sigma2","",19,mubins,100,0,16);
  TH2D *m_hElepRecoVsEv3_cov_sigma2 = new TH2D("m_hElepRecoVsEv3_cov_sigma2","",19,mubins,100,0,16);
  TH2D *e_hElepRecoVsEv0_cov_sigma2 = new TH2D("e_hElepRecoVsEv0_cov_sigma2","",19,mubins,100,0,16);
  TH2D *e_hElepRecoVsEv3_cov_sigma2 = new TH2D("e_hElepRecoVsEv3_cov_sigma2","",19,mubins,100,0,16);

  TH2D *m_hElepRecoVsEv0_cov_sigma4 = new TH2D("m_hElepRecoVsEv0_cov_sigma4","",19,mubins,100,0,16);
  TH2D *m_hElepRecoVsEv3_cov_sigma4 = new TH2D("m_hElepRecoVsEv3_cov_sigma4","",19,mubins,100,0,16);
  TH2D *e_hElepRecoVsEv0_cov_sigma4 = new TH2D("e_hElepRecoVsEv0_cov_sigma4","",19,mubins,100,0,16);
  TH2D *e_hElepRecoVsEv3_cov_sigma4 = new TH2D("e_hElepRecoVsEv3_cov_sigma4","",19,mubins,100,0,16);

  TH2D *m_hElepRecoVsEv0_cov_sigma5 = new TH2D("m_hElepRecoVsEv0_cov_sigma5","",19,mubins,100,0,16);
  TH2D *m_hElepRecoVsEv3_cov_sigma5 = new TH2D("m_hElepRecoVsEv3_cov_sigma5","",19,mubins,100,0,16);
  TH2D *e_hElepRecoVsEv0_cov_sigma5 = new TH2D("e_hElepRecoVsEv0_cov_sigma5","",19,mubins,100,0,16);
  TH2D *e_hElepRecoVsEv3_cov_sigma5 = new TH2D("e_hElepRecoVsEv3_cov_sigma5","",19,mubins,100,0,16);

  TH1D *hEv0_nc        = new TH1D("hEv0_nc","",100,0,15);
  TH1D *hEv3_nc        = new TH1D("hEv3_nc","",100,0,15);
  TH1D *true_hEv3_nc   = new TH1D("true_hEv3_nc","",100,0,15);
  TH1D *hEv0           = new TH1D("hEv0","",100,0,15);
  TH1D *hEv3           = new TH1D("hEv3","",100,0,15);
  TH1D *true_hEv3      = new TH1D("true_hEv3","",100,0,15);

  TH1D *hEv0_nc_z      = new TH1D("hEv0_nc_z","",100,0,5);
  TH1D *hEv3_nc_z      = new TH1D("hEv3_nc_z","",100,0,5);
  TH1D *true_hEv3_nc_z = new TH1D("true_hEv3_nc_z","",100,0,5);
  TH1D *hEv0_z         = new TH1D("hEv0_z","",100,0,5);
  TH1D *hEv3_z         = new TH1D("hEv3_z","",100,0,5);
  TH1D *true_hEv3_z    = new TH1D("true_hEv3_z","",100,0,5);

/*
  TH1D *m_hElep0        = new TH1D("m_hElep0","",100,0,16);
  TH1D *m_hElep0_sigma2 = new TH1D("m_hElep0_sigma2","",100,0,16);
  TH1D *m_hElep0_sigma4 = new TH1D("m_hElep0_sigma4","",100,0,16);
  TH1D *m_hElep0_z        = new TH1D("m_hElep0_z","",100,0,8);
  TH1D *m_hElep0_sigma2_z = new TH1D("m_hElep0_sigma2_z","",100,0,8);
  TH1D *m_hElep0_sigma4_z = new TH1D("m_hElep0_sigma4_z","",100,0,8);

  TH1D *m_hElep3        = new TH1D("m_hElep3","",100,0,16);
  TH1D *m_hElep3_sigma2 = new TH1D("m_hElep3_sigma2","",100,0,16);
  TH1D *m_hElep3_sigma4 = new TH1D("m_hElep3_sigma4","",100,0,16);
  TH1D *m_hElep3_z        = new TH1D("m_hElep3_z","",100,0,8);
  TH1D *m_hElep3_sigma2_z = new TH1D("m_hElep3_sigma2_z","",100,0,8);
  TH1D *m_hElep3_sigma4_z = new TH1D("m_hElep3_sigma4_z","",100,0,8);

  TH1D *e_hElep0        = new TH1D("e_hElep0","",100,0,16);
  TH1D *e_hElep0_sigma2 = new TH1D("e_hElep0_sigma2","",100,0,16);
  TH1D *e_hElep0_sigma4 = new TH1D("e_hElep0_sigma4","",100,0,16);

  TH1D *e_hElep0_z        = new TH1D("e_hElep0_z","",100,0,8);
  TH1D *e_hElep0_sigma2_z = new TH1D("e_hElep0_sigma2_z","",100,0,8);
  TH1D *e_hElep0_sigma4_z = new TH1D("e_hElep0_sigma4_z","",100,0,8);

  TH1D *e_hElep3        = new TH1D("e_hElep3","",100,0,16);
  TH1D *e_hElep3_sigma2 = new TH1D("e_hElep3_sigma2","",100,0,16);
  TH1D *e_hElep3_sigma4 = new TH1D("e_hElep3_sigma4","",100,0,16);
  TH1D *e_hElep3_z        = new TH1D("e_hElep3_z","",100,0,8);
  TH1D *e_hElep3_sigma2_z = new TH1D("e_hElep3_sigma2_z","",100,0,8);
  TH1D *e_hElep3_sigma4_z = new TH1D("e_hElep3_sigma4_z","",100,0,8);
*/

  TH2D *m_hEvRecoVsEv0    = new TH2D("m_hEvRecoVsEv0","",   nbinsX,xEdges,100,0,20);
  TH2D *m_hEvRecoVsEv3    = new TH2D("m_hEvRecoVsEv3","",   nbinsX,xEdges,100,0,20);
  TH2D *nc_m_hEvRecoVsEv0 = new TH2D("nc_m_hEvRecoVsEv0","",nbinsX,xEdges,100,0,20);
  TH2D *nc_m_hEvRecoVsEv3 = new TH2D("nc_m_hEvRecoVsEv3","",nbinsX,xEdges,100,0,20);
  TH2D *e_hEvRecoVsEv0    = new TH2D("e_hEvRecoVsEv0","",   nbinsX,xEdges,100,0,20);
  TH2D *e_hEvRecoVsEv3    = new TH2D("e_hEvRecoVsEv3","",   nbinsX,xEdges,100,0,20);

  TH2D *m_hEvRecoVsEv0_cov = new TH2D("m_hEvRecoVsEv0_cov","",19,mubins,100,0,20);
  TH2D *m_hEvRecoVsEv3_cov = new TH2D("m_hEvRecoVsEv3_cov","",19,mubins,100,0,20);
  TH2D *e_hEvRecoVsEv0_cov = new TH2D("e_hEvRecoVsEv0_cov","",7, ebins,100,0,20);
  TH2D *e_hEvRecoVsEv3_cov = new TH2D("e_hEvRecoVsEv3_cov","",7, ebins,100,0,20);

  TH2D *m_hEvRecoVsEv0_cov_sigma1 = new TH2D("m_hEvRecoVsEv0_cov_sigma1","",19,mubins,100,0,16);
  TH2D *m_hEvRecoVsEv3_cov_sigma1 = new TH2D("m_hEvRecoVsEv3_cov_sigma1","",19,mubins,100,0,16);
  TH2D *e_hEvRecoVsEv0_cov_sigma1 = new TH2D("e_hEvRecoVsEv0_cov_sigma1","",19,mubins,100,0,16);
  TH2D *e_hEvRecoVsEv3_cov_sigma1 = new TH2D("e_hEvRecoVsEv3_cov_sigma1","",19,mubins,100,0,16);

  TH2D *m_hEvRecoVsEv0_cov_sigma2 = new TH2D("m_hEvRecoVsEv0_cov_sigma2","",19,mubins,100,0,16);
  TH2D *m_hEvRecoVsEv3_cov_sigma2 = new TH2D("m_hEvRecoVsEv3_cov_sigma2","",19,mubins,100,0,16);
  TH2D *e_hEvRecoVsEv0_cov_sigma2 = new TH2D("e_hEvRecoVsEv0_cov_sigma2","",19,mubins,100,0,16);
  TH2D *e_hEvRecoVsEv3_cov_sigma2 = new TH2D("e_hEvRecoVsEv3_cov_sigma2","",19,mubins,100,0,16);

  TH2D *m_hEvRecoVsEv0_cov_sigma4 = new TH2D("m_hEvRecoVsEv0_cov_sigma4","",19,mubins,100,0,16);
  TH2D *m_hEvRecoVsEv3_cov_sigma4 = new TH2D("m_hEvRecoVsEv3_cov_sigma4","",19,mubins,100,0,16);
  TH2D *e_hEvRecoVsEv0_cov_sigma4 = new TH2D("e_hEvRecoVsEv0_cov_sigma4","",19,mubins,100,0,16);
  TH2D *e_hEvRecoVsEv3_cov_sigma4 = new TH2D("e_hEvRecoVsEv3_cov_sigma4","",19,mubins,100,0,16);

  TH2D *m_hEvRecoVsEv0_cov_sigma5 = new TH2D("m_hEvRecoVsEv0_cov_sigma5","",19,mubins,100,0,16);
  TH2D *m_hEvRecoVsEv3_cov_sigma5 = new TH2D("m_hEvRecoVsEv3_cov_sigma5","",19,mubins,100,0,16);
  TH2D *e_hEvRecoVsEv0_cov_sigma5 = new TH2D("e_hEvRecoVsEv0_cov_sigma5","",19,mubins,100,0,16);
  TH2D *e_hEvRecoVsEv3_cov_sigma5 = new TH2D("e_hEvRecoVsEv3_cov_sigma5","",19,mubins,100,0,16);

  TH2D *e_ThetaVsElep0 = new TH2D("e_ThetaVsElep0","",400,0.5,60,20,0,1200);
  TH2D *e_ThetaVsElep3 = new TH2D("e_ThetaVsElep3","",400,0.5,60,20,0,1200);
  TH2D *m_ThetaVsElep0 = new TH2D("m_ThetaVsElep0","",400,0.5,60,20,0,1200);
  TH2D *m_ThetaVsElep3 = new TH2D("m_ThetaVsElep3","",400,0.5,60,20,0,1200);

  // Most of them are weights related to systematic uncertainties
  // information about the true neutrino interaction
  double vtx_x, vtx_y, vtx_z; // the position where the neutrino interaction occurred, in cm
  int nuPDG; // PDG code of the neutrino, numu = 14, nue = 12, antineutrinos are negative
  double Ev; // the energy of the neutrino, in GeV
  int LepPDG; // PDG code of the final-state lepton, mu = 13, e = 11
  double LepE; // Total energy of the final-state lepton; note that true nu is not saved but nu = Ev - LepE
  double LepNuAngle; // angle between lepton and neutrino
  double eP, eN, ePip, ePim, ePi0, eOther; // energy in the final state due to different particles
  double Ev_reco, Elep_reco;

  // information about the "reconstruction". We can talk more about what "reconstruction" means for this
  //double Ev_reco, Elep_reco; // the measured neutrino energy and lepton energy. Reco nu = Ev_reco - Elep_reco
  int reco_numu, reco_nue; // = 1 if the reconstruction thinks it's a muon or electron
  int muon_contained, muon_tracker, muon_ecal; // different ways that the muon can be measured
  double Ehad_veto; // hadronic energy near the edge of the detector, which is a hint that we might not have measured all the energy
  
  //std::list <const char *> namelist = {"wgt_MaCCQE"};
  //const char *name[] = {"wgt_MaCCQE"};
  //std::list <const char *> namelist = {"wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};

  //const char *name[] = {"wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};
  
  std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};

  const char *name[] = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};


  int N_wgt = namelist.size();
  std::cout << "N_wgt = " << N_wgt << "\n";

    TH1D** hm01 = new TH1D*[N_wgt];
    TH1D** hm31 = new TH1D*[N_wgt];
    TH1D** hm02 = new TH1D*[N_wgt];
    TH1D** hm32 = new TH1D*[N_wgt];
    TH1D** hm03 = new TH1D*[N_wgt];
    TH1D** hm33 = new TH1D*[N_wgt];
    TH1D** hm04 = new TH1D*[N_wgt];
    TH1D** hm34 = new TH1D*[N_wgt];
    TH1D** hm05 = new TH1D*[N_wgt];
    TH1D** hm35 = new TH1D*[N_wgt];

    TH1D** he01 = new TH1D*[N_wgt];
    TH1D** he31 = new TH1D*[N_wgt];
    TH1D** he02 = new TH1D*[N_wgt];
    TH1D** he32 = new TH1D*[N_wgt];
    TH1D** he03 = new TH1D*[N_wgt];
    TH1D** he33 = new TH1D*[N_wgt];
    TH1D** he04 = new TH1D*[N_wgt];
    TH1D** he34 = new TH1D*[N_wgt];
    TH1D** he05 = new TH1D*[N_wgt];
    TH1D** he35 = new TH1D*[N_wgt];

    TH2D** tgt_m02    = new TH2D*[N_wgt];  
    TH2D** tgt_m32    = new TH2D*[N_wgt];  
    TH2D** tgt_m03    = new TH2D*[N_wgt];  
    TH2D** tgt_m33    = new TH2D*[N_wgt];  
    TH2D** tgt_m04    = new TH2D*[N_wgt];  
    TH2D** tgt_m34    = new TH2D*[N_wgt];  

    TH2D** tgt_nc_m02 = new TH2D*[N_wgt];  
    TH2D** tgt_nc_m32 = new TH2D*[N_wgt];  
    TH2D** tgt_nc_m03 = new TH2D*[N_wgt];  
    TH2D** tgt_nc_m33 = new TH2D*[N_wgt];  
    TH2D** tgt_nc_m04 = new TH2D*[N_wgt];  
    TH2D** tgt_nc_m34 = new TH2D*[N_wgt];  

    TH2D** tgt_e02    = new TH2D*[N_wgt];  
    TH2D** tgt_e32    = new TH2D*[N_wgt];  
    TH2D** tgt_e03    = new TH2D*[N_wgt];  
    TH2D** tgt_e33    = new TH2D*[N_wgt];  
    TH2D** tgt_e04    = new TH2D*[N_wgt];  
    TH2D** tgt_e34    = new TH2D*[N_wgt];  

    TH1D** Ev_hm01 = new TH1D*[N_wgt];
    TH1D** Ev_hm31 = new TH1D*[N_wgt];
    TH1D** Ev_hm02 = new TH1D*[N_wgt];
    TH1D** Ev_hm32 = new TH1D*[N_wgt];
    TH1D** Ev_hm03 = new TH1D*[N_wgt];
    TH1D** Ev_hm33 = new TH1D*[N_wgt];
    TH1D** Ev_hm04 = new TH1D*[N_wgt];
    TH1D** Ev_hm34 = new TH1D*[N_wgt];
    TH1D** Ev_hm05 = new TH1D*[N_wgt];
    TH1D** Ev_hm35 = new TH1D*[N_wgt];

    TH1D** Ev_he01 = new TH1D*[N_wgt];
    TH1D** Ev_he31 = new TH1D*[N_wgt];
    TH1D** Ev_he02 = new TH1D*[N_wgt];
    TH1D** Ev_he32 = new TH1D*[N_wgt];
    TH1D** Ev_he03 = new TH1D*[N_wgt];
    TH1D** Ev_he33 = new TH1D*[N_wgt];
    TH1D** Ev_he04 = new TH1D*[N_wgt];
    TH1D** Ev_he34 = new TH1D*[N_wgt];
    TH1D** Ev_he05 = new TH1D*[N_wgt];
    TH1D** Ev_he35 = new TH1D*[N_wgt];

    TH2D** Ev_tgt_m02    = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_m32    = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_m03    = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_m33    = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_m04    = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_m34    = new TH2D*[N_wgt];  

    TH2D** Ev_tgt_nc_m02 = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_nc_m32 = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_nc_m03 = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_nc_m33 = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_nc_m04 = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_nc_m34 = new TH2D*[N_wgt];  

    TH2D** Ev_tgt_e02    = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_e32    = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_e03    = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_e33    = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_e04    = new TH2D*[N_wgt];  
    TH2D** Ev_tgt_e34    = new TH2D*[N_wgt];  
  
  //std::map <const char*, TH2D*> hmap_m02, hmap_m32, hmap_e02, hmap_e32, hmap_m04, hmap_m34, hmap_e04, hmap_e34;
  for(int i=0; i<N_wgt; i++){
    hm01[i] = new TH1D(Form("%s_m_hElep0_sigma1",name[i]),"",100,0,16);
    hm31[i] = new TH1D(Form("%s_m_hElep3_sigma1",name[i]),"",100,0,16);
    hm02[i] = new TH1D(Form("%s_m_hElep0_sigma2",name[i]),"",100,0,16);
    hm32[i] = new TH1D(Form("%s_m_hElep3_sigma2",name[i]),"",100,0,16);
    hm03[i] = new TH1D(Form("%s_m_hElep0_sigma3",name[i]),"",100,0,16);
    hm33[i] = new TH1D(Form("%s_m_hElep3_sigma3",name[i]),"",100,0,16);
    hm04[i] = new TH1D(Form("%s_m_hElep0_sigma4",name[i]),"",100,0,16);
    hm34[i] = new TH1D(Form("%s_m_hElep3_sigma4",name[i]),"",100,0,16);
    hm05[i] = new TH1D(Form("%s_m_hElep0_sigma5",name[i]),"",100,0,16);
    hm35[i] = new TH1D(Form("%s_m_hElep3_sigma5",name[i]),"",100,0,16);

    he01[i] = new TH1D(Form("%s_e_hElep0_sigma1",name[i]),"",100,0,16);
    he31[i] = new TH1D(Form("%s_e_hElep3_sigma1",name[i]),"",100,0,16);
    he02[i] = new TH1D(Form("%s_e_hElep0_sigma2",name[i]),"",100,0,16);
    he32[i] = new TH1D(Form("%s_e_hElep3_sigma2",name[i]),"",100,0,16);
    he03[i] = new TH1D(Form("%s_e_hElep0_sigma3",name[i]),"",100,0,16);
    he33[i] = new TH1D(Form("%s_e_hElep3_sigma3",name[i]),"",100,0,16);
    he04[i] = new TH1D(Form("%s_e_hElep0_sigma4",name[i]),"",100,0,16);
    he34[i] = new TH1D(Form("%s_e_hElep3_sigma4",name[i]),"",100,0,16);
    he05[i] = new TH1D(Form("%s_e_hElep0_sigma5",name[i]),"",100,0,16);
    he35[i] = new TH1D(Form("%s_e_hElep3_sigma5",name[i]),"",100,0,16);

    tgt_m02[i]  = new TH2D(Form("%s_m_hElepRecoVsEv0_sigma2",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_m32[i]  = new TH2D(Form("%s_m_hElepRecoVsEv3_sigma2",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_m03[i]  = new TH2D(Form("%s_m_hElepRecoVsEv0_sigma3",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_m33[i]  = new TH2D(Form("%s_m_hElepRecoVsEv3_sigma3",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_m04[i]  = new TH2D(Form("%s_m_hElepRecoVsEv0_sigma4",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_m34[i]  = new TH2D(Form("%s_m_hElepRecoVsEv3_sigma4",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);

    tgt_nc_m02[i]  = new TH2D(Form("%s_nc_m_hElepRecoVsEv0_sigma2",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_nc_m32[i]  = new TH2D(Form("%s_nc_m_hElepRecoVsEv3_sigma2",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_nc_m03[i]  = new TH2D(Form("%s_nc_m_hElepRecoVsEv0_sigma3",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_nc_m33[i]  = new TH2D(Form("%s_nc_m_hElepRecoVsEv3_sigma3",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_nc_m04[i]  = new TH2D(Form("%s_nc_m_hElepRecoVsEv0_sigma4",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_nc_m34[i]  = new TH2D(Form("%s_nc_m_hElepRecoVsEv3_sigma4",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);

    tgt_e02[i]  = new TH2D(Form("%s_e_hElepRecoVsEv0_sigma2",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_e32[i]  = new TH2D(Form("%s_e_hElepRecoVsEv3_sigma2",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_e03[i]  = new TH2D(Form("%s_e_hElepRecoVsEv0_sigma3",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_e33[i]  = new TH2D(Form("%s_e_hElepRecoVsEv3_sigma3",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_e04[i]  = new TH2D(Form("%s_e_hElepRecoVsEv0_sigma4",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    tgt_e34[i]  = new TH2D(Form("%s_e_hElepRecoVsEv3_sigma4",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);


    Ev_hm01[i] = new TH1D(Form("%s_m_hEv0_sigma1",name[i]),"",100,0,15);
    Ev_hm31[i] = new TH1D(Form("%s_m_hEv3_sigma1",name[i]),"",100,0,15);
    Ev_hm02[i] = new TH1D(Form("%s_m_hEv0_sigma2",name[i]),"",100,0,15);
    Ev_hm32[i] = new TH1D(Form("%s_m_hEv3_sigma2",name[i]),"",100,0,15);
    Ev_hm03[i] = new TH1D(Form("%s_m_hEv0_sigma3",name[i]),"",100,0,15);
    Ev_hm33[i] = new TH1D(Form("%s_m_hEv3_sigma3",name[i]),"",100,0,15);
    Ev_hm04[i] = new TH1D(Form("%s_m_hEv0_sigma4",name[i]),"",100,0,15);
    Ev_hm34[i] = new TH1D(Form("%s_m_hEv3_sigma4",name[i]),"",100,0,15);
    Ev_hm05[i] = new TH1D(Form("%s_m_hEv0_sigma5",name[i]),"",100,0,15);
    Ev_hm35[i] = new TH1D(Form("%s_m_hEv3_sigma5",name[i]),"",100,0,15);

    Ev_he01[i] = new TH1D(Form("%s_e_hEv0_sigma1",name[i]),"",100,0,15);
    Ev_he31[i] = new TH1D(Form("%s_e_hEv3_sigma1",name[i]),"",100,0,15);
    Ev_he02[i] = new TH1D(Form("%s_e_hEv0_sigma2",name[i]),"",100,0,15);
    Ev_he32[i] = new TH1D(Form("%s_e_hEv3_sigma2",name[i]),"",100,0,15);
    Ev_he03[i] = new TH1D(Form("%s_e_hEv0_sigma3",name[i]),"",100,0,15);
    Ev_he33[i] = new TH1D(Form("%s_e_hEv3_sigma3",name[i]),"",100,0,15);
    Ev_he04[i] = new TH1D(Form("%s_e_hEv0_sigma4",name[i]),"",100,0,15);
    Ev_he34[i] = new TH1D(Form("%s_e_hEv3_sigma4",name[i]),"",100,0,15);
    Ev_he05[i] = new TH1D(Form("%s_e_hEv0_sigma5",name[i]),"",100,0,15);
    Ev_he35[i] = new TH1D(Form("%s_e_hEv3_sigma5",name[i]),"",100,0,15);

    Ev_tgt_m02[i]  = new TH2D(Form("%s_m_hEvRecoVsEv0_sigma2",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_m32[i]  = new TH2D(Form("%s_m_hEvRecoVsEv3_sigma2",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_m03[i]  = new TH2D(Form("%s_m_hEvRecoVsEv0_sigma3",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_m33[i]  = new TH2D(Form("%s_m_hEvRecoVsEv3_sigma3",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_m04[i]  = new TH2D(Form("%s_m_hEvRecoVsEv0_sigma4",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_m34[i]  = new TH2D(Form("%s_m_hEvRecoVsEv3_sigma4",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);

    Ev_tgt_nc_m02[i]  = new TH2D(Form("%s_nc_m_hEvRecoVsEv0_sigma2",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_nc_m32[i]  = new TH2D(Form("%s_nc_m_hEvRecoVsEv3_sigma2",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_nc_m03[i]  = new TH2D(Form("%s_nc_m_hEvRecoVsEv0_sigma3",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_nc_m33[i]  = new TH2D(Form("%s_nc_m_hEvRecoVsEv3_sigma3",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_nc_m04[i]  = new TH2D(Form("%s_nc_m_hEvRecoVsEv0_sigma4",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_nc_m34[i]  = new TH2D(Form("%s_nc_m_hEvRecoVsEv3_sigma4",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);

    Ev_tgt_e02[i]  = new TH2D(Form("%s_e_hEvRecoVsEv0_sigma2",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_e32[i]  = new TH2D(Form("%s_e_hEvRecoVsEv3_sigma2",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_e03[i]  = new TH2D(Form("%s_e_hEvRecoVsEv0_sigma3",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_e33[i]  = new TH2D(Form("%s_e_hEvRecoVsEv3_sigma3",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_e04[i]  = new TH2D(Form("%s_e_hEvRecoVsEv0_sigma4",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);
    Ev_tgt_e34[i]  = new TH2D(Form("%s_e_hEvRecoVsEv3_sigma4",name[i]),"",nbinsX,xEdges,nbinsY,yEdges);

  }

  double wgt[N_wgt][6];
  //double wgt[5];

  tree->SetBranchAddress( "vtx_x", &vtx_x );
  tree->SetBranchAddress( "vtx_y", &vtx_y );
  tree->SetBranchAddress( "vtx_z", &vtx_z );
  tree->SetBranchAddress( "Ev", &Ev );
  tree->SetBranchAddress( "nuPDG", &nuPDG );
  tree->SetBranchAddress( "LepPDG", &LepPDG );
  tree->SetBranchAddress( "LepE", &LepE );
  tree->SetBranchAddress( "LepNuAngle", &LepNuAngle );
  tree->SetBranchAddress( "Ev_reco", &Ev_reco );
  tree->SetBranchAddress( "Elep_reco", &Elep_reco );
  tree->SetBranchAddress( "reco_numu", &reco_numu );
  tree->SetBranchAddress( "reco_nue", &reco_nue );
  tree->SetBranchAddress( "muon_contained", &muon_contained );
  tree->SetBranchAddress( "muon_tracker", &muon_tracker );
  tree->SetBranchAddress( "muon_ecal", &muon_ecal );
  tree->SetBranchAddress( "Ehad_veto", &Ehad_veto );
  tree->SetBranchAddress( "eP", &eP );
  tree->SetBranchAddress( "eN", &eN );
  tree->SetBranchAddress( "ePip", &ePip );
  tree->SetBranchAddress( "ePim", &ePim );
  tree->SetBranchAddress( "ePi0", &ePi0 );
  tree->SetBranchAddress( "eOther", &eOther );
  tree->SetBranchAddress( "LepNuAngle", &LepNuAngle );

  double E_max=0;

  for(int i=0; i<N_wgt; i++){
    //std::cout << i+1 << "\t" << name[i] << "\n";
    tree->SetBranchAddress( name[i] , &wgt[i] );
  //tree->SetBranchAddress( name[i] , &wgt );
  }

  double nu, true_nu;
  double LepE_sm, Ev_sm;

  TRandom *r1 = new TRandom(8888);
  TRandom3 *rando = new TRandom3(8888); 
 
  const int N = tree->GetEntries();
  //const int N = 1000;
  //double scalePOT = 1.;
  for( int ii = 0; ii < N; ++ii )
  {
    if( ii % 10000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
    tree->GetEntry(ii);

    // Skip events that occur outside the "Fiducial Volume" which is a region in the middle of the detector
//    // Basically we can't measure neutrinos that interact right next to the edge very well
//    wnu.cxx.swp
    // These numbers are in cm; the detector goes from -357 to +357 in x, -150 to +150 in y, and 0 to 507 in z
    if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z > 350. ) continue;
    //if( Ev < 0.2 ) continue;

    true_nu = eP + eN + ePip + ePim + ePi0 + eOther;
      
    if(LepPDG == 13){

        nu = Ev_reco - Elep_reco;

        hEv0_nc->Fill(Ev);
        hEv0_nc_z->Fill(Ev);

      //templates: muon without cuts
      //if(nu<10.0) { 
        nc_m_hElepRecoVsEv0->Fill(Ev,Elep_reco); 
        nc_m_hEvRecoVsEv0->Fill(Ev,Ev_reco); 
  
        for(int i=0; i<N_wgt; i++) {
          tgt_nc_m02[i]->Fill(Ev,Elep_reco,wgt[i][2]);
          tgt_nc_m03[i]->Fill(Ev,Elep_reco,1.);
          tgt_nc_m04[i]->Fill(Ev,Elep_reco,wgt[i][4]);

          Ev_tgt_nc_m02[i]->Fill(Ev,Ev_reco,wgt[i][2]);
          Ev_tgt_nc_m03[i]->Fill(Ev,Ev_reco,1.);
          Ev_tgt_nc_m04[i]->Fill(Ev,Ev_reco,wgt[i][4]);
        }
      //}   
      
      if(true_nu<0.3) {
        true_hEv3_nc->Fill(Ev);
        true_hEv3_nc_z->Fill(Ev);
      }


      if(nu<0.3) { 
        hEv3_nc->Fill(Ev);
        hEv3_nc_z->Fill(Ev);

        nc_m_hElepRecoVsEv3->Fill(Ev,Elep_reco); 
        nc_m_hEvRecoVsEv3->Fill(Ev,Ev_reco); 

        for(int i=0; i<N_wgt; i++) {
          tgt_nc_m32[i]->Fill(Ev,Elep_reco,wgt[i][2]);
          tgt_nc_m33[i]->Fill(Ev,Elep_reco,1.);
          tgt_nc_m34[i]->Fill(Ev,Elep_reco,wgt[i][4]);

          Ev_tgt_nc_m32[i]->Fill(Ev,Ev_reco,wgt[i][2]);
          Ev_tgt_nc_m33[i]->Fill(Ev,Ev_reco,1.);
          Ev_tgt_nc_m34[i]->Fill(Ev,Ev_reco,wgt[i][4]);
        }
      }

    
      if( reco_numu && (muon_contained || muon_tracker || muon_ecal))
      {
        hEv0->Fill(Ev);
        hEv0_z->Fill(Ev);
        //templates: muon with cuts
        //if(nu<10.0) { 
          for(int i=0; i<N_wgt; i++) {
            hm01[i]->Fill(Elep_reco,wgt[i][1]);
            hm02[i]->Fill(Elep_reco,wgt[i][2]);
            hm03[i]->Fill(Elep_reco,1.);
            hm04[i]->Fill(Elep_reco,wgt[i][4]);
            hm05[i]->Fill(Elep_reco,wgt[i][5]);

            tgt_m02[i]->Fill(Ev,Elep_reco,wgt[i][2]);
            tgt_m03[i]->Fill(Ev,Elep_reco,1.);
            tgt_m04[i]->Fill(Ev,Elep_reco,wgt[i][4]);

            Ev_hm01[i]->Fill(Ev_reco,wgt[i][1]);
            Ev_hm02[i]->Fill(Ev_reco,wgt[i][2]);
            Ev_hm03[i]->Fill(Ev_reco,1.);
            Ev_hm04[i]->Fill(Ev_reco,wgt[i][4]);
            Ev_hm05[i]->Fill(Ev_reco,wgt[i][5]);

            Ev_tgt_m02[i]->Fill(Ev,Ev_reco,wgt[i][2]);
            Ev_tgt_m03[i]->Fill(Ev,Ev_reco,1.);
            Ev_tgt_m04[i]->Fill(Ev,Ev_reco,wgt[i][4]);
          }
                      //m_hElep0_sigma2->Fill(LepE_sm,wgt[2]);
                      //m_hElep0_sigma4->Fill(LepE_sm,wgt[4]);

          m_hElepRecoVsEv0->Fill(Ev,Elep_reco); 
          m_hElepRecoVsEv0_cov->Fill(Ev,Elep_reco);

          m_hEvRecoVsEv0->Fill(Ev,Ev_reco);
          m_hEvRecoVsEv0_cov->Fill(Ev,Ev_reco); 
        //}    

        if(true_nu<0.3)  { 
          true_hEv3->Fill(Ev);
          true_hEv3_z->Fill(Ev);
        }

        if(nu<0.3)  { 

          hEv3->Fill(Ev);
          hEv3_z->Fill(Ev);

          for(int i=0; i<N_wgt; i++) {
            hm31[i]->Fill(Elep_reco,wgt[i][1]);
            hm32[i]->Fill(Elep_reco,wgt[i][2]);
            hm33[i]->Fill(Elep_reco,1.);
            hm34[i]->Fill(Elep_reco,wgt[i][4]);
            hm35[i]->Fill(Elep_reco,wgt[i][5]);

            tgt_m32[i]->Fill(Ev,Elep_reco,wgt[i][2]);
            tgt_m33[i]->Fill(Ev,Elep_reco,1.);
            tgt_m34[i]->Fill(Ev,Elep_reco,wgt[i][4]);

            Ev_hm31[i]->Fill(Ev_reco,wgt[i][1]);
            Ev_hm32[i]->Fill(Ev_reco,wgt[i][2]);
            Ev_hm33[i]->Fill(Ev_reco,1.);
            Ev_hm34[i]->Fill(Ev_reco,wgt[i][4]);
            Ev_hm35[i]->Fill(Ev_reco,wgt[i][5]);

            Ev_tgt_m32[i]->Fill(Ev,Ev_reco,wgt[i][2]);
            Ev_tgt_m33[i]->Fill(Ev,Ev_reco,1.);
            Ev_tgt_m34[i]->Fill(Ev,Ev_reco,wgt[i][4]);
          }

          m_hElepRecoVsEv3->Fill(Ev,Elep_reco); 
          m_hElepRecoVsEv3_cov->Fill(Ev,Elep_reco);

          m_hEvRecoVsEv3->Fill(Ev,Ev_reco); 
	  m_hEvRecoVsEv3_cov->Fill(Ev,Ev_reco); 
        }
      }  
    }

    if(LepPDG == 11){

      LepE_sm = r1->Gaus(LepE,0.05);
      Ev_sm   = r1->Gaus(Ev,0.05);
      nu = Ev_sm - LepE_sm; //reco_nu
      
      if(LepE_sm*LepNuAngle*LepNuAngle*1E3>3.){
      if( reco_nue ){
      
      //templates: eletron
      //if(nu<10.0) {
        for(int i=0; i<N_wgt; i++) {
          he01[i]->Fill(LepE_sm,wgt[i][1]);
          he02[i]->Fill(LepE_sm,wgt[i][2]);
          he03[i]->Fill(LepE_sm,1.);
          he04[i]->Fill(LepE_sm,wgt[i][4]);
          he05[i]->Fill(LepE_sm,wgt[i][5]);
            
          tgt_e02[i]->Fill(Ev,LepE_sm,wgt[i][2]);
          tgt_e03[i]->Fill(Ev,LepE_sm,1.);
          tgt_e04[i]->Fill(Ev,LepE_sm,wgt[i][4]);

          Ev_he01[i]->Fill(Ev_sm,wgt[i][1]);
          Ev_he02[i]->Fill(Ev_sm,wgt[i][2]);
          Ev_he03[i]->Fill(Ev_sm,1.);
          Ev_he04[i]->Fill(Ev_sm,wgt[i][4]);
          Ev_he05[i]->Fill(Ev_sm,wgt[i][5]);
            
          Ev_tgt_e02[i]->Fill(Ev,Ev_sm,wgt[i][2]);
          Ev_tgt_e03[i]->Fill(Ev,Ev_sm,1.);
          Ev_tgt_e04[i]->Fill(Ev,Ev_sm,wgt[i][4]);
        }

        e_hElepRecoVsEv0->Fill(Ev,LepE_sm);
        e_hElepRecoVsEv0_cov->Fill(Ev,LepE_sm);

	e_hEvRecoVsEv0->Fill(Ev,Ev_sm);
	e_hEvRecoVsEv0_cov->Fill(Ev,Ev_sm); 
      //}

      if(nu<0.3) {  
        for(int i=0; i<N_wgt; i++) {
          he31[i]->Fill(LepE_sm,wgt[i][1]);
          he32[i]->Fill(LepE_sm,wgt[i][2]);
          he33[i]->Fill(LepE_sm,1.);
          he34[i]->Fill(LepE_sm,wgt[i][4]);
          he35[i]->Fill(LepE_sm,wgt[i][5]);

          tgt_e32[i]->Fill(Ev,LepE_sm,wgt[i][2]);
          tgt_e33[i]->Fill(Ev,LepE_sm,1.);
          tgt_e34[i]->Fill(Ev,LepE_sm,wgt[i][4]);

          Ev_he31[i]->Fill(Ev_sm,wgt[i][1]);
          Ev_he32[i]->Fill(Ev_sm,wgt[i][2]);
          Ev_he33[i]->Fill(Ev_sm,1.);
          Ev_he34[i]->Fill(Ev_sm,wgt[i][4]);
          Ev_he35[i]->Fill(Ev_sm,wgt[i][5]);

          Ev_tgt_e32[i]->Fill(Ev,Ev_sm,wgt[i][2]);
          Ev_tgt_e33[i]->Fill(Ev,Ev_sm,1.);
          Ev_tgt_e34[i]->Fill(Ev,Ev_sm,wgt[i][4]);
        }

	e_hElepRecoVsEv3->Fill(Ev,LepE_sm);
	e_hElepRecoVsEv3_cov->Fill(Ev,LepE_sm);

	e_hEvRecoVsEv3->Fill(Ev,Ev_sm);
	e_hEvRecoVsEv3_cov->Fill(Ev,Ev_sm); 
      }
      }
      }
    }
  }
  
  e_hElepRecoVsEv0->Scale(scalePOT);    e_hElepRecoVsEv3->Scale(scalePOT);
  m_hElepRecoVsEv0->Scale(scalePOT);    m_hElepRecoVsEv3->Scale(scalePOT);
  nc_m_hElepRecoVsEv0->Scale(scalePOT); nc_m_hElepRecoVsEv3->Scale(scalePOT);

  e_hEvRecoVsEv0->Scale(scalePOT);    e_hEvRecoVsEv3->Scale(scalePOT);
  m_hEvRecoVsEv0->Scale(scalePOT);    m_hEvRecoVsEv3->Scale(scalePOT);
  nc_m_hEvRecoVsEv0->Scale(scalePOT); nc_m_hEvRecoVsEv3->Scale(scalePOT);

  e_hElepRecoVsEv0_cov->Scale(scalePOT); e_hElepRecoVsEv3_cov->Scale(scalePOT);
  m_hElepRecoVsEv0_cov->Scale(scalePOT); m_hElepRecoVsEv3_cov->Scale(scalePOT);

  e_hEvRecoVsEv0_cov->Scale(scalePOT); e_hEvRecoVsEv3_cov->Scale(scalePOT);
  m_hEvRecoVsEv0_cov->Scale(scalePOT); m_hEvRecoVsEv3_cov->Scale(scalePOT);

  hEv0_nc->Scale(scalePOT);
  hEv3_nc->Scale(scalePOT);
  true_hEv3_nc->Scale(scalePOT);

  hEv0_nc_z->Scale(scalePOT);
  hEv3_nc_z->Scale(scalePOT);
  true_hEv3_nc_z->Scale(scalePOT);

  hEv0->Scale(scalePOT);
  hEv3->Scale(scalePOT);
  true_hEv3->Scale(scalePOT);

  hEv0_z->Scale(scalePOT);
  hEv3_z->Scale(scalePOT);
  true_hEv3_z->Scale(scalePOT);

  for(int i=0; i<N_wgt; i++){  
  hm01[i]->Scale(scalePOT);
  hm02[i]->Scale(scalePOT);
  hm03[i]->Scale(scalePOT);
  hm04[i]->Scale(scalePOT);
  hm05[i]->Scale(scalePOT);
  hm31[i]->Scale(scalePOT);
  hm32[i]->Scale(scalePOT);
  hm33[i]->Scale(scalePOT);
  hm34[i]->Scale(scalePOT);
  hm35[i]->Scale(scalePOT);

  he01[i]->Scale(scalePOT);
  he02[i]->Scale(scalePOT);
  he03[i]->Scale(scalePOT);
  he04[i]->Scale(scalePOT);
  he05[i]->Scale(scalePOT);
  he31[i]->Scale(scalePOT);
  he32[i]->Scale(scalePOT);
  he33[i]->Scale(scalePOT);
  he34[i]->Scale(scalePOT);
  he35[i]->Scale(scalePOT);

  tgt_m02[i]->Scale(scalePOT);
  tgt_m03[i]->Scale(scalePOT);
  tgt_m04[i]->Scale(scalePOT);
  tgt_m32[i]->Scale(scalePOT);
  tgt_m33[i]->Scale(scalePOT);
  tgt_m34[i]->Scale(scalePOT);

  tgt_nc_m02[i]->Scale(scalePOT);
  tgt_nc_m03[i]->Scale(scalePOT);
  tgt_nc_m04[i]->Scale(scalePOT);
  tgt_nc_m32[i]->Scale(scalePOT);
  tgt_nc_m33[i]->Scale(scalePOT);
  tgt_nc_m34[i]->Scale(scalePOT);

  tgt_e02[i]->Scale(scalePOT);
  tgt_e03[i]->Scale(scalePOT);
  tgt_e04[i]->Scale(scalePOT);
  tgt_e32[i]->Scale(scalePOT);
  tgt_e33[i]->Scale(scalePOT);
  tgt_e34[i]->Scale(scalePOT);

  Ev_hm01[i]->Scale(scalePOT);
  Ev_hm02[i]->Scale(scalePOT);
  Ev_hm03[i]->Scale(scalePOT);
  Ev_hm04[i]->Scale(scalePOT);
  Ev_hm05[i]->Scale(scalePOT);
  Ev_hm31[i]->Scale(scalePOT);
  Ev_hm32[i]->Scale(scalePOT);
  Ev_hm33[i]->Scale(scalePOT);
  Ev_hm34[i]->Scale(scalePOT);
  Ev_hm35[i]->Scale(scalePOT);

  Ev_he01[i]->Scale(scalePOT);
  Ev_he02[i]->Scale(scalePOT);
  Ev_he03[i]->Scale(scalePOT);
  Ev_he04[i]->Scale(scalePOT);
  Ev_he05[i]->Scale(scalePOT);
  Ev_he31[i]->Scale(scalePOT);
  Ev_he32[i]->Scale(scalePOT);
  Ev_he33[i]->Scale(scalePOT);
  Ev_he34[i]->Scale(scalePOT);
  Ev_he35[i]->Scale(scalePOT);

  Ev_tgt_m02[i]->Scale(scalePOT);
  Ev_tgt_m03[i]->Scale(scalePOT);
  Ev_tgt_m04[i]->Scale(scalePOT);
  Ev_tgt_m32[i]->Scale(scalePOT);
  Ev_tgt_m33[i]->Scale(scalePOT);
  Ev_tgt_m34[i]->Scale(scalePOT);

  Ev_tgt_nc_m02[i]->Scale(scalePOT);
  Ev_tgt_nc_m03[i]->Scale(scalePOT);
  Ev_tgt_nc_m04[i]->Scale(scalePOT);
  Ev_tgt_nc_m32[i]->Scale(scalePOT);
  Ev_tgt_nc_m33[i]->Scale(scalePOT);
  Ev_tgt_nc_m34[i]->Scale(scalePOT);

  Ev_tgt_e02[i]->Scale(scalePOT);
  Ev_tgt_e03[i]->Scale(scalePOT);
  Ev_tgt_e04[i]->Scale(scalePOT);
  Ev_tgt_e32[i]->Scale(scalePOT);
  Ev_tgt_e33[i]->Scale(scalePOT);
  Ev_tgt_e34[i]->Scale(scalePOT);
  }


/*
  gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  m_hElep0->SetTitle("CC numu Elep_reco (nu<10.0GeV) (0sigma)");
  m_hElep0->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep0->GetYaxis()->SetTitle("counts");
  m_hElep0_sigma2->SetTitle("CC numu Elep_reco (nu<10.0GeV) (-1sigma)");
  m_hElep0_sigma2->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep0_sigma2->GetYaxis()->SetTitle("counts");
  m_hElep0_sigma4->SetTitle("CC numu Elep_reco (nu<10.0GeV) (+1sigma)");
  m_hElep0_sigma4->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep0_sigma4->GetYaxis()->SetTitle("counts");
  TCanvas *m_cElep0 = new TCanvas("m_cElep0","",2700,700);
  m_cElep0->Divide(3,1);
  m_cElep0->cd(1);
  m_hElep0->Draw();
  m_cElep0->cd(2);
  m_hElep0_sigma2->Draw();
  m_cElep0->cd(3);
  m_hElep0_sigma4->Draw();
  m_cElep0->SaveAs("CC_m_cElep0.png");

  m_hElep0_z->SetTitle("CC numu Elep_reco (nu<10.0GeV) (0sigma)");
  m_hElep0_z->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep0_z->GetYaxis()->SetTitle("counts");
  m_hElep0_sigma2_z->SetTitle("CC numu Elep_reco (nu<10.0GeV) (-1sigma)");
  m_hElep0_sigma2_z->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep0_sigma2_z->GetYaxis()->SetTitle("counts");
  m_hElep0_sigma4_z->SetTitle("CC numu Elep_reco (nu<10.0GeV) (+1sigma)");
  m_hElep0_sigma4_z->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep0_sigma4_z->GetYaxis()->SetTitle("counts");
  TCanvas *m_cElep0_z = new TCanvas("m_cElep0_z","",2700,700);
  m_cElep0_z->Divide(3,1);
  m_cElep0_z->cd(1);
  m_hElep0_z->Draw();
  m_cElep0_z->cd(2);
  m_hElep0_sigma2_z->Draw();
  m_cElep0_z->cd(3);
  m_hElep0_sigma4_z->Draw();
  m_cElep0_z->SaveAs("CC_m_cElep0_z.png");

  m_hElep3->SetTitle("CC numu Elep_reco (nu<0.3GeV) (0sigma)");
  m_hElep3->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep3->GetYaxis()->SetTitle("counts");
  m_hElep3_sigma2->SetTitle("CC numu Elep_reco (nu<0.3GeV) (-1sigma)");
  m_hElep3_sigma2->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep3_sigma2->GetYaxis()->SetTitle("counts");
  m_hElep3_sigma4->SetTitle("CC numu Elep_reco (nu<0.3GeV) (+1sigma)");
  m_hElep3_sigma4->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep3_sigma4->GetYaxis()->SetTitle("counts");
  TCanvas *m_cElep3 = new TCanvas("m_cElep3","",2700,700);
  m_cElep3->Divide(3,1);
  m_cElep3->cd(1);
  m_hElep3->Draw();
  m_cElep3->cd(2);
  m_hElep3_sigma2->Draw();
  m_cElep3->cd(3);
  m_hElep3_sigma4->Draw();
  m_cElep3->SaveAs("CC_m_cElep3.png");

  m_hElep3_z->SetTitle("CC numu Elep_reco (nu<0.3GeV) (0sigma)");
  m_hElep3_z->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep3_z->GetYaxis()->SetTitle("counts");
  m_hElep3_sigma2_z->SetTitle("CC numu Elep_reco (nu<0.3GeV) (-1sigma)");
  m_hElep3_sigma2_z->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep3_sigma2_z->GetYaxis()->SetTitle("counts");
  m_hElep3_sigma4_z->SetTitle("CC numu Elep_reco (nu<0.3GeV) (+1sigma)");
  m_hElep3_sigma4_z->GetXaxis()->SetTitle("Elep_reco (GeV)");
  m_hElep3_sigma4_z->GetYaxis()->SetTitle("counts");
  TCanvas *m_cElep3_z = new TCanvas("m_cElep3_z","",2700,700);
  m_cElep3_z->Divide(3,1);
  m_cElep3_z->cd(1);
  m_hElep3_z->Draw();
  m_cElep3_z->cd(2);
  m_hElep3_sigma2_z->Draw();
  m_cElep3_z->cd(3);
  m_hElep3_sigma4_z->Draw();
  m_cElep3_z->SaveAs("CC_m_cElep3_z.png");

  e_hElep0->SetTitle("CC nue Elep_reco (nu<10.0GeV) (0sigma)");
  e_hElep0->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep0->GetYaxis()->SetTitle("counts");
  e_hElep0_sigma2->SetTitle("CC nue Elep_reco (nu<10.0GeV) (-1sigma)");
  e_hElep0_sigma2->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep0_sigma2->GetYaxis()->SetTitle("counts");
  e_hElep0_sigma4->SetTitle("CC nue Elep_reco (nu<10.0GeV) (+1sigma)");
  e_hElep0_sigma4->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep0_sigma4->GetYaxis()->SetTitle("counts");
  TCanvas *e_cElep0 = new TCanvas("e_cElep0","",2700,700);
  e_cElep0->Divide(3,1);
  e_cElep0->cd(1);
  e_hElep0->Draw();
  e_cElep0->cd(2);
  e_hElep0_sigma2->Draw();
  e_cElep0->cd(3);
  e_hElep0_sigma4->Draw();
  e_cElep0->SaveAs("CC_e_cElep0.png");

  e_hElep0_z->SetTitle("CC nue Elep_reco (nu<10.0GeV) (0sigma)");
  e_hElep0_z->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep0_z->GetYaxis()->SetTitle("counts");
  e_hElep0_sigma2_z->SetTitle("CC nue Elep_reco (nu<10.0GeV) (-1sigma)");
  e_hElep0_sigma2_z->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep0_sigma2_z->GetYaxis()->SetTitle("counts");
  e_hElep0_sigma4_z->SetTitle("CC nue Elep_reco (nu<10.0GeV) (+1sigma)");
  e_hElep0_sigma4_z->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep0_sigma4_z->GetYaxis()->SetTitle("counts");
  TCanvas *e_cElep0_z = new TCanvas("e_cElep0_z","",2700,700);
  e_cElep0_z->Divide(3,1);
  e_cElep0_z->cd(1);
  e_hElep0_z->Draw();
  e_cElep0_z->cd(2);
  e_hElep0_sigma2_z->Draw();
  e_cElep0_z->cd(3);
  e_hElep0_sigma4_z->Draw();
  e_cElep0_z->SaveAs("CC_e_cElep0_z.png");

  e_hElep3->SetTitle("CC nue Elep_reco (nu<0.3GeV) (0sigma)");
  e_hElep3->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep3->GetYaxis()->SetTitle("counts");
  e_hElep3_sigma2->SetTitle("CC nue Elep_reco (nu<0.3GeV) (-1sigma)");
  e_hElep3_sigma2->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep3_sigma2->GetYaxis()->SetTitle("counts");
  e_hElep3_sigma4->SetTitle("CC nue Elep_reco (nu<0.3GeV) (+1sigma)");
  e_hElep3_sigma4->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep3_sigma4->GetYaxis()->SetTitle("counts");
  TCanvas *e_cElep3 = new TCanvas("e_cElep3","",2700,700);
  e_cElep3->Divide(3,1);
  e_cElep3->cd(1);
  e_hElep3->Draw();
  e_cElep3->cd(2);
  e_hElep3_sigma2->Draw();
  e_cElep3->cd(3);
  e_hElep3_sigma4->Draw();
  e_cElep3->SaveAs("CC_e_cElep3.png");

  e_hElep3_z->SetTitle("CC nue Elep_reco (nu<0.3GeV) (0sigma)");
  e_hElep3_z->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep3_z->GetYaxis()->SetTitle("counts");
  e_hElep3_sigma2_z->SetTitle("CC nue Elep_reco (nu<0.3GeV) (-1sigma)");
  e_hElep3_sigma2_z->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep3_sigma2_z->GetYaxis()->SetTitle("counts");
  e_hElep3_sigma4_z->SetTitle("CC nue Elep_reco (nu<0.3GeV) (+1sigma)");
  e_hElep3_sigma4_z->GetXaxis()->SetTitle("Elep_reco (GeV)");
  e_hElep3_sigma4_z->GetYaxis()->SetTitle("counts");
  TCanvas *e_cElep3_z = new TCanvas("e_cElep3_z","",2700,700);
  e_cElep3_z->Divide(3,1);
  e_cElep3_z->cd(1);
  e_hElep3_z->Draw();
  e_cElep3_z->cd(2);
  e_hElep3_sigma2_z->Draw();
  e_cElep3_z->cd(3);
  e_hElep3_sigma4_z->Draw();
  e_cElep3_z->SaveAs("CC_e_cElep3_z.png");

  m_hElepRecoVsEv0->SetTitle("CC numu Elep_reco Vs true Ev (nu<10.0)");
  m_hElepRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv0->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepRecoVsEv3->SetTitle("CC numu Elep_reco Vs true Ev (nu<0.3)");
  m_hElepRecoVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv3->GetYaxis()->SetTitle("Elep_reco (GeV)");
  nc_m_hElepRecoVsEv0->SetTitle("CC numu Elep_reco Vs true Ev (nu<10.0) without reconstruction cuts");
  nc_m_hElepRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  nc_m_hElepRecoVsEv0->GetYaxis()->SetTitle("Elep_reco (GeV)");
  nc_m_hElepRecoVsEv3->SetTitle("CC numu Elep_reco Vs true Ev (nu<0.3) without reconstruction cuts");
  nc_m_hElepRecoVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  nc_m_hElepRecoVsEv3->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv0->SetTitle("CC nue Elep_reco Vs true Ev (nu<10.0)");
  e_hElepRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv0->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv3->SetTitle("CC nue Elep_reco Vs true Ev (nu<0.3)");
  e_hElepRecoVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv3->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv = new TCanvas("cElepVsEv","",2700,1400);
  cElepVsEv->Divide(3,2);
  cElepVsEv->cd(1);
  m_hElepRecoVsEv0->Draw("colz");
  cElepVsEv->cd(2);
  nc_m_hElepRecoVsEv0->Draw("colz");
  cElepVsEv->cd(3);
  e_hElep3ecoVsEv0->Draw("colz");
  cElepVsEv->cd(4);
  m_hElepRecoVsEv3->Draw("colz");
  cElepVsEv->cd(5);
  nc_m_hElepRecoVsEv3->Draw("colz");
  cElepVsEv->cd(6);
  e_hElepRecoVsEv3->Draw("colz");
  cElepVsEv->SaveAs("CC_Elep_templates.png");

  m_hEvRecoVsEv0->SetTitle("CC numu Ev_reco Vs true Ev (nu<10.0)");
  m_hEvRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv0->GetYaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvRecoVsEv3->SetTitle("CC numu Ev_reco Vs true Ev (nu<0.3)");
  m_hEvRecoVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv3->GetYaxis()->SetTitle("Ev_reco (GeV)");
  nc_m_hEvRecoVsEv0->SetTitle("CC numu Ev_reco Vs true Ev (nu<10.0) without reconstruction cuts");
  nc_m_hEvRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  nc_m_hEvRecoVsEv0->GetYaxis()->SetTitle("Ev_reco (GeV)");
  nc_m_hEvRecoVsEv3->SetTitle("CC numu Ev_reco Vs true Ev (nu<0.3) without reconstruction cuts");
  nc_m_hEvRecoVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  nc_m_hEvRecoVsEv3->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv0->SetTitle("CC nue Ev_reco Vs true Ev (nu<10.0)");
  e_hEvRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv0->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv3->SetTitle("CC nue Ev_reco Vs true Ev (nu<0.3)");
  e_hEvRecoVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv3->GetYaxis()->SetTitle("Ev_reco (GeV)");
  TCanvas *cEvVsEv = new TCanvas("cEvVsEv","",2700,1400);
  cEvVsEv->Divide(3,2);
  cEvVsEv->cd(1);
  m_hEvRecoVsEv0->Draw("colz");
  cEvVsEv->cd(2);
  nc_m_hEvRecoVsEv0->Draw("colz");
  cEvVsEv->cd(3);
  e_hEvRecoVsEv0->Draw("colz");
  cEvVsEv->cd(4);
  m_hEvRecoVsEv3->Draw("colz");
  cEvVsEv->cd(5);
  nc_m_hEvRecoVsEv3->Draw("colz");
  cEvVsEv->cd(6);
  e_hEvRecoVsEv3->Draw("colz");
  cEvVsEv->SaveAs("CC_Ev_templates.png");

  m_hElepRecoVsEv0_cov->SetTitle("CC numu Elep_reco Vs true Ev (nu<10.0)");
  m_hElepRecoVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv0_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepRecoVsEv3_cov->SetTitle("CC numu Elep_reco Vs true Ev (nu<0.3)");
  m_hElepRecoVsEv3_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv3_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv0_cov->SetTitle("CC nue Elep_reco Vs true Ev (nu<10.0)");
  e_hElepRecoVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv0_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv3_cov->SetTitle("CC nue Elep_reco Vs true Ev (nu<0.3)");
  e_hElepRecoVsEv3_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv3_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv_cov = new TCanvas("cElepVsEv_cov","",1800,1400);
  cElepVsEv_cov->Divide(2,2);
  cElepVsEv_cov->cd(1);
  m_hElepRecoVsEv0_cov->Draw("colz");
  cElepVsEv_cov->cd(2);
  e_hElepRecoVsEv0_cov->Draw("colz");
  cElepVsEv_cov->cd(3);
  m_hElepRecoVsEv3_cov->Draw("colz");
  cElepVsEv_cov->cd(4);
  e_hElepRecoVsEv3_cov->Draw("colz");
  cElepVsEv_cov->SaveAs("cov_CC_Elep.png");

  m_hEvRecoVsEv0_cov->SetTitle("CC numu Ev_reco Vs true Ev (nu<10.0)");
  m_hEvRecoVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv0_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  m_hEvRecoVsEv3_cov->SetTitle("CC numu Ev_reco Vs true Ev (nu<0.3)");
  m_hEvRecoVsEv3_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hEvRecoVsEv3_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv0_cov->SetTitle("CC nue Ev_reco Vs true Ev (nu<10.0)");
  e_hEvRecoVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv0_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  e_hEvRecoVsEv3_cov->SetTitle("CC nue Ev_reco Vs true Ev (nu<0.3)");
  e_hEvRecoVsEv3_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hEvRecoVsEv3_cov->GetYaxis()->SetTitle("Ev_reco (GeV)");
  TCanvas *cEvVsEv_cov = new TCanvas("cEvVsEv_cov","",1800,1400);
  cEvVsEv_cov->Divide(2,2);
  cEvVsEv_cov->cd(1);
  m_hEvRecoVsEv0_cov->Draw("colz");
  cEvVsEv_cov->cd(2);
  e_hEvRecoVsEv0_cov->Draw("colz");
  cEvVsEv_cov->cd(3);
  m_hEvRecoVsEv3_cov->Draw("colz");
  cCCQE", "wgt_VecFFCCQEshape"};EvVsEv_cov->cd(4);
  e_hEvRecoVsEv3_cov->Draw("colz");
  cEvVsEv_cov->SaveAs("cov_CC_Ev.png");

  m_hElepRecoVsEv0_cov_sigma2->SetTitle("CC numu Elep_reco Vs true Ev (nu<10.0) (-1sigma)");
  m_hElepRecoVsEv0_cov_sigma2->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv0_cov_sigma2->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepRecoVsEv3_cov_sigma2->SetTitle("CC numu Elep_reco Vs true Ev (nu<0.3) (-1sigma)") ;
  m_hElepRecoVsEv3_cov_sigma2->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv3_cov_sigma2->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv0_cov_sigma2->SetTitle("CC nue Elep_reco Vs true Ev (nu<10.0) (-1sigma)");
  e_hElepRecoVsEv0_cov_sigma2->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv0_cov_sigma2->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv3_cov_sigma2->SetTitle("CC nue Elep_reco Vs true Ev (nu<0.3) (-1sigma)");
  e_hElepRecoVsEv3_cov_sigma2->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv3_cov_sigma2->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv_cov_sigma2 = new TCanvas("cElepVsEv_cov_sigma2","",1800,1400);
  cElepVsEv_cov_sigma2->Divide(2,2);
  cElepVsEv_cov_sigma2->cd(1);
  m_hElepRecoVsEv0_cov_sigma2->Draw("colz");
  cElepVsEv_cov_sigma2->cd(2);
  e_hElepRecoVsEv0_cov_sigma2->Draw("colz");
  cElepVsEv_cov_sigma2->cd(3);
  m_hElepRecoVsEv3_cov_sigma2->Draw("colz");
  cElepVsEv_cov_sigma2->cd(4);
  e_hElepRecoVsEv3_cov_sigma2->Draw("colz");
  cElepVsEv_cov_sigma2->SaveAs("cov_sigma2_CC_Elep.png");

  m_hElepRecoVsEv0_cov_sigma4->SetTitle("CC numu Elep_reco Vs true Ev (nu<10.0) (+1sigma)");
  m_hElepRecoVsEv0_cov_sigma4->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv0_cov_sigma4->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepRecoVsEv3_cov_sigma4->SetTitle("CC numu Elep_reco Vs true Ev (nu<0.3) (+1sigma)") ;
  m_hElepRecoVsEv3_cov_sigma4->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv3_cov_sigma4->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv0_cov_sigma4->SetTitle("CC nue Elep_reco Vs true Ev (nu<10.0) (+1sigma)");
  e_hElepRecoVsEv0_cov_sigma4->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv0_cov_sigma4->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv3_cov_sigma4->SetTitle("CC nue Elep_reco Vs true Ev (nu<0.3) (+1sigma)");
  e_hElepRecoVsEv3_cov_sigma4->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv3_cov_sigma4->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv_cov_sigma4 = new TCanvas("cElepVsEv_cov_sigma4","",1800,1400);
  cElepVsEv_cov_sigma4->Divide(2,2);
  cElepVsEv_cov_sigma4->cd(1);
  m_hElepRecoVsEv0_cov_sigma4->Draw("colz");
  cElepVsEv_cov_sigma4->cd(2);
  e_hElepRecoVsEv0_cov_sigma4->Draw("colz");
  cElepVsEv_cov_sigma4->cd(3);
  m_hElepRecoVsEv3_cov_sigma4->Draw("colz");
  cElepVsEv_cov_sigma4->cd(4);
  e_hElepRecoVsEv3_cov_sigma4->Draw("colz");
  cElepVsEv_cov_sigma4->SaveAs("cov_sigma4_CC_Elep.png");
*/

  //std::cout << E_max;

  TFile *out = new TFile("/dune/app/users/qvuong/data/lownu/CC_output.root","RECREATE");

  m_hElepRecoVsEv0->Write();
  m_hElepRecoVsEv3->Write();
  nc_m_hElepRecoVsEv0->Write();
  nc_m_hElepRecoVsEv3->Write();
  e_hElepRecoVsEv0->Write();
  e_hElepRecoVsEv3->Write();
  m_hElepRecoVsEv0_cov->Write();
  m_hElepRecoVsEv3_cov->Write();
  e_hElepRecoVsEv0_cov->Write();
  e_hElepRecoVsEv3_cov->Write();

  hEv0_nc->Write();
  hEv3_nc->Write();
  true_hEv3_nc->Write();

  hEv0_nc_z->Write();
  hEv3_nc_z->Write();
  true_hEv3_nc_z->Write();

  hEv0->Write();
  hEv3->Write();
  true_hEv3->Write();

  hEv0_z->Write();
  hEv3_z->Write();
  true_hEv3_z->Write();

  for(int i=0; i<N_wgt; i++){
  hm01[i]->Write();
  hm02[i]->Write();
  hm03[i]->Write();
  hm04[i]->Write();
  hm05[i]->Write();
  hm31[i]->Write();
  hm32[i]->Write();
  hm33[i]->Write();
  hm34[i]->Write();
  hm35[i]->Write();

  he01[i]->Write();
  he02[i]->Write();
  he03[i]->Write();
  he04[i]->Write();
  he05[i]->Write();
  he31[i]->Write();
  he32[i]->Write();
  he33[i]->Write();
  he34[i]->Write();
  he35[i]->Write();

  tgt_m02[i]->Write();
  tgt_m03[i]->Write();
  tgt_m04[i]->Write();
  tgt_m32[i]->Write();
  tgt_m33[i]->Write();
  tgt_m34[i]->Write();

  tgt_nc_m02[i]->Write();
  tgt_nc_m03[i]->Write();
  tgt_nc_m04[i]->Write();
  tgt_nc_m32[i]->Write();
  tgt_nc_m33[i]->Write();
  tgt_nc_m34[i]->Write();

  tgt_e02[i]->Write();
  tgt_e03[i]->Write();
  tgt_e04[i]->Write();
  tgt_e32[i]->Write();
  tgt_e33[i]->Write();
  tgt_e34[i]->Write();

  Ev_hm01[i]->Write();
  Ev_hm02[i]->Write();
  Ev_hm03[i]->Write();
  Ev_hm04[i]->Write();
  Ev_hm05[i]->Write();
  Ev_hm31[i]->Write();
  Ev_hm32[i]->Write();
  Ev_hm33[i]->Write();
  Ev_hm34[i]->Write();
  Ev_hm35[i]->Write();

  Ev_he01[i]->Write();
  Ev_he02[i]->Write();
  Ev_he03[i]->Write();
  Ev_he04[i]->Write();
  Ev_he05[i]->Write();
  Ev_he31[i]->Write();
  Ev_he32[i]->Write();
  Ev_he33[i]->Write();
  Ev_he34[i]->Write();
  Ev_he35[i]->Write();

  Ev_tgt_m02[i]->Write();
  Ev_tgt_m03[i]->Write();
  Ev_tgt_m04[i]->Write();
  Ev_tgt_m32[i]->Write();
  Ev_tgt_m33[i]->Write();
  Ev_tgt_m34[i]->Write();

  Ev_tgt_nc_m02[i]->Write();
  Ev_tgt_nc_m03[i]->Write();
  Ev_tgt_nc_m04[i]->Write();
  Ev_tgt_nc_m32[i]->Write();
  Ev_tgt_nc_m33[i]->Write();
  Ev_tgt_nc_m34[i]->Write();

  Ev_tgt_e02[i]->Write();
  Ev_tgt_e03[i]->Write();
  Ev_tgt_e04[i]->Write();
  Ev_tgt_e32[i]->Write();
  Ev_tgt_e33[i]->Write();
  Ev_tgt_e34[i]->Write();
  }
  out->Close();

  //gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  return(0);

}



