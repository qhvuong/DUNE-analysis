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
#include <TParameter.h>

void test()
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
  //double scalePOT = 1;

  const int N_nucut = 5;
  double nucut[N_nucut] = {50., 5., 1., 0.5, 0.3};

  const Int_t nbinsX = 430;
  const Int_t nbinsCC = 56;
  double xEdges[nbinsX+1], CCEdges[nbinsCC+1];
  xEdges[0]=CCEdges[0]=0.;
  CCEdges[1] = 0.3;

  //neutrino energy binning
  for(int i=0; i<nbinsX+1; i++)
  {
    if(i<200)                xEdges[i+1] = xEdges[i] + 0.02;
    else if(i>=200 && i<240) xEdges[i+1] = xEdges[i] + 0.1;
    else if(i>=240 && i<400) xEdges[i+1] = xEdges[i] + 0.2;
    else if(i>=400 && i<420) xEdges[i+1] = xEdges[i] + 1.0;
    else                     xEdges[i+1] = xEdges[i] + 4.0;
  }

  //lepton energy binning
  for(int i=1; i<nbinsCC+1; i++)
  {
    if(i<38)               CCEdges[i+1] = CCEdges[i] + 0.1;
    else if(i>=38 && i<43) CCEdges[i+1] = CCEdges[i] + 0.2;
    else if(i>=43 && i<48) CCEdges[i+1] = CCEdges[i] + 0.4;
    else if(i>=48 && i<53) CCEdges[i+1] = CCEdges[i] + 0.8;
    else if(i>=53 && i<55) CCEdges[i+1] = CCEdges[i] + 1.5;
    else                   CCEdges[i+1] = CCEdges[i] + 2.0;
  }

  const double mubins[20] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,8.,12.,16.,20.,40.,100.};
  const double ebins[8] = {0.,2.,4.,6.,8.,10.,20.,100.};

  TH2D** mElepRecoVsEv = new TH2D*[N_nucut];
  TH2D** nc_mElepRecoVsEv = new TH2D*[N_nucut];
  TH2D** eElepRecoVsEv = new TH2D*[N_nucut];

  TH2D** mElepRecoVsEv_cov = new TH2D*[N_nucut];
  TH2D** eElepRecoVsEv_cov = new TH2D*[N_nucut];

  for(int j=0; j<N_nucut; j++){
    mElepRecoVsEv[j] = new TH2D(Form("mElepRecoVsEv%d",j),"",nbinsX,xEdges,nbinsCC,CCEdges);
    nc_mElepRecoVsEv[j] = new TH2D(Form("nc_mElepRecoVsEv%d",j),"",nbinsX,xEdges,nbinsCC,CCEdges);
    eElepRecoVsEv[j] = new TH2D(Form("eElepRecoVsEv%d",j),"",nbinsX,xEdges,nbinsCC,CCEdges);

    mElepRecoVsEv_cov[j] = new TH2D(Form("mElepRecoVsEv%d_cov",j),"",19,mubins,nbinsCC,CCEdges);
    eElepRecoVsEv_cov[j] = new TH2D(Form("eElepRecoVsEv%d_cov",j),"",7, ebins,nbinsCC,CCEdges);
  }

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
  
  
  std::list <const char *> namelist = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};

  const char *name[] = {"wgt_MaCCQE", "wgt_VecFFCCQEshape", "wgt_MaNCEL", "wgt_EtaNCEL", "wgt_MaCCRES", "wgt_MvCCRES", "wgt_MaNCRES", "wgt_MvNCRES", "wgt_RDecBR1gamma", "wgt_RDecBR1eta", "wgt_Theta_Delta2Npi", "wgt_AhtBY", "wgt_BhtBY", "wgt_CV1uBY", "wgt_CV2uBY", "wgt_FormZone", "wgt_MFP_pi", "wgt_FrCEx_pi", "wgt_FrElas_pi", "wgt_FrInel_pi", "wgt_FrAbs_pi", "wgt_FrPiProd_pi", "wgt_MFP_N", "wgt_FrCEx_N", "wgt_FrElas_N", "wgt_FrInel_N", "wgt_FrAbs_N", "wgt_FrPiProd_N", "wgt_CCQEPauliSupViaKF", "wgt_Mnv2p2hGaussEnhancement", "wgt_MKSPP_ReWeight", "wgt_E2p2h_A_nu", "wgt_E2p2h_B_nu", "wgt_E2p2h_A_nubar", "wgt_E2p2h_B_nubar", "wgt_NR_nu_n_CC_2Pi", "wgt_NR_nu_n_CC_3Pi", "wgt_NR_nu_p_CC_2Pi", "wgt_NR_nu_p_CC_3Pi", "wgt_NR_nu_np_CC_1Pi", "wgt_NR_nu_n_NC_1Pi", "wgt_NR_nu_n_NC_2Pi", "wgt_NR_nu_n_NC_3Pi", "wgt_NR_nu_p_NC_1Pi", "wgt_NR_nu_p_NC_2Pi", "wgt_NR_nu_p_NC_3Pi", "wgt_NR_nubar_n_CC_1Pi", "wgt_NR_nubar_n_CC_2Pi", "wgt_NR_nubar_n_CC_3Pi", "wgt_NR_nubar_p_CC_1Pi", "wgt_NR_nubar_p_CC_2Pi", "wgt_NR_nubar_p_CC_3Pi", "wgt_NR_nubar_n_NC_1Pi", "wgt_NR_nubar_n_NC_2Pi", "wgt_NR_nubar_n_NC_3Pi", "wgt_NR_nubar_p_NC_1Pi", "wgt_NR_nubar_p_NC_2Pi", "wgt_NR_nubar_p_NC_3Pi", "wgt_BeRPA_A", "wgt_BeRPA_B", "wgt_BeRPA_D", "wgt_BeRPA_E", "wgt_C12ToAr40_2p2hScaling_nu", "wgt_C12ToAr40_2p2hScaling_nubar", "wgt_nuenuebar_xsec_ratio", "wgt_nuenumu_xsec_ratio", "wgt_SPPLowQ2Suppression", "wgt_FSILikeEAvailSmearing"};


  int N_wgt = namelist.size();
  std::cout << "N_wgt = " << N_wgt << "\n";
  
  TH1D*** hmp = new TH1D**[N_nucut];
  TH1D*** hmo = new TH1D**[N_nucut];
  TH1D*** hmn = new TH1D**[N_nucut];

  TH1D*** hep = new TH1D**[N_nucut];
  TH1D*** heo = new TH1D**[N_nucut];
  TH1D*** hen = new TH1D**[N_nucut];

  for(int j=0; j<N_nucut; j++){
    hmp[j] = new TH1D*[N_wgt];
    hmo[j] = new TH1D*[N_wgt];
    hmn[j] = new TH1D*[N_wgt];

    hep[j] = new TH1D*[N_wgt];
    heo[j] = new TH1D*[N_wgt];
    hen[j] = new TH1D*[N_wgt];

    for(int k=0; k<N_wgt; k++){
      hmp[j][k] = new TH1D(Form("%s_mElep%d_p",name[k],j),"",nbinsCC,CCEdges);
      hmo[j][k] = new TH1D(Form("%s_mElep%d_o",name[k],j),"",nbinsCC,CCEdges);
      hmn[j][k] = new TH1D(Form("%s_mElep%d_n",name[k],j),"",nbinsCC,CCEdges);

      hep[j][k] = new TH1D(Form("%s_eElep%d_p",name[k],j),"",nbinsCC,CCEdges);
      heo[j][k] = new TH1D(Form("%s_eElep%d_o",name[k],j),"",nbinsCC,CCEdges);
      hen[j][k] = new TH1D(Form("%s_eElep%d_n",name[k],j),"",nbinsCC,CCEdges);
    }
  }
  

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


  double wgt[N_wgt][7];
  for(int i=0; i<N_wgt; i++){
    tree->SetBranchAddress( name[i] , &wgt[i] );
  }

  double nu, true_nu;
  double LepE_sm, Ev_sm;

  TRandom3 *rando = new TRandom3(12345);
 
  const int N = tree->GetEntries();
  //const int N = 1000;
  //double scalePOT = 1.;
  for( int ii = 0; ii < N; ++ii ) {
    if( ii % 10000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
    tree->GetEntry(ii);

    // Skip events that occur outside the "Fiducial Volume" which is a region in the middle of the detector
    // Basically we can't measure neutrinos that interact right next to the edge very well wnu.cxx.swp
    // These numbers are in cm; the detector goes from -357 to +357 in x, -150 to +150 in y, and 0 to 507 in z
    if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z > 350. ) continue;

    true_nu = eP + eN + ePip + ePim + ePi0 + eOther;
    LepE_sm = rando->Gaus(LepE,0.05);
    Ev_sm   = rando->Gaus(Ev,0.05);
    nu = Ev_sm - LepE_sm; //reco_nu
      
    if(LepPDG == 13){
      for(int j=0; j<N_nucut; j++){
        if(nu<nucut[j]){
          nc_mElepRecoVsEv[j]->Fill(Ev, Elep_reco);
        }
      }
    }

    if( reco_numu && (muon_contained || muon_tracker || muon_ecal ))
    {
      for(int j=0; j<N_nucut; j++){
        if(nu<nucut[j]){
          for(int k=0; k<N_wgt; k++){
            hmn[j][k]->Fill(Elep_reco, wgt[k][2]);
            hmo[j][k]->Fill(Elep_reco, wgt[k][3]);
            hmp[j][k]->Fill(Elep_reco, wgt[k][4]);
          }

          mElepRecoVsEv[j]->Fill(Ev, Elep_reco);
          mElepRecoVsEv_cov[j]->Fill(Ev, Elep_reco);
        }
      }
    }     

    //if(LepPDG == 11){
    if(reco_nue){
      if(LepE_sm*LepNuAngle*LepNuAngle*1E3>3.){
        for(int j=0; j<N_nucut; j++){
          if(nu<nucut[j]){
            for(int k=0; k<N_wgt; k++){
              hen[j][k]->Fill(LepE_sm, wgt[k][2]);
              heo[j][k]->Fill(LepE_sm, wgt[k][3]);
              hep[j][k]->Fill(LepE_sm, wgt[k][4]);
            }

            eElepRecoVsEv[j]->Fill(Ev, LepE_sm);
            eElepRecoVsEv_cov[j]->Fill(Ev, LepE_sm);
          }
        }
    
      }
    }
  }
  
  for(int j=0; j<N_nucut; j++){
    nc_mElepRecoVsEv[j]->Scale(scalePOT);
    mElepRecoVsEv[j]->Scale(scalePOT);
    eElepRecoVsEv[j]->Scale(scalePOT);

    mElepRecoVsEv_cov[j]->Scale(scalePOT);
    eElepRecoVsEv_cov[j]->Scale(scalePOT);

    for(int k=0; k<N_wgt; k++){
      hmn[j][k]->Scale(scalePOT);
      hmo[j][k]->Scale(scalePOT);
      hmp[j][k]->Scale(scalePOT);

      hen[j][k]->Scale(scalePOT);
      heo[j][k]->Scale(scalePOT);
      hep[j][k]->Scale(scalePOT);
    }
  }



  TFile *out = new TFile("/exp/dune/app/users/qvuong/data/lownu/input_dfiles/CCtest_output.root","RECREATE");

  for(int j=0; j<N_nucut; j++){
    nc_mElepRecoVsEv[j]->Write();
    mElepRecoVsEv[j]->Write();
    eElepRecoVsEv[j]->Write();

    mElepRecoVsEv_cov[j]->Write();
    eElepRecoVsEv_cov[j]->Write();

    for(int k=0; k<N_wgt; k++){
      hmn[j][k]->Write();
      hmo[j][k]->Write();
      hmp[j][k]->Write();

      hen[j][k]->Write();
      heo[j][k]->Write();
      hep[j][k]->Write();
    }
  }

  TParameter<double> totalPOT("total_pot", total_pot);
  totalPOT.Write();

  //out->ls();

  out->Close();

}



