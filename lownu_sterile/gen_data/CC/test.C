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
  
  for(int i = 0; i<4; i++){
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


  const Int_t nbinsX = 430; 
  const Int_t nbinsCC = 56;
  double xEdges[nbinsX+1], CCEdges[nbinsCC+1];
  xEdges[0]=CCEdges[0]=0.;
  CCEdges[1] = 0.3;

  for(int i=0; i<nbinsX+1; i++)
  {
    if(i<200)                xEdges[i+1] = xEdges[i] + 0.02;
    else if(i>=200 && i<240) xEdges[i+1] = xEdges[i] + 0.1;
    else if(i>=240 && i<400) xEdges[i+1] = xEdges[i] + 0.2;
    else if(i>=400 && i<420) xEdges[i+1] = xEdges[i] + 1.0;
    else                     xEdges[i+1] = xEdges[i] + 4.0;
  }

  for(int i=1; i<nbinsCC+1; i++)
  {
    if(i<38)               CCEdges[i+1] = CCEdges[i] + 0.1;
    else if(i>=38 && i<43) CCEdges[i+1] = CCEdges[i] + 0.2;
    else if(i>=43 && i<48) CCEdges[i+1] = CCEdges[i] + 0.4;
    else if(i>=48 && i<53) CCEdges[i+1] = CCEdges[i] + 0.8;
    else if(i>=53 && i<55) CCEdges[i+1] = CCEdges[i] + 1.5;
    else                   CCEdges[i+1] = CCEdges[i] + 2.0;
  }


/*
  for(int i=1; i<nbinsCC+1; i++)
  {
    if(i<38)               CCEdges[i+1] = CCEdges[i] + 0.1;
    else if(i>=38 && i<43) CCEdges[i+1] = CCEdges[i] + 0.2;
    else if(i>=43 && i<53) CCEdges[i+1] = CCEdges[i] + 0.5;
    else                   CCEdges[i+1] = CCEdges[i] + 1.2;
  }
*/
  TH1D *hm0 = new TH1D("hm0","",nbinsCC,CCEdges);
  TH1D *hm3 = new TH1D("hm3","",nbinsCC,CCEdges);
  TH1D *he0 = new TH1D("he0","",nbinsCC,CCEdges);
  TH1D *he3 = new TH1D("he3","",nbinsCC,CCEdges);

  const double mubins[20] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,8.,12.,16.,20.,40.,100.};
  const double ebins[8] = {0.,2.,4.,6.,8.,10.,20.,100.};


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
    // Basically we can't measure neutrinos that interact right next to the edge very well wnu.cxx.swp
    // These numbers are in cm; the detector goes from -357 to +357 in x, -150 to +150 in y, and 0 to 507 in z
    if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z > 350. ) continue;

    true_nu = eP + eN + ePip + ePim + ePi0 + eOther;
      
    if(LepPDG == 13){

      nu = Ev_reco - Elep_reco;

      if( reco_numu && (muon_contained || muon_tracker || muon_ecal))
      {
        //templates: muon with cuts
        if(nu<10.0) { 
          hm0->Fill(Elep_reco);
        }    


        if(nu<0.3)  { 
          hm3->Fill(Elep_reco);
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
      if(nu<10.0) {
        he0->Fill(LepE_sm);
      }

      if(nu<0.3) {  
        he3->Fill(LepE_sm);
      }
      }
      }
    }
  }


  hm3->Scale(scalePOT);
  he3->Scale(scalePOT);

  for(int i=0; i<hm0->GetNbinsX(); i++){
    std::cout << i << "\t" << hm0->GetBinContent(i+1) << "\t" << hm3->GetBinContent(i+1) << "\n";}

  TCanvas *c = new TCanvas("c","",1200,500);
  c->Divide(2,1);
  c->cd(1);
  hm0->Draw();
  c->cd(2);
  hm3->Draw();
  c->SaveAs("hm.png");
 
  TCanvas *c1 = new TCanvas("c1","",1200,500);
  c1->Divide(2,1);
  c1->cd(1);
  he0->Draw();
  c1->cd(2);
  he3->Draw();
  c1->SaveAs("he.png");

/* 
  e_hElepRecoVsEv0->Scale(scalePOT);    
  e_hElepRecoVsEv3->Scale(scalePOT);

  m_hElepRecoVsEv0->Scale(scalePOT);    
  m_hElepRecoVsEv3->Scale(scalePOT);

  nc_m_hElepRecoVsEv0->Scale(scalePOT);
  nc_m_hElepRecoVsEv3->Scale(scalePOT);

  e_hElepRecoVsEv0_cov->Scale(scalePOT);
  e_hElepRecoVsEv3_cov->Scale(scalePOT);

  m_hElepRecoVsEv0_cov->Scale(scalePOT); 
  m_hElepRecoVsEv3_cov->Scale(scalePOT);

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
  }


  TFile *out = new TFile("/dune/app/users/qvuong/data/lownu/CC_output_50.root","RECREATE");

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
  }

  TParameter<double> totalPOT("total_pot", total_pot);
  totalPOT.Write();

  out->Close();

  return(0);
*/
}



