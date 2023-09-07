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

int eff()
{

  // Load the CAF file = Common Analysis Format, a standard TTree that we use in DUNE
  TChain * tree = new TChain( "cafTree", "cafTree" );
  TChain * meta = new TChain( "meta", "meta" );
  
  for(int i = 0; i<400; i++){
  if(i<10) {
  tree->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_0%d.root",i) );
  meta->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_0%d.root",i) ); // make certain this is the exact same file(s)
  }
  else {
  tree->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%d.root",i) );
  meta->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%d.root",i) ); // make certain this is the exact same file(s)
  }
  std::cout << "File number:" << i << "\n";
  } 

  double total_pot = 0.;
  double pot;
  meta->SetBranchAddress( "pot", &pot );
  const int Nfiles = meta->GetEntries();
  for( int ii = 0; ii < Nfiles; ++ii ) {
    meta->GetEntry(ii);
  }
  double yrPOT = 1.1E21;
  double scalePOT = yrPOT/total_pot;

  TH1D *hEv0 = new TH1D("hEv0","",50,0,30);
  TH1D *hEv1 = new TH1D("hEv1","",50,0,30);
  TH1D *hEv2 = new TH1D("hEv2","",50,0,30);
  TH1D *hEv3 = new TH1D("hEv3","",50,0,30);
  TH1D *true_hEv0 = new TH1D("true_hEv0","",50,0,30);
  TH1D *true_hEv1 = new TH1D("true_hEv1","",50,0,30);
  TH1D *true_hEv2 = new TH1D("true_hEv2","",50,0,30);
  TH1D *true_hEv3 = new TH1D("true_hEv3","",50,0,30);
  TH1D *reco_hEv0 = new TH1D("reco_hEv0","",50,0,30);
  TH1D *reco_hEv1 = new TH1D("reco_hEv1","",50,0,30);
  TH1D *reco_hEv2 = new TH1D("reco_hEv2","",50,0,30);
  TH1D *reco_hEv3 = new TH1D("reco_hEv3","",50,0,30);
  
  TH1D *hv0 = new TH1D("hv0","",50,0,30);
  TH1D *hv1 = new TH1D("hv1","",50,0,30);
  TH1D *hv2 = new TH1D("hv2","",50,0,30);
  TH1D *hv3 = new TH1D("hv3","",50,0,30);
  TH1D *true_hv0 = new TH1D("true_hv0","",50,0,30);
  TH1D *true_hv1 = new TH1D("true_hv1","",50,0,30);
  TH1D *true_hv2 = new TH1D("true_hv2","",50,0,30);
  TH1D *true_hv3 = new TH1D("true_hv3","",50,0,30);

  TH1D *hEm0 = new TH1D("hEm0","",50,0,30);
  TH1D *hEm1 = new TH1D("hEm1","",50,0,30);
  TH1D *hEm2 = new TH1D("hEm2","",50,0,30);
  TH1D *hEm3 = new TH1D("hEm3","",50,0,30);
  TH1D *true_hEm0 = new TH1D("true_hEm0","",50,0,30);
  TH1D *true_hEm1 = new TH1D("true_hEm1","",50,0,30);
  TH1D *true_hEm2 = new TH1D("true_hEm2","",50,0,30);
  TH1D *true_hEm3 = new TH1D("true_hEm3","",50,0,30);
  
  TH1D *hA0 = new TH1D("hA0","",50,0,TMath::Pi());
  TH1D *hA1 = new TH1D("hA1","",50,0,TMath::Pi());
  TH1D *hA2 = new TH1D("hA2","",50,0,TMath::Pi());
  TH1D *hA3 = new TH1D("hA3","",50,0,TMath::Pi());
  TH1D *true_hA0 = new TH1D("true_hA0","",50,0,TMath::Pi());
  TH1D *true_hA1 = new TH1D("true_hA1","",50,0,TMath::Pi());
  TH1D *true_hA2 = new TH1D("true_hA2","",50,0,TMath::Pi());
  TH1D *true_hA3 = new TH1D("true_hA3","",50,0,TMath::Pi());


  TH2D *hEvRecoVsEv0 = new TH2D("hEvRecoVsEv0","",50,0,30,50,0,30);
  TH2D *hEvRecoVsEv3 = new TH2D("hEvRecoVsEv3","",50,0,30,50,0,30);
  TH2D *hvRecoVsv0 = new TH2D("hvRecoVsv0","",50,0,30,50,0,30);
  TH2D *hvRecoVsv3 = new TH2D("hvRecoVsv3","",50,0,30,50,0,30);
  TH1D *hEvRes0 = new TH1D("hEvRes0","",50,-1,1);
  TH1D *hEvRes1 = new TH1D("hEvRes1","",50,-1,1);
  TH1D *hEvRes2 = new TH1D("hEvRes2","",50,-1,1);
  TH1D *hEvRes3 = new TH1D("hEvRes3","",50,-1,1);


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

  double nu, true_nu;
  double LepE_sm, Ev_sm;

  const int N = tree->GetEntries();
  for( int ii = 0; ii < N; ++ii )
  {
    if( ii % 100000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
    tree->GetEntry(ii);

     // Skip events that occur outside the "Fiducial Volume" which is a region in the middle of the detector
    // Basically we can't measure neutrinos that interact right next to the edge very well
    // These numbers are in cm; the detector goes from -357 to +357 in x, -150 to +150 in y, and 0 to 507 in z
    if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z > 350. ) continue;

    true_nu = eP + eN + ePip + ePim + ePi0 + eOther;
      
    if(LepPDG == 13){

      nu = Ev_reco - Elep_reco;

      if(true_nu<30.0) { hEv0->Fill(Ev); hv0->Fill(true_nu); hEm0->Fill(LepE); hA0->Fill(LepNuAngle); }
      if(true_nu<2.0)  { hEv1->Fill(Ev); hv1->Fill(true_nu); hEm1->Fill(LepE); hA1->Fill(LepNuAngle); }
      if(true_nu<0.8)  { hEv2->Fill(Ev); hv2->Fill(true_nu); hEm2->Fill(LepE); hA2->Fill(LepNuAngle); }
      if(true_nu<0.3)  { hEv3->Fill(Ev); hv3->Fill(true_nu); hEm3->Fill(LepE); hA3->Fill(LepNuAngle); }

      if( reco_numu && (muon_contained || muon_tracker || muon_ecal))
      {
        hEvRecoVsEv0->Fill(Ev,Ev_reco); 
        hvRecoVsv0->Fill(true_nu,nu);

        if(true_nu<30.0) { true_hEv0->Fill(Ev); true_hv0->Fill(true_nu); true_hEm0->Fill(LepE); true_hA0->Fill(LepNuAngle); hEvRes0->Fill((Ev_reco-Ev)/Ev); }
        if(true_nu<2.0)  { true_hEv1->Fill(Ev); true_hv1->Fill(true_nu); true_hEm1->Fill(LepE); true_hA1->Fill(LepNuAngle); hEvRes1->Fill((Ev_reco-Ev)/Ev); }
        if(true_nu<0.8)  { true_hEv2->Fill(Ev); true_hv2->Fill(true_nu); true_hEm2->Fill(LepE); true_hA2->Fill(LepNuAngle); hEvRes2->Fill((Ev_reco-Ev)/Ev); }
        if(true_nu<0.3)  { true_hEv3->Fill(Ev); true_hv3->Fill(true_nu); true_hEm3->Fill(LepE); true_hA3->Fill(LepNuAngle); hEvRes3->Fill((Ev_reco-Ev)/Ev); hEvRecoVsEv3->Fill(Ev,Ev_reco); }

        if(nu<30.0) reco_hEv0->Fill(Ev_reco);
        if(nu<2.0)  reco_hEv1->Fill(Ev_reco);
        if(nu<0.8)  reco_hEv2->Fill(Ev_reco);
        if(nu<0.3)  reco_hEv3->Fill(Ev_reco);
      }  
    }

/*
    if(LepPDG == 11){

      LepE_sm = r1->Gaus(LepE,0.05);
      Ev_sm   = r1->Gaus(Ev,0.05);
      nu = Ev_sm - LepE_sm; //reco_nu
      
      if(reco_nue && LepE_sm*LepNuAngle*LepNuAngle*1E3>3.)
      {
        
      }
    }
*/
  }

  true_hEv0->SetTitle("Efficiency vs true Ev");
  true_hv0->SetTitle("Efficiency vs true nu");
  true_hEm0->SetTitle("Efficiency vs true E_muon");
  true_hA0->SetTitle("Efficiency vs true LepNuAngle");

  true_hEv0->SetStats(0);
  true_hEv1->SetStats(0);
  true_hEv2->SetStats(0);
  true_hEv3->SetStats(0);
  true_hEm0->SetStats(0);
  true_hEm1->SetStats(0);
  true_hEm2->SetStats(0);
  true_hEm3->SetStats(0);
  true_hA0->SetStats(0);
  true_hA1->SetStats(0);
  true_hA2->SetStats(0);
  true_hA3->SetStats(0);

  true_hEv0->SetLineColor(kBlack);
  true_hEv1->SetLineColor(kRed);
  true_hEv2->SetLineColor(kGreen);
  true_hEv3->SetLineColor(kBlue);
  true_hEm0->SetLineColor(kBlack);
  true_hEm1->SetLineColor(kRed);
  true_hEm2->SetLineColor(kGreen);
  true_hEm3->SetLineColor(kBlue);
  true_hA0->SetLineColor(kBlack);
  true_hA1->SetLineColor(kRed);
  true_hA2->SetLineColor(kGreen);
  true_hA3->SetLineColor(kBlue);
 
  true_hEv0->SetLineWidth(2);
  true_hEv1->SetLineWidth(2);
  true_hEv2->SetLineWidth(2);
  true_hEv3->SetLineWidth(2);
  true_hEm0->SetLineWidth(2);
  true_hEm1->SetLineWidth(2);
  true_hEm2->SetLineWidth(2);
  true_hEm3->SetLineWidth(2);
  true_hA0->SetLineWidth(2);
  true_hA1->SetLineWidth(2);
  true_hA2->SetLineWidth(2);
  true_hA3->SetLineWidth(2);

  true_hEv0->Divide(hEv0);
  true_hEv1->Divide(hEv1);
  true_hEv2->Divide(hEv2);
  true_hEv3->Divide(hEv3);
  true_hEm0->Divide(hEm0);
  true_hEm1->Divide(hEm1);
  true_hEm2->Divide(hEm2);
  true_hEm3->Divide(hEm3);
  true_hA0->Divide(hA0);
  true_hA1->Divide(hA1);
  true_hA2->Divide(hA2);
  true_hA3->Divide(hA3);


  true_hEv0->GetXaxis()->SetTitle("Ev (GeV)");
  true_hEv0->GetYaxis()->SetTitle("Event Selection Efficiency");
  true_hEv0->SetMaximum(1.1);
  TCanvas *c = new TCanvas("c","",800,600);
  true_hEv0->Draw();
  true_hEv1->Draw("same");
  true_hEv2->Draw("same");
  true_hEv3->Draw("same");
  TLegend *lEv = new TLegend(0.60,0.10,0.90,0.30);
  lEv->AddEntry(true_hEv0,"no nu cut");
  lEv->AddEntry(true_hEv1,"true_nu < 2.0GeV");
  lEv->AddEntry(true_hEv2,"true_nu < 0.8GeV");
  lEv->AddEntry(true_hEv3,"true_nu < 0.3GeV");
  lEv->Draw();
  c->SaveAs("eff_Ev.png");

  true_hEm0->GetXaxis()->SetTitle("Em (GeV)");
  true_hEm0->GetYaxis()->SetTitle("Event Selection Efficiency");
  true_hEm0->SetMaximum(1.1);
  TCanvas *cEm = new TCanvas("cEm","",800,600);
  true_hEm0->Draw();
  true_hEm1->Draw("same");
  true_hEm2->Draw("same");
  true_hEm3->Draw("same");
  TLegend *lEm = new TLegend(0.60,0.10,0.90,0.30);
  lEm->AddEntry(true_hEm0,"no nu cut");
  lEm->AddEntry(true_hEm1,"true_nu < 2.0GeV");
  lEm->AddEntry(true_hEm2,"true_nu < 0.8GeV");
  lEm->AddEntry(true_hEm3,"true_nu < 0.3GeV");
  lEm->Draw();
  cEm->SaveAs("eff_Em.png");
  
  true_hA0->GetXaxis()->SetTitle("LepNuAngle");
  true_hA0->GetYaxis()->SetTitle("Event Selection Efficiency");
  true_hA0->SetMaximum(1.1);
  TCanvas *cA = new TCanvas("cA","",800,600);
  true_hA0->Draw();
  true_hA1->Draw("same");
  true_hA2->Draw("same");
  true_hA3->Draw("same");
  TLegend *lA = new TLegend(0.60,0.70,0.90,0.90);
  lA->AddEntry(true_hA0,"no nu cut");
  lA->AddEntry(true_hA1,"true_nu < 2.0GeV");
  lA->AddEntry(true_hA2,"true_nu < 0.8GeV");
  lA->AddEntry(true_hA3,"true_nu < 0.3GeV");
  lA->Draw();
  cA->SaveAs("eff_A.png");

  
  gStyle->SetNumberContours(999); gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

  hEvRecoVsEv0->SetStats(0);
  hEvRecoVsEv0->SetTitle("reco Ev vs true Ev (no nu cut)");
  hEvRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  hEvRecoVsEv0->GetYaxis()->SetTitle("reco Ev (GeV)");
  TCanvas *cEvRecoVsEv0 = new TCanvas("cEvRecoVsEv0","",800,600);
  cEvRecoVsEv0->SetLogz();
  hEvRecoVsEv0->Draw("colz");
  cEvRecoVsEv0->SaveAs("EvRecoVsEv0.png");
  
  hEvRecoVsEv3->SetStats(0);
  hEvRecoVsEv3->SetTitle("reco Ev vs true Ev (true_nu<0.3 GeV)");
  hEvRecoVsEv3->GetXaxis()->SetTitle("true Ev (GeV)");
  hEvRecoVsEv3->GetYaxis()->SetTitle("reco Ev (GeV)");
  TCanvas *cEvRecoVsEv3 = new TCanvas("cEvRecoVsEv3","",800,600);
  cEvRecoVsEv3->SetLogz();
  hEvRecoVsEv3->Draw("colz");
  cEvRecoVsEv3->SaveAs("EvRecoVsEv3.png");

  hvRecoVsv0->SetStats(0);
  hvRecoVsv0->GetXaxis()->SetTitle("true nu (GeV)");
  hvRecoVsv0->GetYaxis()->SetTitle("reco nu (GeV)");
  TCanvas *cvRecoVsv0 = new TCanvas("cvRecoVsv0","",800,600);
  cvRecoVsv0->SetLogz();
  hvRecoVsv0->Draw("colz");
  cvRecoVsv0->SaveAs("vRecoVsv0.png");


  hEvRes0->SetStats(0);
  hEvRes1->SetStats(0);
  hEvRes2->SetStats(0);
  hEvRes3->SetStats(0);

  hEvRes0->SetLineColor(kBlack);
  hEvRes1->SetLineColor(kRed);
  hEvRes2->SetLineColor(kGreen);
  hEvRes3->SetLineColor(kBlue);

  hEvRes0->SetLineWidth(2);
  hEvRes1->SetLineWidth(2);
  hEvRes2->SetLineWidth(2);
  hEvRes3->SetLineWidth(2);

  double norm = hEvRes0->Integral();
  hEvRes0->Scale(1./norm);
  hEvRes1->Scale(1./norm);
  hEvRes2->Scale(1./norm);
  hEvRes3->Scale(1./norm);

  hEvRes0->SetTitle("Feeddown Ev");
  hEvRes0->GetXaxis()->SetTitle("(E_{reco}-E_{true})/E_{true}");
  hEvRes0->GetYaxis()->SetTitle("entries");
  TCanvas *cRes = new TCanvas("cRes","",800,600);
  hEvRes0->Draw("hist");
  hEvRes1->Draw("hist same");
  hEvRes2->Draw("hist same");
  hEvRes3->Draw("hist same");
  TLegend *lRes = new TLegend(0.60,0.70,0.90,0.90);
  lRes->AddEntry(hEvRes0,"no nu cut");
  lRes->AddEntry(hEvRes1,"true_nu < 2.0GeV");
  lRes->AddEntry(hEvRes2,"true_nu < 0.8GeV");
  lRes->AddEntry(hEvRes3,"true_nu < 0.3GeV");
  lRes->Draw();
  cRes->SaveAs("EvRes.png"); 

  return(0);
}



