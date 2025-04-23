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



/*
#include "DUNEStyle.h"


// save a lot of useless repetitive typing
TLegend * MakeLegend(float left=0.7, float bottom=0.5, float right=0.9, float top=0.85)
{
  auto leg = new TLegend(left, bottom, right, top);
  leg->SetFillStyle(0);  // unfortunately can't set this in TStyle :(

  return leg;
}


void plot1D(TCanvas *c, TH1D** hists, double a, double b, const char *name)
{
  c->Clear();
  c->cd();
  double nu_cuts[4] = {30., 2., 0.8, 0.3};
  TLegend * leg = MakeLegend(0.60, a, 0.90, b);
  TH1 * hFirst = nullptr;
  for (std::size_t histIdx = 0; histIdx < 4; histIdx++)
  {
    TH1D *h = hists[histIdx];
    auto color = dunestyle::colors::NextColor(dunestyle::colors::Cycle::OkabeIto, histIdx==0 ? 0 : -1);
    h->SetLineColor(color);
    h->SetFillStyle(0);
    dunestyle::CenterTitles(h);
    auto newh = h->DrawCopy(histIdx == 0 ? "" : "same");  // need to leak it so it doesn't disappear
    if (!hFirst)
      hFirst = newh;

    // we do this the hard way so the legend has the top-most histogram in the stack first
    leg->GetListOfPrimitives()->AddFirst(new TLegendEntry(newh, Form("true #nu < %.1f GeV", nu_cuts[histIdx]), "l"));
  }
  c->RedrawAxis();  // otherwise the last histogram drawn overlaps with the frame
  leg->Draw();
  hFirst->SetMaximum(hFirst->GetMaximum()*1.35); // make some space for the watermark
  dunestyle::Simulation();
  c->SaveAs(Form("%s.png",name));
}



void plot2D(TCanvas *c, TH2D* hists, double valMin, double valMax, const char *name)
{
  c->Clear();
  c->cd();
  dunestyle::CenterTitles(hists);
  hists->SetMinimum(valMin);
  hists->SetMaximum(valMax);
  hists->Draw("colz");
  dunestyle::Simulation();
  c->SaveAs(Form("%s.png",name));
}

*/


void xS_data()
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
  const int NEntries = meta->GetEntries();
  for( int ii = 0; ii < NEntries; ++ii ) {
    meta->GetEntry(ii);
    total_pot += pot;
  }
  double yrPOT = 1.1E21;
  double scalePOT = yrPOT/total_pot;
  std::cout << "scalePOT = " << scalePOT << "\n";

  TH2D *hmEvRecoVsEv[5], *heEvRecoVsEv[5];
  TH2D *hThetaVsEe[5], *hThetaVsEeReco[5];
  TH2D *hThetaVsEe_z[5], *hThetaVsEeReco_z[5];

  // No cuts
  TH1D *hmEv[5], *heEv[5];

  // With reconstruction cuts
  TH1D *hm[5], *he[5], *he_wc[5];

  TH1D *hmResE[5], *heResE[5], *heResE_wc[5];
  

  for (std::size_t histIdx = 0; histIdx < 5; histIdx++) {
    hmEv[histIdx] = new TH1D(Form("hmEv%zu", histIdx), "", 50, 0, 16);
    heEv[histIdx] = new TH1D(Form("heEv%zu", histIdx), "", 50, 0, 16);

    hm[histIdx] = new TH1D(Form("hm%zu", histIdx),"",50,0,16);
    he[histIdx] = new TH1D(Form("he%zu", histIdx),"",50,0,16);
    he_wc[histIdx] = new TH1D(Form("he_wc%zu", histIdx),"",50,0,16); 

    hmResE[histIdx] = new TH1D(Form("hmResE%zu",histIdx),"", 50, -1., 1.);
    heResE[histIdx] = new TH1D(Form("heResE%zu",histIdx),"", 50, -1., 1.);
    heResE_wc[histIdx] = new TH1D(Form("heResE_wc%zu",histIdx),"", 50, -1., 1.);

    hmEvRecoVsEv[histIdx] = new TH2D(Form("hmEvRecoVsEv%zu", histIdx), "", 50,0,16,50,0,16);
    heEvRecoVsEv[histIdx] = new TH2D(Form("heEvRecoVsEv%zu", histIdx), "", 50,0,16,50,0,16);

    hThetaVsEe[histIdx] = new TH2D(Form("hThetaVsEe%zu", histIdx), "", 50,0,16,100,0,1500);
    hThetaVsEeReco[histIdx] = new TH2D(Form("hThetaVsEeReco%zu", histIdx), "", 50,0,16,100,0,1500);

    hThetaVsEe_z[histIdx] = new TH2D(Form("hThetaVsEe_z%zu", histIdx), "", 50,0,16,50,0,150);
    hThetaVsEeReco_z[histIdx] = new TH2D(Form("hThetaVsEeReco_z%zu", histIdx), "", 50,0,16,50,0,150);
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

  double nu_t, nu_r;
  double nu_cuts[5] = {50., 5., 1., 0.5, 0.3};
  double LepE_sm, Ev_sm;

  TRandom3 *rando = new TRandom3(12345);

  const int N = tree->GetEntries();
  //const int N = 10000;
  for( int ii = 0; ii < N; ++ii )
  {
    if( ii % 100000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
    tree->GetEntry(ii);

     // Skip events that occur outside the "Fiducial Volume" which is a region in the middle of the detector
    // Basically we can't measure neutrinos that interact right next to the edge very well
    // These numbers are in cm; the detector goes from -357 to +357 in x, -150 to +150 in y, and 0 to 507 in z
    if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z > 350. ) continue;

    nu_t = eP + eN + ePip + ePim + ePi0 + eOther;
    LepE_sm = rando->Gaus(LepE,0.05);
    Ev_sm   = rando->Gaus(Ev,0.05);
    
      
    if(LepPDG == 13){
      double res = (Ev_reco-Ev)/Ev;

      for(int j=0; j<5; j++){
        if(nu_t < nu_cuts[j]) hmEv[j]->Fill(Ev); 
      }

      if( reco_numu && (muon_contained || muon_tracker || muon_ecal)){
        for(int j=0; j<5; j++){
          if(nu_t < nu_cuts[j]) {
            hm[j]->Fill(Ev);
            hmEvRecoVsEv[j]->Fill(Ev, LepE_sm);
            hmResE[j]->Fill(res);
            //std::cout << (Ev_sm-Ev)/Ev << "\n";
          }
        }    
      }
    }


    else if(LepPDG == 11){
      double res_sm = (Ev_sm-Ev)/Ev;
      double res = (Ev_reco-Ev)/Ev;
      for(int j=0; j<5; j++){
        if(nu_t < nu_cuts[j]) heEv[j]->Fill(Ev); 
      }
    

      // Without Etheta2 cut
      if( reco_nue ){
        //std::cout << LepNuAngle << "\t" << LepE << "\n";
        for(int j=0; j<5; j++){
          if(nu_t < nu_cuts[j]) {
            he[j]->Fill(Ev);

            heResE[j]->Fill(res);
            heResE_wc[j]->Fill(res_sm);

            heEvRecoVsEv[j]->Fill(Ev, LepE_sm);

            hThetaVsEe[j]->Fill(LepE, LepNuAngle*1E3);
            hThetaVsEe_z[j]->Fill(LepE, LepNuAngle*1E3);
            hThetaVsEeReco[j]->Fill(LepE_sm, LepNuAngle*1E3);
            hThetaVsEeReco_z[j]->Fill(LepE_sm, LepNuAngle*1E3);
          }
        }    
      }

      // With Etheta2 cut
      if( reco_nue && LepE_sm*LepNuAngle*LepNuAngle*1E3>3. ){
        for(int j=0; j<5; j++){
          if(nu_t < nu_cuts[j]){
            he_wc[j]->Fill(Ev);
            
          } 
        }    
      }
    }


  }

   
  for(int j=0; j<5; j++){
    hmEv[j]->Scale(scalePOT);
    heEv[j]->Scale(scalePOT);
    hm[j]->Scale(scalePOT);
    he[j]->Scale(scalePOT);
    he_wc[j]->Scale(scalePOT);

    hmResE[j]->Scale(scalePOT);
    heResE[j]->Scale(scalePOT);
    heResE_wc[j]->Scale(scalePOT);

    hmEvRecoVsEv[j]->Scale(scalePOT);
    heEvRecoVsEv[j]->Scale(scalePOT);

    hThetaVsEe[j]->Scale(scalePOT);
    hThetaVsEeReco[j]->Scale(scalePOT);
    hThetaVsEe_z[j]->Scale(scalePOT);
    hThetaVsEeReco_z[j]->Scale(scalePOT);
  }

  //hmResE[0]->Draw();
  //hmResE[4]->Draw("same");


  TFile *out = new TFile("outFile_forPlots.root","RECREATE");
  for(int j=0; j<5; j++){
    hmEv[j]->Write();
    heEv[j]->Write();
    hm[j]->Write();
    he[j]->Write();
    he_wc[j]->Write();

    hmResE[j]->Write();
    heResE[j]->Write();
    heResE_wc[j]->Write();

    hmEvRecoVsEv[j]->Write();
    heEvRecoVsEv[j]->Write();

    hThetaVsEe[j]->Write();
    hThetaVsEeReco[j]->Write();
    hThetaVsEe_z[j]->Write();
    hThetaVsEeReco_z[j]->Write();
  }

  out->Close();

}




