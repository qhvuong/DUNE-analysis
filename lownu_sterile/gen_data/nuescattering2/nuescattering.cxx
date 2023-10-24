#include <iostream>
#include <string>
#include <TF1.h>
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
#include <TParameter.h>

int main()
{

  // Load the CAF file = Common Analysis Format, a standard TTree that we use in DUNE
  TChain * tree = new TChain( "tree", "tree" );
  TChain * meta = new TChain( "meta", "meta" );

  for(int i = 40; i<50; i++){
  //if(i==39) continue;
  tree->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/marshalc/nue_study/FHC/nueFHC_%03d.root",i) );
  meta->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/marshalc/nue_study/FHC/nueFHC_%03d.root",i) );
  std::cout << "Nue File number:" << i << "\n";
  }
  
  double total_pot = 0.;
  double pot;
  meta->SetBranchAddress( "pot", &pot );
  const int Nfiles = meta->GetEntries();
  //const int Nfiles = 10000;
  for( int ii = 0; ii < Nfiles; ++ii ) {
    meta->GetEntry(ii);
    total_pot += pot;
  }
  double yrPOT = 1.1E21;
  //double scalePOT = yrPOT/total_pot;
  double scalePOT = 1.0;

  //std::cout << total_pot << "\t" << yrPOT*total_pot/1.1E21 << "\t" << scalePOT << "\n";

  const double mubins[20] = {0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,8.,12.,16.,20.,40.,100.};
  const double ebins[8] = {0.,2.,4.,6.,8.,10.,20.,100.};

  const Int_t nbinsX = 430;
  Double_t xEdges[nbinsX+1];
  xEdges[0]=0;
  for(int i=0; i<nbinsX+1; i++)
  {
    if(i<200)                xEdges[i+1] = xEdges[i] + 0.02;
    else if(i>=200 && i<240) xEdges[i+1] = xEdges[i] + 0.1;
    else if(i>=240 && i<400) xEdges[i+1] = xEdges[i] + 0.2;
    else if(i>=400 && i<420) xEdges[i+1] = xEdges[i] + 1.0;
    else                     xEdges[i+1] = xEdges[i] + 4.0;
  }
   
  TH2D *m_hElepRecoVsEv0   = new TH2D("m_hElepRecoVsEv0","",nbinsX,xEdges,50,0,16);
  TH2D *m_hElepRecoVsEv0_w = new TH2D("m_hElepRecoVsEv0_w","",nbinsX,xEdges,50,0,16);
  TH2D *e_hElepRecoVsEv0   = new TH2D("e_hElepRecoVsEv0","",nbinsX,xEdges,50,0,16);
  TH2D *e_hElepRecoVsEv0_w = new TH2D("e_hElepRecoVsEv0_w","",nbinsX,xEdges,50,0,16);

  TH2D *m_hElepRecoVsEv0_cov = new TH2D("m_hElepRecoVsEv0_cov","",19,mubins,50,0,16);
  TH2D *e_hElepRecoVsEv0_cov = new TH2D("e_hElepRecoVsEv0_cov","",7,ebins,50,0,16);

  TH1D *m_hElep0 = new TH1D("m_hElep0","",50,0,16);
  TH1D *e_hElep0 = new TH1D("e_hElep0","",50,0,16);
  TH1D *hElep0 = new TH1D("hElep0","",50,0,16);


  // information about the true neutrino interaction
  double vtx_x, vtx_y, vtx_z; // the position where the neutrino interaction occurred, in cm
  double Enu, E[2], p[2]; // the energy of the neutrino, in GeV
  double px[2], py[2], pz[2];
  double best_px[2], best_py[2], best_pz[2];
  int pdg[2];
  double BthetaX, BthetaY, Btheta;
  double PthetaX, PthetaY, Ptheta;
  double BthetaX_sm, BthetaY_sm, Btheta_sm;
  double wgt_MaCCQE[5];

  tree->SetBranchAddress( "Enu", &Enu );
  tree->SetBranchAddress( "E", &E );
  tree->SetBranchAddress( "pdg", &pdg );
  tree->SetBranchAddress( "px", &px );
  tree->SetBranchAddress( "py", &py );
  tree->SetBranchAddress( "pz", &pz );
  tree->SetBranchAddress( "best_px", &best_px );
  tree->SetBranchAddress( "best_py", &best_py );
  tree->SetBranchAddress( "best_pz", &best_pz );
  //tree->SetBranchAddress( "wgt_MaCCQE", &wgt_MaCCQE );

  TF1 *tsmear1 = new TF1( "tsmear1", "3.29 + 3.485*pow(x,-1.)", 0., 999.9 );
  TF1 *tsmear2 = new TF1( "tsmear2", "10.287 + 4.889*pow(x,-1.)", 0., 999.9 );
  TF1 *tsmearRatio = new TF1( "tsmearRatio", "0.039 + 0.551*pow(x,-1.) - 0.268*pow(x,-0.5)", 0., 999.9 );
  TF1 *doubleGaus = new TF1( "dg", "[0]*TMath::Exp(-0.5*pow(x/[1],2)) + [2]*TMath::Exp(-0.5*pow(x/[3],2))", -1000., 1000. );

  TRandom3 *rando = new TRandom3(8888);

  double s2tW = 0.23;
  double m_C_LL = -1./2. + s2tW; 
  double m_C_LR = s2tW;
  double e_C_LL = 1./2. + s2tW;
  double e_C_LR = s2tW;
  double y, sigma_m, sigma_e;

  double me = 510;  //keV
  double Ev_reco;

  const int N = tree->GetEntries();
  //const int N = 10000;
  //scalePOT = 10;

  double E_max;

  for( int ii = 0; ii < N; ++ii )
  {
    if( ii % 1000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
    tree->GetEntry(ii);

    if( pdg[1]==14 || pdg[1]==12 ){

    if(Enu > E_max) E_max = Enu;

    BthetaX = atan2(best_px[0],best_pz[0])*1E3;
    BthetaY = atan2(best_py[0],best_pz[0])*1E3;
    Btheta = atan2(sqrt(best_py[0]*best_py[0]+best_px[0]*best_px[0]),best_pz[0])*1E3;

    double p1 = tsmear1->Eval(E[0]);
    double p3 = tsmear2->Eval(E[0]);
    double ratio = tsmearRatio->Eval(E[0]);

    doubleGaus->SetParameter( 1, p1 );
    doubleGaus->SetParameter( 3, p3 );
    doubleGaus->SetParameter( 0, 1. );
    doubleGaus->SetParameter( 2, ratio );

    BthetaX_sm = BthetaX + doubleGaus->GetRandom();
    BthetaY_sm = BthetaY + doubleGaus->GetRandom();
    Btheta_sm = atan(sqrt(tan(BthetaX_sm/1E3)*tan(BthetaX_sm/1E3) + tan(BthetaY_sm/1E3)*tan(BthetaY_sm/1E3)))*1E3;    
      
    Ev_reco = E[0]/(1-E[0]*Btheta_sm*Btheta_sm/(2*me));

    for(int i=0; i<2; i++){
      p[i]=sqrt(px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i]);
    }
    y = E[0]/Enu;
    sigma_m = (m_C_LL*m_C_LL + m_C_LR*m_C_LR*(1.-y)*(1.-y));
    sigma_e = (e_C_LL*e_C_LL + e_C_LR*e_C_LR*(1.-y)*(1.-y));

  
    if(E[0]*Btheta_sm*Btheta_sm/1E3<3.){
      if(pdg[1] == 14){
        m_hElep0->Fill(E[0]); 

        //muon template weighted     
        m_hElepRecoVsEv0_w->Fill(Enu,E[0],sigma_e/sigma_m);

        //muon template     
        m_hElepRecoVsEv0->Fill(Enu,E[0]);
        m_hElepRecoVsEv0_cov->Fill(Enu,E[0]);

      }

      if(pdg[1] == 12){
        e_hElep0->Fill(E[0]); 
      
        //electron template weighted
        e_hElepRecoVsEv0_w->Fill(Enu,E[0],sigma_m/sigma_e);

        //electron template
        e_hElepRecoVsEv0->Fill(Enu,E[0]);
        e_hElepRecoVsEv0_cov->Fill(Enu,E[0]);

      }
    }
  }
  }


  m_hElepRecoVsEv0_w->Scale(scalePOT); 
  e_hElepRecoVsEv0_w->Scale(scalePOT);
  m_hElepRecoVsEv0->Scale(scalePOT);   
  e_hElepRecoVsEv0->Scale(scalePOT);

  m_hElepRecoVsEv0_cov->Scale(scalePOT);   
  e_hElepRecoVsEv0_cov->Scale(scalePOT);

  m_hElep0->Scale(scalePOT); 
  e_hElep0->Scale(scalePOT);

  hElep0->Add(m_hElep0); 
  hElep0->Add(e_hElep0);




  //gStyle->SetPalette(kColorPrintableOnGrey); TColor::InvertPalette();

/*
  m_hElepRecoVsEv0->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 3.0MeV)");
  m_hElepRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv0->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv0->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 3.0MeV)");
  e_hElepRecoVsEv0->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv0->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv = new TCanvas("cElepVsEv","",1800,1400);
  cElepVsEv->Divide(2,1);
  cElepVsEv->cd(1);
  e_hElepRecoVsEv0->Draw("colz");
  cElepVsEv->cd(2);
  m_hElepRecoVsEv0->Draw("colz");
  cElepVsEv->SaveAs("nue_Elep_templates.png");

  m_hElepRecoVsEv0_cov->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 3.0MeV)");
  m_hElepRecoVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv0_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv0_cov->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 3.0MeV)");
  e_hElepRecoVsEv0_cov->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv0_cov->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv_cov = new TCanvas("cElepVsEv_cov","",1800,1400);
  cElepVsEv_cov->Divide(2,2);
  cElepVsEv_cov->cd(1);
  e_hElepRecoVsEv0_cov->Draw("colz");
  cElepVsEv_cov->cd(2);
  e_hElepRecoVsEv2_cov->Draw("colz");
  cElepVsEv_cov->cd(3);
  m_hElepRecoVsEv0_cov->Draw("colz");
  cElepVsEv_cov->cd(4);
  m_hElepRecoVsEv2_cov->Draw("colz");
  cElepVsEv_cov->SaveAs("cov_nue_Elep.png");

  m_hElepRecoVsEv0_cov_sigma2->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 3.0MeV) (-1sigma)");
  m_hElepRecoVsEv0_cov_sigma2->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv0_cov_sigma2->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv0_cov_sigma2->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 3.0MeV) (-1sigma)");
  e_hElepRecoVsEv0_cov_sigma2->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv0_cov_sigma2->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepRecoVsEv2_cov_sigma2->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 0.5MeV) (-1sigma)");
  m_hElepRecoVsEv2_cov_sigma2->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv2_cov_sigma2->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv2_cov_sigma2->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 0.5MeV) (-1sigma)");
  e_hElepRecoVsEv2_cov_sigma2->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv2_cov_sigma2->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv_cov_sigma2 = new TCanvas("cElepVsEv_cov_sigma2","",1800,1400);
  cElepVsEv_cov_sigma2->Divide(2,2);
  cElepVsEv_cov_sigma2->cd(1);
  e_hElepRecoVsEv0_cov_sigma2->Draw("colz");
  cElepVsEv_cov_sigma2->cd(2);
  e_hElepRecoVsEv2_cov_sigma2->Draw("colz");
  cElepVsEv_cov_sigma2->cd(3);
  m_hElepRecoVsEv0_cov_sigma2->Draw("colz");
  cElepVsEv_cov_sigma2->cd(4);
  m_hElepRecoVsEv2_cov_sigma2->Draw("colz");
  cElepVsEv_cov_sigma2->SaveAs("cov_sigma2_nue_Elep.png");

  m_hElepRecoVsEv0_cov_sigma4->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 3.0MeV) (+1sigma)");
  m_hElepRecoVsEv0_cov_sigma4->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv0_cov_sigma4->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv0_cov_sigma4->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 3.0MeV) (+1sigma)");
  e_hElepRecoVsEv0_cov_sigma4->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv0_cov_sigma4->GetYaxis()->SetTitle("Elep_reco (GeV)");
  m_hElepRecoVsEv2_cov_sigma4->SetTitle("nu+e numu Elep_reco Vs true Ev (Etheta2 < 0.5MeV) (+1sigma)");
  m_hElepRecoVsEv2_cov_sigma4->GetXaxis()->SetTitle("true Ev (GeV)");
  m_hElepRecoVsEv2_cov_sigma4->GetYaxis()->SetTitle("Elep_reco (GeV)");
  e_hElepRecoVsEv2_cov_sigma4->SetTitle("nu+e nue Elep_reco Vs true Ev (Etheta2 < 0.5MeV) (+1sigma)");
  e_hElepRecoVsEv2_cov_sigma4->GetXaxis()->SetTitle("true Ev (GeV)");
  e_hElepRecoVsEv2_cov_sigma4->GetYaxis()->SetTitle("Elep_reco (GeV)");
  TCanvas *cElepVsEv_cov_sigma4 = new TCanvas("cElepVsEv_cov_sigma4","",1800,1400);
  cElepVsEv_cov_sigma4->Divide(2,2);
  cElepVsEv_cov_sigma4->cd(1);
  e_hElepRecoVsEv0_cov_sigma4->Draw("colz");
  cElepVsEv_cov_sigma4->cd(2);
  e_hElepRecoVsEv2_cov_sigma4->Draw("colz");
  cElepVsEv_cov_sigma4->cd(3);
  m_hElepRecoVsEv0_cov_sigma4->Draw("colz");
  cElepVsEv_cov_sigma4->cd(4);
  m_hElepRecoVsEv2_cov_sigma4->Draw("colz");
  cElepVsEv_cov_sigma4->SaveAs("cov_sigma4_nue_Elep.png");
*/


  TFile *out = new TFile("/dune/app/users/qvuong/data/lownu/nue_output_4.root","RECREATE");
  m_hElepRecoVsEv0->Write();
  m_hElepRecoVsEv0_w->Write();
  e_hElepRecoVsEv0->Write();
  e_hElepRecoVsEv0_w->Write();

  m_hElepRecoVsEv0_cov->Write();
  e_hElepRecoVsEv0_cov->Write();

  m_hElep0->Write();
  e_hElep0->Write();
  
  hElep0->Write();

  TParameter<double> totalPOT("total_pot", total_pot);
  totalPOT.Write();

  out->Close();

  return(0);
}

