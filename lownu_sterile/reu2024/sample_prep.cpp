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
#include <cmath>

int sample_prep()
{
// CC and nu+e ROOT files ---------------------------------------------------------------------------------------------
    TChain * treeCC = new TChain( "cafTree", "cafTree" );
    TChain * metaCC = new TChain( "meta", "meta" );
    TChain * treeNue = new TChain( "tree", "tree" );
    TChain * metaNue = new TChain( "meta", "meta" );

    for(int i = 50; i<60;i++){
        treeNue->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/marshalc/nue_study/FHC/nueFHC_%03d.root",i) );
        metaNue->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/marshalc/nue_study/FHC/nueFHC_%03d.root",i) );
        treeCC->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%02d.root", i) );
        metaCC->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/LBL_TDR/CAFs/v4/ND_FHC_FV_%02d.root", i) );
    }

//protons on target for CC and nu+e------------------------------------------------------------------------------------
    double totPOTCC = 0.;
    double potCC;
    metaCC->SetBranchAddress("pot", &potCC);
    const int NfilesCC = metaCC->GetEntries();
    for( int ii = 0; ii < NfilesCC; ++ii ) {
        metaCC->GetEntry(ii);
        totPOTCC += potCC;
    }

    double totPOTNue = 0.;
    double potNue;
    metaNue->SetBranchAddress("pot", &potNue);
    const int NfilesNue = metaNue->GetEntries();
    for( int ii = 0; ii < NfilesNue; ++ii ) {
        metaNue->GetEntry(ii);
        totPOTNue += potNue;
    }

//CC branches ---------------------------------------------------------------------------------------------------------
    double vtx_xCC, vtx_yCC, vtx_zCC;
    double Ev_reco, Elep_reco;
    int reco_numu, reco_nue;
    int muon_contained, muon_tracker;
    double LepE, Ev;

    treeCC->SetBranchAddress( "vtx_x", &vtx_xCC );
    treeCC->SetBranchAddress( "vtx_y", &vtx_yCC );
    treeCC->SetBranchAddress( "vtx_z", &vtx_zCC );
    treeCC->SetBranchAddress( "Ev_reco", &Ev_reco );
    treeCC->SetBranchAddress( "Elep_reco", &Elep_reco );
    treeCC->SetBranchAddress( "reco_numu", &reco_numu );
    treeCC->SetBranchAddress( "reco_nue", &reco_nue );
    treeCC->SetBranchAddress( "muon_contained", &muon_contained );
    treeCC->SetBranchAddress( "muon_tracker", &muon_tracker );
    treeCC->SetBranchAddress( "LepE", &LepE );
    treeCC->SetBranchAddress( "Ev", &Ev );
//---------------------------------------------------------------------------------------------------------------------

//nu+e branches -------------------------------------------------------------------------------------------------------
    double vtx_xNue, vtx_yNue, vtx_zNue;
    int pdg[2]; 
    double best_px[2], best_py[2], best_pz[2];
    double E[2];

    treeNue->SetBranchAddress( "vtx_x", &vtx_xNue );
    treeNue->SetBranchAddress( "vtx_y", &vtx_yNue );
    treeNue->SetBranchAddress( "vtx_z", &vtx_zNue );
    treeNue->SetBranchAddress( "pdg", &pdg );
    treeNue->SetBranchAddress( "best_px", &best_px );
    treeNue->SetBranchAddress( "best_py", &best_py );
    treeNue->SetBranchAddress( "best_pz", &best_pz );
    treeNue->SetBranchAddress( "E", &E );
//---------------------------------------------------------------------------------------------------------------------

//CC binning and smear setup. small hist ------------------------------------------------------------------------------
    const Int_t nbinsCC = 56;
    double CCEdges[nbinsCC+1];
    CCEdges[0] = 0.;
    CCEdges[1] = 0.3;

    for(int i=1; i<nbinsCC+1; i++)
    {
      if(i<38)               CCEdges[i+1] = CCEdges[i] + 0.1;
      else if(i>=38 && i<43) CCEdges[i+1] = CCEdges[i] + 0.2;
      else if(i>=43 && i<48) CCEdges[i+1] = CCEdges[i] + 0.4;
      else if(i>=48 && i<53) CCEdges[i+1] = CCEdges[i] + 0.8;
      else if(i>=53 && i<55) CCEdges[i+1] = CCEdges[i] + 1.5;
      else                   CCEdges[i+1] = CCEdges[i] + 2.0;
    }

    TRandom *r = new TRandom(12345);

    //I use "o" to denote nominal/original energies, "p" for positive, "n" for negative
    // "c" is for contained and "t" is for tracker
    TH1D* hcmo = new TH1D("hcmo", "nominal contained muCC", nbinsCC, CCEdges); 
    TH1D* hcmp = new TH1D("hcmp", "+2% contained muCC", nbinsCC, CCEdges);
    TH1D* hcmn = new TH1D("hcmn", "-2% contained muCC", nbinsCC, CCEdges);

    TH1D* htmo = new TH1D("htmo", "nominal tracker muCC", nbinsCC, CCEdges);
    TH1D* htmp = new TH1D("htmp", "+2% tracker muCC", nbinsCC, CCEdges);
    TH1D* htmn = new TH1D("htmn", "-2% tracker muCC", nbinsCC, CCEdges);

    //I will add contained and tracker to these histograms. For the sake of completion I fill the type of mu separately
    //combined mu histograms are denoted with "mu"
    TH1D* hmo = new TH1D("hmo", "nominal muCC", nbinsCC, CCEdges);
    TH1D* hmp = new TH1D("hmp", "+2% muCC", nbinsCC, CCEdges);
    TH1D* hmn = new TH1D("hmn", "-2% muCC", nbinsCC, CCEdges);

    TH1D* heo = new TH1D("heo", "nominal eCC", nbinsCC, CCEdges);
    TH1D* hep = new TH1D("hep", "+2.5% eCC", nbinsCC, CCEdges);
    TH1D* hen = new TH1D("hen", "-2.5% eCC", nbinsCC, CCEdges);
//---------------------------------------------------------------------------------------------------------------------

//nu+e binning and smear setup for smaller histograms (nu+e sample only) ----------------------------------------------
    const double ebins[8] = {0.,2.,4.,6.,8.,10.,20.,100.};
    const Int_t nbinsY = 8;
    const double yEdges[9] = {0., 0.3, 0.6, 0.92, 1.3, 1.75, 2.45, 3.9, 16.0};

    TF1 *tsmear1 = new TF1( "tsmear1", "3.29 + 3.485*pow(x,-1.)", 0., 999.9 );
    TF1 *tsmear2 = new TF1( "tsmear2", "10.287 + 4.889*pow(x,-1.)", 0., 999.9 );
    TF1 *tsmearRatio = new TF1( "tsmearRatio", "0.039 + 0.551*pow(x,-1.) - 0.268*pow(x,-0.5)", 0., 999.9 );
    TF1 *doubleGaus = new TF1( "dg", "[0]*TMath::Exp(-0.5*pow(x/[1],2)) + [2]*TMath::Exp(-0.5*pow(x/[3],2))", -1000., 1000. );

    TRandom3 *rando = new TRandom3(12345);
    double me = 510;

    //same naming scheme. "o" for nominal, "p" for positive, "n" for negative
    //for clarity, maybe adding "e" to denote electron would be good but these are added to a larger electron histogram
    TH1D* hno = new TH1D("hno", "nominal nu+e", nbinsY, yEdges);
    TH1D* hnp = new TH1D("hnp", "+2.5% nu+e", nbinsY, yEdges);
    TH1D* hnn = new TH1D("hnn", "-2.5% nu+e", nbinsY, yEdges);
//---------------------------------------------------------------------------------------------------------------------

//define large histograms (all three samples go here) -----------------------------------------------------------------
    const int nbins = 2*nbinsCC + nbinsY;
    double Edges[nbins + 1];
    //bin edges all end on 16, so I shift eCC sample by 16 and nu+e by 32
    //this should then place muCC first, eCC second, and nu+e last in terms of bin numbers with correct bin edges
    for(int i = 0; i < nbins + 1; i++){
        if(i < nbinsCC)             Edges[i] = CCEdges[i];  //for muCC
        else if(i < 2*nbinsCC)      Edges[i] = 16. + CCEdges[i - nbinsCC];  //for eCC
        else                        Edges[i] = 32. + yEdges[i - 2*nbinsCC];  //for nu+e
    }

    TH1D* hnom = new TH1D("hnom", "nominal sample", nbins, Edges);  //nominal energy for all

    TH1D* hpmu = new TH1D("hpmu", "+2% mu sample", nbins, Edges);  //positive mu reco energy shift
    TH1D* hnmu = new TH1D("hnmu", "-2% mu sample", nbins, Edges);  //negative mu reco energy shift

    TH1D* hpel = new TH1D("hpel", "+2.5% electron sample", nbins, Edges);  //positive electron energy shift (CC and nu+e)
    TH1D* hnel = new TH1D("hnel", "-2.5% electron sample", nbins, Edges);  //negative electron energy shift (CC and nu+e)
//---------------------------------------------------------------------------------------------------------------------

    //for mu sample
    double psh_m = 1.02; 
    double nsh_m = 0.98;

    //for e samples
    double psh_e = 1.025;
    double nsh_e = 0.975;

//fill CC samples (mu and eCC) ----------------------------------------------------------------------------------------
    const int NCC = treeCC->GetEntries();
    for(int ii = 0; ii < NCC; ++ii){
        treeCC->GetEntry(ii);
        //FV cut is first
        if( abs(vtx_xCC) <= 300. && abs(vtx_yCC) <= 100. && vtx_zCC >= 50. && vtx_zCC <= 350. ){
            //check for numu
            if(reco_numu == 1){
                double m_nu = Ev_reco - Elep_reco;
                //low nu cut for mu 
                if(m_nu < 0.3){
                    if(muon_contained == 1){
                        hcmo->Fill(Elep_reco);
                        hcmp->Fill(psh_m * Elep_reco);
                        hcmn->Fill(nsh_m * Elep_reco);
                    }
                    else if(muon_tracker == 1){
                        htmo->Fill(Elep_reco);
                        htmp->Fill(psh_m * Elep_reco);
                        htmn->Fill(nsh_m * Elep_reco);
                    }
                }
            }
            //now do nue 
            else if(reco_nue == 1){
                //these steps come from the script you sent me to smear electron energy
                double LepE_sm = r->Gaus(LepE, 0.05);
                double Ev_sm = r->Gaus(Ev, 0.05);
                double e_nu = Ev_sm - LepE_sm;
                //low nu cut for e
                if(e_nu < 0.3){
                    heo->Fill(LepE_sm);
                    hep->Fill(psh_e * LepE_sm);
                    hen->Fill(nsh_e * LepE_sm);
                }
            }
        }
    }

    //I add contained and tracker mu energies to the histograms that hold all mu energies
    //alternatively get rid of the if statements with tracker and contained and just fill these histograms directly
    hmo->Add(hcmo);
    hmo->Add(htmo);

    hmp->Add(hcmp);
    hmp->Add(htmp);

    hmn->Add(hcmn);
    hmn->Add(htmn);

//---------------------------------------------------------------------------------------------------------------------

    //from your script for smearing electron energies in nu+e sample
    double BthetaX, BthetaY, BthetaX_sm, BthetaY_sm, Btheta_sm;

//fill nu+e sample ----------------------------------------------------------------------------------------------------
    const int Nnue = treeNue->GetEntries();
    for(int ii = 0; ii < Nnue; ++ii){
        treeNue->GetEntry(ii);
        //FV cut
        if( abs(vtx_xNue) <= 300. && abs(vtx_yNue) <= 100. && vtx_zNue >= 50. && vtx_zNue <= 350. ){
            //find numu and nue
            if( pdg[1]==14 || pdg[1]==12 ){
                //from your script
                BthetaX = atan2(best_px[0],best_pz[0])*1E3;
                BthetaY = atan2(best_py[0],best_pz[0])*1E3;

                double p1 = tsmear1->Eval(E[0]);
                double p3 = tsmear2->Eval(E[0]);
                double ratio = tsmearRatio->Eval(E[0]);

                doubleGaus->SetParameter( 1, p1 );
                doubleGaus->SetParameter( 3, p3 );
                doubleGaus->SetParameter( 0, 1. );
                doubleGaus->SetParameter( 2, ratio );

                BthetaX_sm = BthetaX + doubleGaus->GetRandom();  //smears the true
                BthetaY_sm = BthetaY + doubleGaus->GetRandom();
                Btheta_sm = atan(sqrt(tan(BthetaX_sm/1E3)*tan(BthetaX_sm/1E3) + tan(BthetaY_sm/1E3)*tan(BthetaY_sm/1E3)))*1E3;    
          
                double Ev_reco = E[0]/(1-E[0]*Btheta_sm*Btheta_sm/(2*me));

                //low nu cut
                if( E[0]*Btheta_sm*Btheta_sm/1E3<3. ){
                    hno->Fill(E[0]);
                    hnp->Fill(psh_e * E[0]);
                    hnn->Fill(nsh_e * E[0]);
                } 
            }
        }
    }
//---------------------------------------------------------------------------------------------------------------------

//fill large histograms -----------------------------------------------------------------------------------------------
    for(int ii = 1; ii < nbins + 1; ++ii){
        //bins 1 through 56 are for muCC
        if(ii < nbinsCC + 1){
            hnom->SetBinContent(ii, hmo->GetBinContent(ii));

            hpmu->SetBinContent(ii, hmp->GetBinContent(ii));
            hnmu->SetBinContent(ii, hmn->GetBinContent(ii));

            hpel->SetBinContent(ii, hmo->GetBinContent(ii));
            hnel->SetBinContent(ii, hmo->GetBinContent(ii));
        }
        //bins 57 through 112 are for eCC
        else if(ii < 2*nbinsCC + 1){
            hnom->SetBinContent(ii, heo->GetBinContent(ii - nbinsCC));

            hpmu->SetBinContent(ii, heo->GetBinContent(ii - nbinsCC));
            hnmu->SetBinContent(ii, heo->GetBinContent(ii - nbinsCC));

            hpel->SetBinContent(ii, hep->GetBinContent(ii - nbinsCC));
            hnel->SetBinContent(ii, hen->GetBinContent(ii - nbinsCC));            
        }
        //bins 113 through 120 are for e in nu+e scattering
        else{
            hnom->SetBinContent(ii, hno->GetBinContent(ii - 2*nbinsCC));

            hpmu->SetBinContent(ii, hno->GetBinContent(ii - 2*nbinsCC));
            hnmu->SetBinContent(ii, hno->GetBinContent(ii - 2*nbinsCC));

            hpel->SetBinContent(ii, hnp->GetBinContent(ii - 2*nbinsCC));
            hnel->SetBinContent(ii, hnn->GetBinContent(ii - 2*nbinsCC));
        }
    }

//save to new ROOT file -----------------------------------------------------------------------------------------------
    TFile *out = new TFile("./CCnue_output_1.root","RECREATE");
    hnom->Write();  //unshifted energies
    hpmu->Write();  //positive mu shift
    hnmu->Write();  //negative mu shift
    hpel->Write();  //positive electron shift
    hnel->Write();  //negative electron shift

    //total POT is saved for CC and nu+e samples to the file
    TParameter<double> total_POTCC("totPOTCC", totPOTCC);
    total_POTCC.Write();

    TParameter<double> total_POTNue("totPOTNue", totPOTNue);
    total_POTNue.Write();

    out->Close();

//---------------------------------------------------------------------------------------------------------------------
    return(0);
}