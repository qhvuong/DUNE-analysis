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

static const char data_path[] = "/exp/dune/app/users/qvuong/data/lownu/uncertainties/det_covmtr";

void plot1D(TH1D* hp, TH1D *hn, TH1D *ho, const char* name){
    ho->SetLineColor(kBlack);
    hp->SetLineColor(kBlue);
    hn->SetLineColor(kRed);
    ho->SetLineWidth(2.);
    hp->SetLineWidth(2.);
    hn->SetLineWidth(2.);
    ho->SetStats(0);
    hp->SetStats(0);
    hn->SetStats(0);
    TCanvas *c = new TCanvas("c", "", 800, 600);
    c->SetLogy();
    c->SetGrid();
    hp->Draw();
    hn->Draw("same");
    ho->Draw("same");
    TLegend *lg = new TLegend(0.65,0.75,0.9,0.9);
    lg->AddEntry(ho, "nominal");
    lg->AddEntry(hp, "positive");
    lg->AddEntry(hn, "negative");
    lg->Draw("same");
    c->SaveAs(Form("%s.png",name));
}


int main()
{
// CC and nu+e ROOT files ---------------------------------------------------------------------------------------------
    TChain * tree = new TChain( "tree", "tree" );
    TChain * meta = new TChain( "meta", "meta" );

    for(int i = 0; i<69; i++){
        tree->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/marshalc/nue_study/FHC/nueFHC_%03d.root",i) );
        meta->Add( Form("root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/marshalc/nue_study/FHC/nueFHC_%03d.root",i) );
    }

//protons on target for CC and nu+e------------------------------------------------------------------------------------
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

//nu+e branches -------------------------------------------------------------------------------------------------------
    double vtx_x, vtx_y, vtx_z;
    int pdg[2]; 
    double best_px[2], best_py[2], best_pz[2];
    double E[2];

    tree->SetBranchAddress( "vtx_x", &vtx_x );
    tree->SetBranchAddress( "vtx_y", &vtx_y );
    tree->SetBranchAddress( "vtx_z", &vtx_z );
    tree->SetBranchAddress( "pdg", &pdg );
    tree->SetBranchAddress( "best_px", &best_px );
    tree->SetBranchAddress( "best_py", &best_py );
    tree->SetBranchAddress( "best_pz", &best_pz );
    tree->SetBranchAddress( "E", &E );
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
    TH1D* hn = new TH1D("hn", "nominal nu+e", nbinsY, yEdges);
    TH1D* hnp = new TH1D("hnp", "+2.5% nu+e", nbinsY, yEdges);
    TH1D* hnn = new TH1D("hnn", "-2.5% nu+e", nbinsY, yEdges);

//---------------------------------------------------------------------------------------------------------------------

    //for e samples
    double psh_e = 1.025;
    double nsh_e = 0.975;

//from your script for smearing electron energies in nu+e sample
double BthetaX, BthetaY, BthetaX_sm, BthetaY_sm, Btheta_sm;

//fill nu+e sample ----------------------------------------------------------------------------------------------------
    const int N = tree->GetEntries();
    //const int N = 10000;
    for(int ii = 0; ii < N; ++ii){
        if( ii % 1000 == 0 ) printf( "%.2f percent of %d Events...\n", ii*100.0/N, N );
        tree->GetEntry(ii);
        //FV cut
        if( abs(vtx_x) > 300. || abs(vtx_y) > 100. || vtx_z < 50. || vtx_z > 350. ) continue;
            
        //find numu and nue
        if( pdg[1]==14 || pdg[1]==12 ){
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
          
            //double Ev_reco = E[0]/(1-E[0]*Btheta_sm*Btheta_sm/(2*me));

            // Filling
            double e_Ep = psh_e * E[0];
            double e_En = nsh_e * E[0];
            if( E[0]*Btheta_sm*Btheta_sm/1E3<=3. ) hn->Fill(E[0]);
            if( e_Ep*Btheta_sm*Btheta_sm/1E3<=3. ) hnp->Fill(e_Ep);
            if( e_En*Btheta_sm*Btheta_sm/1E3<=3. ) hnn->Fill(e_En);

            // Energy resolution scaling 10%
            //double Eresp = E[0] + (1.+0.1)*(E[0]-LepE);
            //double Eresn = LepE_sm + (1.-0.1)*(LepE_sm-LepE);
        }
    }

//---------------------------------------------------------------------------------------------------------------------

//save to new ROOT file -----------------------------------------------------------------------------------------------
    TFile *out = new TFile(Form("%s/output_nue.root",data_path),"RECREATE");
    hn->Write();  //unshifted energies
    hnp->Write();  //positive nue shift
    hnn->Write();  //negative nue shift


    //total POT is saved for CC and nu+e samples to the file
    TParameter<double> totalPOT("total_pot", total_pot);
    totalPOT.Write();
    out->Close();

//---------------------------------------------------------------------------------------------------------------------


    plot1D(hnp, hnn, hn, "nueE");

    return(0);
}